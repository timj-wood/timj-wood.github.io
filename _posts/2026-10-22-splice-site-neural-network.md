---
title: "Creating a (mini) neural network to identify splice regions"
category: machine-learning
tags: [neural-networks, genomics]
---

And so I return to the lands of familiarity (bioinformatics). For my next post, I wanted to develop something that consists of my favourite [buzzwords](https://en.wikipedia.org/wiki/Buzzword), **machine learning** and **genomics**, which leads to the *holy grail*: ***neural networks***. I hadn't touched on machine learning before this year, but the more I read into it, the more fascinating the science behind it turns out to be, so developing (_some_) sort of neural network skillset has become a priority.

Fortunately, I will not be going in depth with the mathematical theory for neural networks. If you want more theory, I'd really recommend 3Blue1Brown's incredible series:

<iframe align="center" width="560" height="315" src="https://www.youtube.com/embed/aircAruvnKk" frameborder="0" allowfullscreen></iframe>

The code follows the Jean Golding Institute's [Introduction to AI](https://bristol-training.github.io/intro-to-ai/) course, with their [Data Analysis in Python (Part 2)](https://bristol-training.github.io/data-analysis-python-2/) course for the baseline. All figures and values in this post come from the repo [splice-site-nn](https://github.com/timj-wood/splice-site-nn), so be sure to have a look to get a better idea of the methods used.

## What is splicing?

When a gene is transcribed, the result is a pre-mRNA that still contains the introns (the non-coding regions). These have to be cut out before the message can move on to the next step, translation. So how do your cells know which parts of the transcript are introns, you might ask? Well, (almost) every intron begins with the two bases "GT" and ends with "AG". It's very close to a universal rule, which makes it a good target for my first neural network: if the model has learnt what it should, I can spot it immediately, and if it has learnt something else then I'll also be able to tell (because it will be wrong).

So the question I'll try to answer with my model is: given a window of 60 bases, is there a splice junction in the middle, and if so which kind? Along the way, I can find out whether the model rediscovers the GT…AG pattern all by itself. I'll compare the results against a logistic-regression baseline, using classical machine-learning methods first.

## The data

Conveniently, there exists a 35-year-old dataset that suits my needs nicely. The UCI [splice-junction dataset](https://doi.org/10.24432/C5M888) has 3,190 windows of 60 bases from primate genes, each labelled `EI` (exon to intron, a donor site), `IE` (intron to exon, an acceptor site) or `N` (neither). This data was originally developed for a paper using neural networks to recognise DNA sequences ([Noordewier *et al.*, 1991](https://proceedings.neurips.cc/paper/1990/file/8efb100a295c0c690931222ff4467bb8-Paper.pdf)), which is rather fitting.

However, there is a pressing issue with the data. A typical model cannot read a DNA sequence, e.g. `CCAGTCAT`. The standard fix for this is [**one-hot encoding**](https://www.sciencedirect.com/topics/computer-science/one-hot-encoding), where each position becomes four columns (one per base), with the base that is present assigned a value of 1 and the other three a 0. So one sequence becomes 240 numbers, and "position +1 is a G" is something tangible which a model can put a weight on.

Certain positions (60 out of 191,400) use [IUPAC ambiguity codes](https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html). For example, a `D` indicates "A, G, or T". Instead of discarding this data, you can spread the 1 evenly across the possible bases, so `D` becomes $\frac{1}{3}$, $\frac{1}{3}$, $\frac{1}{3}$ for A, G and T. This felt necessary to ensure that no data was excluded.

```python
X = one_hot(df["sequence"])   # (3190, 240) array of 0s, 1s and the odd 1/3
y = df["label"].values        # "EI", "IE" or "N"
```

## Creating a model baseline

Before we get to the fun neural-related stuff, it is good practice to fit some simpler classifiers first. We hold back a quarter of the data (798 sequences) to act as a test set which the models won't see in training.

```python
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42, stratify=y
)
```

```python
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report

logreg = LogisticRegression(max_iter=2000)
logreg.fit(X_train, y_train)

pred = logreg.predict(X_test)
acc_logreg = accuracy_score(y_test, pred)
print(f"Logistic regression test accuracy: {acc_logreg:.3f}")
```

**Logistic regression** got **93.7%** of the test sequences correct. For comparison, $k$-nearest neighbours, with $k$ chosen by five-fold cross-validation, scored 85.1%:

```python
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GridSearchCV

hyperparameters = {"n_neighbors": [1, 3, 5, 7, 11, 15, 21, 31]}
model = GridSearchCV(KNeighborsClassifier(), hyperparameters)
model.fit(X_train, y_train)

print(model.best_estimator_)
print(f"kNN test accuracy: {model.score(X_test, y_test):.3f}")
```

This disparity seems in order once you account for the shape of one-hot vectors. Two random 60-base sequences will rarely come close to one another in a 240-dimensional space, so asking the neighbours isn't much help.

Interestingly, logistic regression produces one weight per (position, base) for each class, meaning that the weights can be laid out as a 4 × 60 heatmap:

```python
W = logreg.coef_[list(logreg.classes_).index("EI")]   # 240 weights for the EI class
grid = W.reshape(60, 4).T                               # 4 bases × 60 positions
```

![Logistic regression weights for the EI and IE classes by position and base](/images/splice/logreg_weights.png)

For `EI`, the largest weights are T at +2 and G at +1. For `IE`, they're G at −1 and A at −2, with a stretch of C and T preferred further upstream (this would be the [polypyrimidine tract](https://en.wikipedia.org/wiki/Polypyrimidine_tract) that sits before acceptor sites). In sooth, the *humble* linear model was able to identify the GT…AG splicing pattern. Keep this result in mind, because we don't want anything below this baseline for our neural network.

## Building a neuron from scratch

A single artificial neuron does two things. First, it **aggregates**: a weighted sum of its inputs plus a bias. Then it **activates**: it pushes that sum through a non-linear function.

$$
z = \sum_j w_j x_j + b,
\qquad a = \sigma(z) = \frac{1}{1 + e^{-z}}.
$$

The sigmoid $\sigma$ compresses any number into the range 0 to 1, so the output can be interpreted as a probability. Keeping things simple, I started with a yes or no version of this scenario, basically asking "is there *any* splice junction in this section of bases?" (`EI` or `IE` gives 1, whereas `N` gives 0).

To train it, we need a **loss**, which is essentially a number that states how wrong the prediction was. For a probability the natural choice is binary cross-entropy,

$$
\ell = -[y \log a + (1-y)\log(1-a)],
$$

which punishes a confident wrong answer much more harshly than a hesitant one.

Then we need to know how the loss changes when a weight is nudged. The weight $w_j$ affects $z$, which affects $a$, which affects $\ell$, so by the chain rule we multiply the three sensitivities together. Each one is easy: $\partial z/\partial w_j = x_j$, $\partial a/\partial z = a(1-a)$ for the sigmoid, and $\partial \ell/\partial a = (a-y)/(a(1-a))$ for cross-entropy. The last two cancel almost completely, leaving

$$
\frac{\partial \ell}{\partial w_j} = (a - y)\, x_j, \qquad
\frac{\partial \ell}{\partial b} = a - y.
$$

So the gradient for each weight comes down to "how wrong was I?" multiplied by "was that feature switched on?". If the neuron was confidently wrong about a sequence that had a G at +1, the weight for G at +1 gets a big push; if the feature wasn't there, then nothing happens to it. The cancellation is not by chance. Sigmoid and cross-entropy complement each other nicely, and are naturally paired together. Average that over all 2,392 training sequences, write it as a matrix product, and the entire training loop is just a few lines of code:

```python
def sigmoid(z):
    return 1 / (1 + torch.exp(-z))

w = torch.zeros(240, 1)   # one weight per (position, base)
b = torch.zeros(1)
lr, n = 0.5, len(X_train)

for epoch in range(400):
    z = X_train @ w + b            # aggregate
    a = sigmoid(z)                 # activate
    grad_z = (a - y_train) / n     # how wrong, on average
    w -= lr * (X_train.T @ grad_z) # the gradient step
    b -= lr * grad_z.sum()
```

After 400 steps of plain gradient descent, the neuron classified 93.4% of the test set correctly.

Subsequently, I built the same neuron with PyTorch, using `nn.Linear(240, 1)`, `BCEWithLogitsLoss` and `loss.backward()`: same learning rate, same number of steps and same zero initialisation.

```python
neuron = nn.Linear(240, 1)
nn.init.zeros_(neuron.weight); nn.init.zeros_(neuron.bias)
loss_fn = nn.BCEWithLogitsLoss()               # sigmoid + cross-entropy in one
optimiser = torch.optim.SGD(neuron.parameters(), lr=0.5)

for epoch in range(400):
    optimiser.zero_grad()
    loss = loss_fn(neuron(X_train), y_train)
    loss.backward()                            # PyTorch works out the gradients
    optimiser.step()                           # w -= lr * grad
```

![Training loss for the hand-written neuron and the PyTorch one, overlaid](/images/splice/neuron_two_ways.png)

The two loss curves sit exactly on top of each other, and the final weights agree to within $2 \times 10^{-7}$. What I realised doing this is that there's no magic in `loss.backward()`: it's the chain rule applied automatically, and backpropagation is the same thing done layer by layer, with the intermediate values stored to prevent repeat calculations.

Whilst this is an insightful exercise, PyTorch does everything for you... so once you've done it by hand once, just use PyTorch.

## What the neuron learned

A single neuron has one weight per (position, base), so its 240 weights are quite literally a picture of what it thinks a splice junction looks like:

![The neuron's weights by position and base; red is evidence for a junction, blue against](/images/splice/neuron_weights.png)

The four largest weights, in order: G at −1, A at −2, T at +2, G at +1. That's AG on the left of the junction and GT on the right. Because this neuron was asked whether there's a junction of *either* kind, it has found both signals and stacked them. Nobody told it about introns. It was handed 240 columns of zeros and ones and a label, and it worked out the consensus sequence that took molecular biologists years to establish. I knew this was the result I was hoping for, and it still landed harder than I expected when the plot came up.

## Adding a hidden layer

One neuron can only add up independent pieces of evidence. A hidden layer lets the network combine features ("G at +1 *and* T at +2 *and* a pyrimidine-rich stretch before it"), which a plain sum can't express. So I built a small multi-layer perceptron: 240 inputs, 32 hidden units with ReLU activation ($f(z) = \max(0, z)$), and 3 outputs for `EI`, `IE` and `N`, trained with mini-batch gradient descent for 40 epochs.

```python
mlp = nn.Sequential(
    nn.Linear(240, 32),
    nn.ReLU(),
    nn.Linear(32, 3),
)
loss_fn = nn.CrossEntropyLoss()
optimiser = torch.optim.SGD(mlp.parameters(), lr=0.1)

for epoch in range(40):
    for idx in torch.randperm(len(X_train)).split(64):   # shuffled mini-batches
        optimiser.zero_grad()
        loss = loss_fn(mlp(X_train[idx]), y_train_idx[idx])
        loss.backward()
        optimiser.step()
```

Two practical details are worth mentioning. Weights can't all start at zero, otherwise every hidden unit would compute the same thing, receive the same gradient and stay identical forever (and with ReLU they'd all be stuck outputting zero), so PyTorch's random initialisation does that job. And instead of one gradient step over all 2,392 sequences per epoch, the network takes a step after every 64, which is the *stochastic* in stochastic gradient descent.

![Training and test loss by epoch, and test accuracy against the logistic regression line](/images/splice/mlp_training.png)

It reached **94.0%** on the test set. Logistic regression had 93.7%. That, I'm afraid, is a tie.

The loss plot is the more educational part. Training loss keeps falling for all 40 epochs while test loss bottoms out around epoch 7 and then creeps back upward. That's overfitting, happening live: the 32 hidden units have enough capacity to start memorising the training sequences, and the extra confidence is mostly spent on examples the model already gets right, which is why test accuracy barely moves. It's the textbook reason for early stopping and regularisation, and it's far more convincing to watch than to read about.

As for the tie: Noordewier *et al.* found much the same in 1991, with a plain backprop network reaching 94.5% once it had enough training examples. A splice signal is mostly additive, a few positions each contributing independently, and a linear model captures that almost completely. Extra capacity buys very little without more data or a model with more structure. I'd been half-expecting the neural network to win simply because it was a neural network, and it was a useful correction that it didn't.

## Further reading

- [Morse (2017)](https://stmorse.github.io/journal/Beginner-NN-Part-I.html): much of the theory in this post is derived from here, including the full backpropagation recursion (which I skipped).
- [Introduction to Data Analysis in Python (Part 2)](https://bristol-training.github.io/data-analysis-python-2/) and [Introduction to AI](https://bristol-training.github.io/intro-to-ai/), Jean Golding Institute: two free courses which all the code is based on.
- [Nielsen, _Neural Networks and Deep Learning_](http://neuralnetworksanddeeplearning.com/chap1.html): a longer, free treatment if you want the next step up.
- [Noordewier, Towell & Shavlik (1991)](https://proceedings.neurips.cc/paper/1990/file/8efb100a295c0c690931222ff4467bb8-Paper.pdf): the NIPS paper the dataset was built for.
- [splice-site-nn](https://github.com/timj-wood/splice-site-nn): the notebooks, figures and data (UCI, [CC BY 4.0](https://doi.org/10.24432/C5M888)).

This will be my last post for a while; coursework beckons me.

Thanks for reading!