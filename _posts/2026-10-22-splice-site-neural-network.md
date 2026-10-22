---
title: "Creating a (mini) neural network to identify splice regions"
category: machine-learning
tags: [neural-networks, genomics]
---

And so I return to the lands of familiarity (bioinformatics). For my next post, I wanted to have a go at something involving my favourite [buzzwords](https://en.wikipedia.org/wiki/Buzzword), **machine learning** and **genomics**, which leads to the *holy grail*: ***neural networks***. I hadn't touched machine learning before this year, but the more I read into it the more interesting it gets, so I've been wanting to build (_some_) sort of neural network skillset.

I'm not going to go into much of the mathematical theory here, partly because I'm still learning it myself. If you want the theory done properly, 3Blue1Brown's series is the best thing I've found:

<iframe align="center" width="560" height="315" src="https://www.youtube.com/embed/aircAruvnKk" frameborder="0" allowfullscreen></iframe>

For the code, I followed the Jean Golding Institute's [Introduction to AI](https://bristol-training.github.io/intro-to-ai/) course, plus their [Data Analysis in Python (Part 2)](https://bristol-training.github.io/data-analysis-python-2/) course for the baseline. All the figures and numbers here come from the notebooks in [splice-site-nn](https://github.com/timj-wood/splice-site-nn), so have a look there if you want to see exactly what I did (and tell me if I've got something wrong!).

## What is splicing?

When a gene is transcribed, the result is a pre-mRNA that still contains the introns (the non-coding regions). These have to be cut out before the message can move on to the next step, translation. So how do your cells know which parts of the transcript are introns, you might ask? Well, (almost) every intron begins with the two bases "GT" and ends with "AG". It's very close to a universal rule, which is why I picked it for a first neural network: I already know what the answer should look like, so if the model learns something sensible I should be able to spot it, and if it learns something weird I should be able to spot that too.

So the question I'm asking the model is: given a window of 60 bases, is there a splice junction in the middle, and if so which kind? And the thing I'm actually curious about is whether it finds the GT…AG pattern on its own. Before any of that, though, I'll fit some ordinary classifiers so there's something to compare against.

## The data

Conveniently, there's a dataset that suits this nicely, and it's 35 years old. The UCI [splice-junction dataset](https://doi.org/10.24432/C5M888) has 3,190 windows of 60 bases from primate genes, each labelled `EI` (exon to intron, a donor site), `IE` (intron to exon, an acceptor site) or `N` (neither). It was put together for a paper that used neural networks to recognise DNA sequences ([Noordewier *et al.*, 1991](https://proceedings.neurips.cc/paper/1990/file/8efb100a295c0c690931222ff4467bb8-Paper.pdf)), which seemed fitting.

There's one problem: a model can't do anything with a string like `CCAGTCAT`. The usual fix is [**one-hot encoding**](https://www.sciencedirect.com/topics/computer-science/one-hot-encoding), where each position becomes four columns (one per base), with a 1 in the column for the base that's there and 0 in the other three. So each sequence becomes 240 numbers, and "there's a G at position +1" turns into a single feature the model can put a weight on.

A few positions (60 out of 191,400) use [IUPAC ambiguity codes](https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html), e.g. `D` means "A, G or T". I didn't want to throw those rows away, so I spread the 1 evenly over the possible bases: a `D` becomes $\frac{1}{3}$, $\frac{1}{3}$, $\frac{1}{3}$ for A, G and T. I'm not sure that's the standard thing to do, but it seemed reasonable and it only affects a handful of entries.

```python
X = one_hot(df["sequence"])   # (3190, 240) array of 0s, 1s and the odd 1/3
y = df["label"].values        # "EI", "IE" or "N"
```

## A baseline first

Before the neural network stuff, I fitted a couple of simpler classifiers. The first thing to do is hold back a quarter of the data (798 sequences) as a test set that the models never see during training; otherwise you have no idea whether they've learnt anything or just memorised the answers.

```python
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42, stratify=y
)
```

```python
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score

logreg = LogisticRegression(max_iter=2000)
logreg.fit(X_train, y_train)

pred = logreg.predict(X_test)
print(f"Logistic regression test accuracy: {accuracy_score(y_test, pred):.3f}")
```

Logistic regression got **93.7%** of the test sequences right, which was higher than I expected for something this simple. I also tried $k$-nearest neighbours, picking $k$ by cross-validation as the course does:

```python
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GridSearchCV

hyperparameters = {"n_neighbors": [1, 3, 5, 7, 11, 15, 21, 31]}
model = GridSearchCV(KNeighborsClassifier(), hyperparameters)
model.fit(X_train, y_train)

print(model.best_estimator_)
print(f"kNN test accuracy: {model.score(X_test, y_test):.3f}")
```

That managed 85.1%. I think the gap is because of what one-hot vectors look like: two random 60-base sequences are hardly ever "close" to each other in a 240-dimensional space, so the nearest neighbours aren't especially informative. (If someone has a better explanation I'd like to hear it.)

The part I found more interesting is what logistic regression actually learnt. It has one weight per (position, base) for each class, so you can reshape the 240 weights into a 4 × 60 grid and plot it:

```python
W = logreg.coef_[list(logreg.classes_).index("EI")]   # 240 weights for the EI class
grid = W.reshape(60, 4).T                               # 4 bases × 60 positions
```

![Logistic regression weights for the EI and IE classes by position and base](/images/splice/logreg_weights.png)

For `EI`, the biggest weights are T at +2 and G at +1. For `IE`, they're G at −1 and A at −2, with a stretch of C and T favoured a bit further upstream, which I believe is the [polypyrimidine tract](https://en.wikipedia.org/wiki/Polypyrimidine_tract) that sits before acceptor sites. So the humble linear model has already found GT…AG. Worth keeping in mind, since the neural network has to at least match this to be worth the effort.

## Building a neuron from scratch

Here's my understanding of what a single artificial neuron does. First it **aggregates**: it takes a weighted sum of its inputs and adds a bias. Then it **activates**: it pushes that sum through some non-linear function.

$$
z = \sum_j w_j x_j + b,
\qquad a = \sigma(z) = \frac{1}{1 + e^{-z}}.
$$

The sigmoid $\sigma$ squashes any number into the range 0 to 1, so we can read the output as a probability. To keep things simple I started with a yes/no version of the problem: "is there *any* splice junction in this window?" (`EI` or `IE` gives 1, `N` gives 0).

To train it we need a **loss**, i.e. a number that says how wrong a prediction was. For a probability the usual choice is binary cross-entropy,

$$
\ell = -[y \log a + (1-y)\log(1-a)],
$$

which (as far as I understand it) punishes a confident wrong answer much more than a hesitant one.

Then we need to know how the loss changes when we nudge a weight. The weight $w_j$ affects $z$, which affects $a$, which affects $\ell$, so by the chain rule we multiply the three sensitivities together. Each one on its own is fine: $\partial z/\partial w_j = x_j$, $\partial a/\partial z = a(1-a)$ for the sigmoid, and $\partial \ell/\partial a = (a-y)/(a(1-a))$ for cross-entropy. The last two mostly cancel, which leaves

$$
\frac{\partial \ell}{\partial w_j} = (a - y)\, x_j, \qquad
\frac{\partial \ell}{\partial b} = a - y.
$$

So the gradient for each weight is just "how wrong was I?" times "was that feature on?". If the neuron was confidently wrong about a sequence that had a G at +1, the weight for G at +1 gets a big nudge; if that feature wasn't there, its weight is left alone. I gather the cancellation isn't an accident, and that sigmoid and cross-entropy are usually paired for this reason. Either way, if you average this over all 2,392 training sequences and write it as a matrix product, the whole training loop is only a few lines:

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

After 400 steps of plain gradient descent, this got 93.4% of the test set right.

Then I built the same neuron the way the course does, with `nn.Linear(240, 1)`, `BCEWithLogitsLoss` and `loss.backward()`, using the same learning rate, the same number of steps and the same zero initialisation:

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

The two loss curves sit on top of each other, and the final weights agree to within $2 \times 10^{-7}$. This was quite reassuring. It means `loss.backward()` isn't doing anything mysterious, at least for one neuron: it's the same chain rule I did by hand, just done automatically. My understanding is that backpropagation is this same idea applied layer by layer, storing the intermediate values as it goes so nothing gets computed twice, but I'll admit I haven't tested that part myself.

Anyway, having done it once by hand, I was quite happy to let PyTorch do it from here on.

## What did the neuron learn?

Because a single neuron has one weight per (position, base), its 240 weights can be plotted in the same way as the logistic regression ones:

![The neuron's weights by position and base; red is evidence for a junction, blue against](/images/splice/neuron_weights.png)

The four biggest weights, in order, are G at −1, A at −2, T at +2 and G at +1. That's AG on the left of the junction and GT on the right. Since this neuron was only asked whether there was a junction of *either* kind, it's picked up both signals at once. I was hoping for this, but it was still quite nice to see it appear from nothing more than 240 columns of 0s and 1s.

## Adding a hidden layer

A single neuron can only add up independent bits of evidence. A hidden layer should let the network combine features (something like "G at +1 *and* T at +2 *and* lots of C/T before it"), which a plain weighted sum can't do. So I built a small multi-layer perceptron: 240 inputs, a hidden layer of 32 units with ReLU activation ($f(z) = \max(0, z)$), and 3 outputs for `EI`, `IE` and `N`, trained with mini-batch gradient descent for 40 epochs.

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

Two things I learnt along the way. First, you can't start all the weights at zero this time: every hidden unit would compute the same thing, get the same gradient and stay identical forever (and with ReLU they'd all be stuck at zero), so PyTorch's random initialisation is used instead. Second, rather than one gradient step over all 2,392 sequences per epoch, the network takes a step after every 64, which is the "stochastic" in stochastic gradient descent.

![Training and test loss by epoch, and test accuracy against the logistic regression line](/images/splice/mlp_training.png)

It reached **94.0%** on the test set. Logistic regression had 93.7%. So, a tie, more or less; a difference that small could easily be down to the particular test split.

The loss plot turned out to be the more useful part. Training loss keeps dropping for all 40 epochs, but test loss bottoms out around epoch 7 and then slowly climbs. I think this is overfitting: the 32 hidden units have enough room to start memorising the training sequences, and the extra confidence mostly goes on examples it already gets right, which would explain why test accuracy barely moves. I'd read about early stopping and regularisation before, but it was a lot clearer seeing it happen in my own plot.

As for the tie, I don't think it's that surprising in hindsight. Noordewier *et al.* found something similar back in 1991, with a plain backprop network reaching about 94.5% once it had enough training examples. My guess is that the splice signal is mostly additive (a few positions each contributing on their own), so a linear model gets most of the way there, and the extra capacity doesn't have much left to find. I'd half expected the neural network to win just because it was a neural network, so this was a useful correction.

## Further reading

- [Morse (2017)](https://stmorse.github.io/journal/Beginner-NN-Part-I.html): a lot of my understanding of the theory came from this post, including the full backpropagation recursion that I skipped.
- [Introduction to Data Analysis in Python (Part 2)](https://bristol-training.github.io/data-analysis-python-2/) and [Introduction to AI](https://bristol-training.github.io/intro-to-ai/), Jean Golding Institute: the two free courses all the code is based on.
- [Nielsen, _Neural Networks and Deep Learning_](http://neuralnetworksanddeeplearning.com/chap1.html): a longer free treatment if you want the next step up.
- [Noordewier, Towell & Shavlik (1991)](https://proceedings.neurips.cc/paper/1990/file/8efb100a295c0c690931222ff4467bb8-Paper.pdf): the NIPS paper the dataset was built for.
- [splice-site-nn](https://github.com/timj-wood/splice-site-nn): the notebooks, figures and data (UCI, [CC BY 4.0](https://doi.org/10.24432/C5M888)).

This will be my last post for a while; coursework beckons.

If you spot any mistakes or have questions, feel free to email me. Thanks for reading!