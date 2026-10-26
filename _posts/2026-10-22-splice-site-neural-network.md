---
title: "Creating a (mini) neural network to identify splice regions"
redirect_from: /posts/splice-site-neural-network/
category: machine-learning
tags: [neural-networks, genomics]
---

For my next post, I wanted to have a go at something involving my favourite bioinformatics buzzwords, **machine learning** and **genomics**, which both lead to the *holy grail*: [***neural networks***](https://en.wikipedia.org/wiki/Neural_network_(machine_learning)). I had never touched machine learning before this, but it is more relevant to computer science now than ever, so building some sort of neural network skillset has been a top priority for me.

Note that this is not a post about neural network theory, partly because I'm still learning it myself. If you want the theory done properly, [3Blue1Brown's tutorials](https://www.3blue1brown.com/lessons/neural-networks/) are fantastic.

For the code, I have roughly followed the Jean Golding Institute's [Data Analysis in Python (Part 2)](https://bristol-training.github.io/data-analysis-python-2/) and [Introduction to AI](https://bristol-training.github.io/intro-to-ai/) courses. All the figures and numbers here come from the scripts in the repo [splice-site-nn](https://github.com/timj-wood/splice-site-nn), so have a look there if you want to see exactly what I did (and, I'm sure, what I did wrong).

## Gene splicing and rationale

When a gene is transcribed, the result is a pre-mRNA that still contains the introns (the non-coding regions). These have to be cut out before the message can move on to the next step, translation. So how do your cells know which parts of the transcript are introns, you might ask? Well, (almost) every intron begins with the two bases "GT" and ends with "AG". It's very close to a universal rule, which makes it a great target for my first neural network. Since I already know what the answer should look like, I should be able to see if the model learns something sensible, and spot it if it learns something weird... because it will be wrong.

So the question I'm asking the model is: given a window of 60 bases, is there a splice junction in the middle, and if so, which kind? The end goal is a model that learns the GT…AG pattern by itself. Before any of that, though, I'll fit some ordinary classifiers so there's something to compare the neural network against.

## The data

The UCI [splice-junction dataset](https://doi.org/10.24432/C5M888) has 3,190 windows of 60 bases from primate genes, each labelled `EI` (exon to intron, a donor site), `IE` (intron to exon, an acceptor site) or `N` (neither). It was put together for a paper that used neural networks to recognise DNA sequences ([Noordewier *et al.*, 1991](https://proceedings.neurips.cc/paper/1990/file/8efb100a295c0c690931222ff4467bb8-Paper.pdf)), which is quite fitting.

There is a problem with the data, though: a model can't do anything with a string like `CCAGTCAT`. The standard fix is [one-hot encoding](https://www.geeksforgeeks.org/machine-learning/ml-one-hot-encoding/). Each position becomes four columns (one per base), with a 1 in the column for the base that's there and 0 in the other three, so each sequence becomes a row of 240 1s and 0s. The model gives each column its own weight, and each column now represents one base at one position, so the model can learn how much evidence a G at position +1 provides, independently of every other base and position.

A few positions (60 out of 191,400) use [IUPAC ambiguity codes](https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html); for example, `D` means "A, G or T". I didn't want to discard this data, so I spread the 1 evenly over the possible bases: a `D` becomes $\frac{1}{3}$, $\frac{1}{3}$, $\frac{1}{3}$ for A, G and T. I'm not sure if that's the standard approach, but it seemed reasonable, and it only affects a handful of entries.

```python
X = one_hot(df["sequence"])   # (3190, 240) array of 0s, 1s and the odd 1/3
y = df["label"].values        # "EI", "IE" or "N"
```

## Creating a comparable baseline with classical methods

Before the fun neural network stuff, I fitted a couple of simpler classifiers.

First, hold back a quarter of the data (798 sequences) as a test set that the models never see during training. Otherwise there's no way to tell whether they have learnt anything at all.

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

Logistic regression got **93.7%** of the test sequences right, which was higher than I expected for such a simple model. I also tried $k$-nearest neighbours, picking $k$ by cross-validation as the JGI course does:

```python
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GridSearchCV

hyperparameters = {"n_neighbors": [1, 3, 5, 7, 11, 15, 21, 31]}
model = GridSearchCV(KNeighborsClassifier(), hyperparameters)
model.fit(X_train, y_train)

print(model.best_estimator_)
print(f"kNN test accuracy: {model.score(X_test, y_test):.3f}")
```

This was noticeably worse, at 85.1%. I think the gap comes from how distances behave between one-hot vectors in 240 dimensions: two unrelated 60-base sequences differ at most positions, so almost every pair ends up roughly the same distance apart, and the nearest neighbours aren't especially informative. (If someone has a better explanation, please let me know!)

The logistic regression model is simple, but its weights are easy to read. It has one weight per (position, base) for each class, so you can reshape the 240 weights into a 4 × 60 grid and plot it:

```python
W = logreg.coef_[list(logreg.classes_).index("EI")]   # 240 weights for the EI class
grid = W.reshape(60, 4).T                               # 4 bases × 60 positions
```

![Logistic regression weights for the EI and IE classes by position and base](/images/splice/logreg_weights.png)

For `EI`, the biggest weights are T at +2 and G at +1. For `IE`, they're G at −1 and A at −2, with a stretch of C and T favoured a bit further upstream, which is probably the [polypyrimidine tract](https://en.wikipedia.org/wiki/Polypyrimidine_tract) that sits before acceptor sites. So the *humble* linear model has already found GT…AG. That's worth keeping in mind, since the neural network has to at least match it to be worth the effort (otherwise, why am I even bothering?).

## Building a neuron from scratch

*Preface: we are working with a single neuron, you might find in the literature that $z$ is computed to the $i$-th neuron - we only have one.*

A single artificial neuron does two things. First, it takes a weighted sum of its inputs and adds a bias - known as **aggregation**. Then, it pushes that sum through some non-linear function - **activation**.

$$
z = \sum_j w_j x_j + b,
\qquad a = \sigma(z) = \frac{1}{1 + e^{-z}}.
$$

The sigmoid $\sigma$ compresses any number into the range 0 to 1, so the output can be read as a probability. To keep things simple, I started with a binary approach to this problem, i.e., "is there *any* splice junction in this window?" - where `EI` or `IE` gives 1 and `N` gives 0.

To train the model, we first need to measure the difference between model predictions and our desired values - known as loss (what the model got wrong). In this case, smaller values are more favourable. For a binary classification problem (our problem), the standardised choice is binary cross-entropy,

$$
\ell = -[y \log a + (1-y)\log(1-a)],
$$

where $y$ is the binary label (0 or 1) of the observation, and $p$ is the predicted probability of the observation. Lower values equivalate to a better alignment between a model prediction and the true data labels. 

When training a model, it is common for errors to occur, which impact its final accuracy - how loss impacted when a weight is "nudged". This is where backpropagation comes in, by dividing the error between connections. The weight $w_j$ affects $z$, which affects $a$, which affects $\ell$, using the chain rule the sensitivites can be multiplied together. The calculation for each component is straightforward: $\partial z/\partial w_j = x_j$, $\partial a/\partial z = a(1-a)$ for the sigmoid, and $\partial \ell/\partial a = (a-y)/(a(1-a))$ for cross-entropy. The last two mostly cancel out, which leaves

$$
\frac{\partial \ell}{\partial w_j} = (a - y)\, x_j, \qquad
\frac{\partial \ell}{\partial b} = a - y.
$$

So the gradient for each weight is the partial derivative of the loss function with respect to that specific weight. If the neuron was confidently wrong about a sequence that had a G at +1, the weight for G at +1 receives a big nudge; if that feature wasn't there, the weight is unchanged. I gather the cancellation isn't an accident, and that sigmoid and cross-entropy are usually paired for this reason. Either way, if you average this over all 2,392 training sequences and write it as a matrix product, the whole training loop is only a few lines:

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

After 400 steps of plain gradient descent, this scored **93.4%** (on the binary "yes-or-no" task, so it isn't directly comparable to the three-class 93.7% above).

I then built the same neuron the way the JGI course does, with `nn.Linear(240, 1)`, `BCEWithLogitsLoss` and `loss.backward()`, using the same learning rate, the same number of steps and the same zero initialisation:

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

The two loss curves sit on top of each other, and the final weights agree to within $2 \times 10^{-7}$. A reassuring result! It means `loss.backward()` isn't doing anything weird, at least for one neuron: it's the same chain rule I applied by hand, just automated. My understanding is that backpropagation is the same idea applied layer by layer, storing the intermediate values as it goes so nothing gets computed twice, but I'll admit I haven't tested that part myself.

Building a neuron yourself is an insightful (and difficult) exercise, but PyTorch does practically all of this for you... so just use PyTorch.

## What the neuron learnt

Because a single neuron has one weight per (position, base), its 240 weights can be plotted the same way as the logistic regression ones:

![The neuron's weights by position and base; red is evidence for a junction, blue against](/images/splice/neuron_weights.png)

The four biggest weights, in order, are G at −1, A at −2, T at +2 and G at +1. That's AG on the left of the junction and GT on the right. Since this neuron was only asked whether there was a junction of *either* kind (`EI` or `IE`), it picked up both signals at once. I was hoping for this, but it was still nice to see it emerge from nothing more than 240 columns of (mostly) 1s and 0s.

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

Two things I learnt along the way:

1. You can't start all the weights at zero this time. Every hidden unit would compute the same thing, get the same gradient and stay identical forever (and with ReLU they'd all be stuck at zero), so PyTorch's random initialisation is used instead.
2. Rather than taking one gradient step over all 2,392 sequences per epoch, the network takes a step after every 64. That's the "stochastic" part of stochastic gradient descent.

![Training and test loss by epoch, and test accuracy against the logistic regression line](/images/splice/mlp_training.png)

It reached **94.0%** on the test set, against 93.7% for logistic regression. So, unfortunately, it's more or less a tie; a difference that small could easily be down to the particular test split.

The loss plot turned out to be the most useful part. Training loss keeps dropping for all 40 epochs, but test loss bottoms out around epoch 7 and then slowly climbs. I read this as a classic case of overfitting: the 32 hidden units have enough room to start memorising the training sequences, and the extra confidence mostly goes on examples it already gets right, which would explain why test accuracy barely moves. I'd read about early stopping and regularisation before, but it was much clearer seeing it happen in my own plot.

In hindsight, the tie isn't that surprising. Noordewier *et al.* found something similar back in 1991, with a plain backprop network reaching about 94.5% once it had enough training examples. My guess is that the splice signal is mostly additive (a few positions each contributing on their own), so a linear model gets most of the way there and the extra capacity doesn't have much left to find. I'd half expected the neural network to win simply because it's a ***neural network***, so this was a useful correction.

## Further reading

There are lots of useful (and free!) tutorials online that actually cover the theory and give good coding examples. To name a few: [Steven Morse's blog series](https://stmorse.github.io/journal/Beginner-NN-Part-I.html), [Michael Nielsen's online book](http://neuralnetworksanddeeplearning.com/), and [JGI's online course](https://bristol-training.github.io/intro-to-ai/). If you're not from a maths/CS background (I'm not), the theory can get complicated quickly, but if you keep it simple and stick to the code, you can get a long way. My aim for this post wasn't to be thorough, but to show what neural networks can do in the biological sciences, even with a very simple example.

This will be my last post for a while; coursework beckons. If you have any ideas for applied neural networks in biology, please do let me know. Thank you for reading!
