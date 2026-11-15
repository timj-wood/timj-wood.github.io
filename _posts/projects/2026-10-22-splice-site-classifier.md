---
title: "A small neural network for splice-junction classification"
short: "Splice-junction classifier"
summary: "Logistic regression, a single sigmoid neuron and a small multilayer perceptron, trained to find the GT…AG splice signal in 60-base DNA windows."
when: oct 2026
links:
  - { label: repo, url: "https://github.com/timj-wood/splice-site-nn" }
  - { label: write-up, url: /posts/splice-site-neural-network.html }
---

## Approach

Three classifiers were trained to decide whether a 60-base window of primate DNA contains a splice junction at its centre, and of which kind, using the 3,190 labelled sequences of the [UCI splice-junction dataset](https://doi.org/10.24432/C5M888): logistic regression and *k*-nearest neighbours as baselines, a single sigmoid neuron with the gradient descent written by hand and checked against PyTorch, and a multilayer perceptron with one hidden layer of 32 units. A quarter of the data was held out for testing.

## Findings

Logistic regression reached 93.7% test accuracy, the neuron 93.4% on the binary version of the task, and the perceptron 94.0%; *k*-nearest neighbours managed 85.1%. Every model that can be read position by position found the same thing: GT immediately after the junction and AG immediately before it. The perceptron's small lead over logistic regression is within what a different test split could produce, and it began overfitting after about seven epochs. The splice signal in this data appears to be largely additive, so a linear model captures most of it.

The code and figures are in the [repository](https://github.com/timj-wood/splice-site-nn); how each model was built, and what I learnt doing it, is in the [write-up]({{ '/posts/splice-site-neural-network.html' | relative_url }}).
