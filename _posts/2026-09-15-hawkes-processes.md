---
title: Hawkes processes and malaria modelling
summary: A literature-review style walk through self-exciting point processes, centred on Juliette Unwin's malaria modelling.
category: lit review
math: true
---
Note: I am by no means a mathematician. I am a biologist by background, so there is likely many errors in this post.

# 1. Introduction
Let's say you want to eradicate malaria - modelling the spread of a disease is crucial in doing so. However, it is always probable that a traveller, from a place where the disease is common, will bring it with them. Soon after, if there are no more recorded cases, then nothing has come of it - a dead end. If a handful of new cases appear within the same area/district, then you have a causality - mosquitoes have picked up the parasite and moved on. Deciding between these two outcomes is the decisive *'game'* which you must play, and often the only evidence is an unreliable list of travel dates. 

Therefore, this post will be about a class of models built kind of evidence: **Hawkes processes**.

# 2. Prerequisites
## 2.2 Point processes / counting processes
A point process can be described as a random set of event times $$t_1, t_2, t_3, ...$$ on a timeline. These can be found everywhere, for example reported case times in epidemiology, earthquakes, and neuron firings. The counting process $$N(t)$$ is the running tally up to $$t$$, which is essentially a "staircase" that starts at 0 and increases by 1 at each sequential event (Laub *et al.,* 2015).

Think of point processes as a list of times, and counting processes as the staircase. Everything before time $$t$$ is known as the *history*, which is wrriten as $$\mathcal{H}(t)$$.

A biological example, for providence: imagine the events are mutations fixing along a lineage, where the point process is their positions (or times) and $$N(t)$$ is the cumulative mutation count. 




Hawkes processes are a particularly interesting field of Mathematics - a stochastic process that bridges the gap between statistical and mechanistic models, and is extremely applicable to a diverse range of fields. Hawkes processes are very well deployed in finance, for example they are estimate [transactional data volatility](https://arxiv.org/html/1502.04592v2), but what I will focus on in this post is their applications to epidemiology. 

A Hawkes process has conditional intensity

$$\lambda(t) = \mu + \sum_{t_i < t} \phi(t - t_i),$$

where $\mu$ is the background rate and $\phi$ the excitation kernel.
