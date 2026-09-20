---
title: "Events that cause events: Hawkes processes and their applications"
summary: A literature-review style walk through self-exciting point processes, centred on Juliette Unwin's malaria modelling.
category: lit review
math: true
---
**Note:** I am a biologist by background, so there is likely many errors in this post. I am by no means a mathematician, please bear with my elementary LaTeX skills. 

# 0. Epilogue
Let's say you want to eradicate malaria - modelling the spread of a disease is crucial in doing so. However, it is always probable that a traveller, from a place where the disease is common, will bring it with them. Soon after, if there are no more recorded cases, then nothing has come of it - a dead end. If a handful of new cases appear within the same area/district, then you have a causality - mosquitoes have picked up the parasite and moved on. Deciding between these two outcomes is the decisive *'game'* which you must play, and often the only evidence is an unreliable list of travel dates. 

Therefore, this post will be about a class of models built kind of evidence: **Hawkes processes**.

# 1. Prerequisites
## 1.1 Point processes / counting processes
A point process can be described as a random set of event times $$t_1, t_2, t_3, ...$$ on a timeline. These can be found everywhere, for example reported case times in epidemiology, earthquakes, and neuron firings. The counting process $$N(t)$$ is the running tally up to $$t$$, which is essentially a staircase that starts at 0 and increases by 1 at each sequential event ([Laub *et al.,* 2015](https://arxiv.org/pdf/1507.02822)).

Think of point processes as a list of times, and counting processes as the staircase. Everything before time $$t$$ is known as the *history*, which is wrriten as $$\mathcal{H}(t)$$.

A biological example, for providence: imagine the events are mutations fixing along a lineage, where the point process is their positions (or times) and $$N(t)$$ is the cumulative mutation count. 

## 1.2 Homogeneous Poisson process

The simplest point process and natural null model invovles events that occur at a completely random and constant rate $$\mu$$. There are three concrete facts associated with this process. First, the constant rate occurs in a tiny interval of width $$h$$, the chance of the event occuring is $$\approx \lambda h$$. Second, there is independence between non-overlapping intervals. Finally, exponential waiting times equate to a lack of memory in the process, this is because gaps are exponential ($$\lambda$$) with mean $$1/\lambda$$, and time already waited gives no indication about the remaining wait. The count in a window of length $$T$$ is Poisson with mean $$\lambda T$$. 

Whilst this process is incredibly useful at modelling random events in various fields of STEM, such as queuing theory and reliability engineering; the key limitation is that, in real life, events tend to cluster - one event will lead to other events - which the Poisson process cannot capture. This is where Hawkes processes come into play.

## 1.3 Inhomogeneous Poisson process





