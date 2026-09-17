---
title: Hawkes processes and malaria
summary: A literature-review style walk through self-exciting point processes, centred on Juliette Unwin's malaria modelling.
category: lit review
math: true
---
Let's say you want to eradicate malaria - modelling the spread of a disease is crucial in doing so. However, it is always probable that a traveller, from a place where the disease is common, will bring it with them. Soon after, if there are no more recorded cases, then nothing has come of it - a dead end. If a handful of new cases appear within the same area/district, then you have a causality - mosquitoes have picked up the parasite and moved on. Deciding between these two outcomes is the decisive *'game'* which you must play, and often the only evidence is an unreliable list of travel dates. 

Therefore, this post will be about a class of models built kind of evidence: **Hawkes processes**.

# 1. Background

Hawkes processes are a particularly interesting field of Mathematics - a stochastic process that bridges the gap between statistical and mechanistic models, and is extremely applicable to a diverse range of fields. Hawkes processes are very well deployed in finance, to analyse market volatility, but what I will focus on in this post is their applications to epidemiology. 





A Hawkes process has conditional intensity

$$\lambda(t) = \mu + \sum_{t_i < t} \phi(t - t_i),$$

where $\mu$ is the background rate and $\phi$ the excitation kernel.
