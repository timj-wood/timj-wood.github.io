---
title: Hawkes processes and malaria
summary: A literature-review style walk through self-exciting point processes, centred on Juliette Unwin's malaria modelling.
category: lit review
math: true
---
Hawkes processes are a particularly interesting field of Mathematics - a stochastic process that bridges the gap between statistical and mechanistic models, and are extremely applicable to a diverse range of fields. Hawkes processes are very well deployed in finance, to analyse market volatility, but what I will focus on in this post is their applications to epidemiology. 



A Hawkes process has conditional intensity

$$\lambda(t) = \mu + \sum_{t_i < t} \phi(t - t_i),$$

where $\mu$ is the background rate and $\phi$ the excitation kernel.
