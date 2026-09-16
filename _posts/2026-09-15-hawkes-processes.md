---
title: Hawkes processes and malaria
summary: A literature-review style walk through self-exciting point processes, centred on Juliette Unwin's malaria modelling.
category: lit review
math: true
---
Work in progress.

A Hawkes process has conditional intensity

$$\lambda(t) = \mu + \sum_{t_i < t} \phi(t - t_i),$$

where $\mu$ is the background rate and $\phi$ the excitation kernel.
