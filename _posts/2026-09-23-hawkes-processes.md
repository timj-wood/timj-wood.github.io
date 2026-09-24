---
title: "Hawkes processes in epidemiology"
summary: A beginner friendly walkthrough of self-exciting point processes, working towards Juliette Unwin's method for separating imported and local malaria cases.
category: hawkes-processes
math: true
---
**Note:** I am an MSc Bioinformatics student, and a biologist by background; by no means a mathematician. If you spot any mistakes then please let me know!

Let's suppose you work in a country that has almost eliminated malaria, and a new case is reported. This could be from a traveller, who was infected abroad. If there are no subsequent infections, it's a dead end. However, the case might be the first visible link in a local chain where a mosquito bites the infected person and later bites someone else, and so on. These situations obviouslt require different different responses, but the day the initial case is reported they cannot be distinguished.

This is the exact problem which [Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830) aimed to address using *Hawkes processes*. To do so, they separated imported from locally acquired malaria in Yunnan (China) and Eswatini, only using the timing of cases, demoting travel histories to a simple check to their answers. Subsequently, I wanted to understand how this is possible, so here we are.

The rough plan is to cover: 

1. Enough background on random events (point processes) for understanding.
2. The definition of Hawkes processes themselves.
3. Simulation of a "fake" malaria outbreak.
4. What can occur when you apply real data.

# (Some) Background

# Events in time: point processes

To start off: a point process is a random list of event times $t_1, t_2, t_3, \ldots$ on a timeline. Point processes are observed in everywhere, for example reported disease cases, earthquakes, and neuron firings. The counting process $N(t)$ is the running tally of events up to and including time $t$, essentially a staircase that starts at 0 and sequentially increases by 1 at each event. Throughout, we'll assume that no two events can occur at the same moment (final section will cover why, in epidemiology, this is discredited). 

Everything which occurred before time $t$ is the history, $\mathcal{H}(t) = \\{t_i : t_i < t\\}$. The strict inequality is required, because whatever we predict for time $t$ may use events from the past, but not the event we are trying to predit. 

## The Poisson process

The most simplistic point process model is the homogeneous Poisson process, where events occur randomly at a constant rate $\lambda$. In any small interval of width $h$, the change of one event is roughly $\lambda h$, the chance of two or more is negligable, and the next occurrence in one interval is independent of every other interval. 

Three key points: 

- The gaps between events are exponentially distributed with mean $1/\lambda$. The process is *memoryless*: however long since the last event, the wait for the next one is the same. 
- The number of events in a window of length $T$ is $\text{Poisson}(\lambda T)$, so the variance is equal to the mean. Counts which vary more than the mean (known as overdispersion) are an early indicator that events cluster, though the cause cannot be deciphered from this. 
- If the rate changes overtime but is fixed in advance, for example $\lambda(t)$ rising around the holidays when more people travel, this produces an *inhomogeneous* Poisson process. The expected number of events between $a$ and $b$ is then the area under the rate curve, $\int_a^b \lambda(t)\,dt$.

The key limitation of the Poisson process is taht events essentially ignore each other. During an epidemic, where each case causes more transmission, this is not the applicable.

<figure>
  <img src="/assets/posts/hawkes/poisson_vs_hawkes.png" alt="A Poisson process and a Hawkes process with the same background rate">
  <figcaption>Fig 1. A Poisson process and a Hawkes process with the same background rate. The Hawkes intensity jumps at each event and decays, and its events come in bursts.</figcaption>
</figure>

Fig 1 provides an illustration of the Poisson process, and alludes to the Hawkes process. Both processes have the same background rate of 0.1 events per day, but each event in a Hawkes process model causes an increase in frequency which slowly fades. In sooth, events arrive in bursts. 

## The conditional intensity 

The way around *memorylessness* is to let the rate depend on what has previously occurred. This is the **conditional intensity**, where the expected rate of events at time $t$, given the history: 

$$
\lambda^{\ast}(t) = \lim_{h \downarrow 0} \frac{\mathbb{E}\big[\,N(t+h) - N(t) \;\big|\; \mathcal{H}(t)\,\big]}{h}.
$$



