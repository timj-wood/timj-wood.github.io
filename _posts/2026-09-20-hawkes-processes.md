---
title: "Hawkes processes in epidemiology"
summary: A walkthrough of self-exciting point processes, inspired by Juliette Unwin's paper on malaria modelling.
category: hawkes-processes
math: true
---
Work in progress. 

**Note:** I am a biologist by background, and by no means a mathematician. Please bear with my elementary LaTeX skills and let me know if you spot any mistakes.

Suppose you work in a country that has almost eliminated malaria, and a new case is reported. It might come from a traveller who was infected abroad and brought the parasite home with them. If nothing follows, it was a dead end. On the other hand, it may be the first visible link in a local chain, where a mosquito bites an infected person, subsequently biting someone else, and transmission soon becomes widespread.

These two situations require completely different responses, yet on the day the case is reported they look identical. Often the only evidence is a list of cases, the dates they occurred, and some sketchy travel histories.

However, what you can do is identify a temporal pattern. Imported cases arrive more or less at random. Locally transmitted cases arrive in clusters, because each one makes further cases more likely for a while. Distinguishing these two patterns, soley relying on the timing of events, is precisely what a Hawkes process has been developed for.

This post introduces Hawkes processes from scratch, building up to their applications in epidemiology. 

# 1. Prerequisites

## 1.1. Point processes / counting processes

A point process can be described as a random set of event times $$t_1, t_2, t_3, ...$$ on a timeline. These can be found everywhere, for example reported case times in epidemiology, earthquakes, and neuron firings. The counting process $$N(t)$$ is the running tally up to $$t$$, which is essentially a staircase that starts at 0 and increases by 1 at each sequential event ([Laub *et al.,* 2015](https://arxiv.org/pdf/1507.02822)).

Think of point processes as a list of times, and counting processes as the staircase. Everything before time $$t$$ is known as the *history*, which is wrriten as $$\mathcal{H}(t)$$.

To give a biological example, imagine the events are mutations fixing along a lineage, where the point process is their positions (or times) and $$N(t)$$ is the cumulative mutation count. 

## 1.2. Homogeneous Poisson process

The simplest point process and natural null model invovles events that occur at a completely random and constant rate $$\mu$$. There are three concrete facts associated with this process. First, the constant rate occurs in a small interval of width $$h$$, the chance of the event occuring is $$\approx \lambda h$$. Second, there is independence between non-overlapping intervals. Finally, exponential waiting times equate to a lack of memory in the process, this is because gaps are exponential ($$\lambda$$) with mean $$1/\lambda$$, and time already waited gives no indication about the remaining wait. The count in a window of length $$T$$ is Poisson with mean $$\lambda T$$. 

Whilst this process is incredibly useful at modelling random events in various fields of STEM, such as queuing theory and reliability engineering; the key limitation is that, in real life, events tend to cluster - one event will lead to other events - which the Poisson process cannot capture. This is the core drive for Hawkes processes.

## 1.3. Inhomogeneous Poisson process

In this case of the Poisson process, the rate changes over time, where the function $$\lambda(t)$$ replaces the constant $$\lambda$$. Counts in disjoint intervals remain independent, but the background tendency rises and falls. For example, the number of flu cases varies with the seasons. The number of events between $$a$$ and $$b$$ is Poisson distributed with mean

$$\int_a^b \lambda(t)\,dt.$$

Crucially, the rate varies because of external forces (seasons, weather), not because previous events change it. This distinction matters because a time-varying rate produces clusters of events, and so does self-excitation, so the two are easily confused in data. Fitting a self-exciting model to seasonally driven data, or vice versa, will lead to the incorrect conclusion for the driving processes.


## 1.4. The conditional intensity function 

The most important concept in this post is the **conditional intensity**. This is the instantaneous expected rate of events at time $t$, given the entire history of the process up to that point.

$$
\lambda^*(t) = \lim_{h \downarrow 0} \frac{\mathbb{E}\big[\,\text{events in } (t,\, t+h] \;\big|\; \mathcal{H}(t)\,\big]}{h}
$$

The asterisk is shorthand for "conditional on the history". Breaking down the formula: 

- $\lambda^*(t)$ denotes the present rate, given everything that has previously occured.
- $\mathbb{E}[\,\cdots \mid \mathcal{H}(t)\,]$ is the expected number of events, given the history $\mathcal{H}(t)$.
- $\lim_{h \downarrow 0}$ shrinks the window to zero, thus detailing "the next instant".

In practical terms, $$\lambda^*(t)\,dt$$ is $$\approx$$ the probability of an event occurring in the next instant.

For a homogeneous Poisson process the history is irrelevant, so $$\lambda^*(t) = \lambda$$. For a Hawkes process, $$\lambda^*(t)$$ increases at each event and then decays.

The integral of the intensity,

$$
\Lambda(t) = \int_0^t \lambda^*(s)\,ds,
$$

is known as the **compensator**: the cumulative number of events the model expects to have seen by time $t$.

# 2. What is a Hawkes processes?

## 2.1. History and motivation

Hawkes processes are named after the British statistician Alan G. Hawkes, who introduced them in a pair of papers ([Hawkes, 1971a](https://academic-oup-com.bris.idm.oclc.org/biomet/article/58/1/83/224809); [Hawkes, 1971b](https://academic-oup-com.bris.idm.oclc.org/jrsssb/article/33/3/438/7027167)). His idea was to write down a point process in which events are "self-exciting", where each event temporarily raises the rate of future events. Soon after, [Hawkes and Oakes (1974)](https://www-cambridge-org.bris.idm.oclc.org/core/journals/journal-of-applied-probability/article/abs/cluster-process-representation-of-a-selfexciting-process/E836A3D07D808068E2F9F3E7E366B081) demonstrated that the same process can be viewed as a family tree, in which some events arrive spontaneously and each event goes on to produce "offspring" events of its own. This concept of branching has turned out to be the most intuitive way to think about the model.

The model soon became famous for its applications in seismology. An earthquake triggers aftershocks, which can trigger aftershocks of their own, and this idea was built into the Epidemic-Type Aftershock Sequence (ETAS) model ([Ogata, 1988](https://www-tandfonline-com.bris.idm.oclc.org/doi/abs/10.1080/01621459.1988.10478560?casa_token=W-srgC9krLYAAAAA:3HduSWc9el-inU3aIGRQcXjM6qyHG7mSpMB156f02R4JE8XDi_soQjGcxtBbyL8VQozYUJWXgp6IxQ); [Ogata, 1998](https://link-springer-com.bris.idm.oclc.org/article/10.1023/A:1003403601725)). To describe earthquakes, seisomologists borrowed the language of contagion from epidemiology. Decades later, epidemiologists would adopt this seismological model to describe contagion. 

## 2.2. Core definition 

A Hawkes process is defined by its conditional intensity. In the simplest case, with a single stream of events, it is

$$
\lambda^*(t) = \mu + \sum_{t_i < t} \phi(t - t_i).
$$

$$\mu > 0$$ is the background rate: the rate at which events occur spontaneously, regardless of what has occured previously. For example, in malaria, this would be the reported cases. 

$$\phi(\cdot) \geq 0$$ is the triggering kernel: the extra rate presently contributed by a past event. It's input, $$t-t_i$$, is the time elapsed since the event at $$t_i$$. Typically, a kernel will start high and decay, so an event's influence is strongest immediately after its occurence, and then fades. In the case of malaria, this is local transmission. 

The sum runs over every event before time $$t$$. Each past event adds a decaying spike to the rate, which stack on top of the background. 

This is the definition of **self-exciting**, where an event raises $$\lambda^*(t)$$, increasing the likelihood of another event occuring, thus raising $$\lambda^*(t)$$ again. This feedback is what produces clusters of events in time, and unlike the inhomogeneous Poisson process (Section 1.3), the clustering is generated by the actual events occuring, not external forces. 

## 2.3. Formalising kernels 

The kernel $$\phi$$ describes how an event's influence plays out overtime, and different choices of kernel will produce different model outcomes. In the context of this review, there are three worth stating. 

**The exponential kernel** is the original and most common kernel:

$$
\phi(t) = \alpha e^{-\beta t}, \qquad \alpha, \beta > 0
$$

Each event instantly raises the intensity by $$\alpha$$, where the increase then decays at rate $$\beta$$, so an event's influence lasts for $$\approx$$ $$1/\beta$$ units of time. The popularity of this kernal lies in its practicality - the exponential's lack of memory means that the model can be quickly fitted to data. 

**The power-law kernel** decays at a much slower rate:

$$
\phi(t) = \frac{k}{(c + t)^p}, \qquad k, c > 0, \; p > 1.
$$

Here $$k$$ sets the overall strength of triggering, $$p$$ controls how quickly an event's influence fades (larger $$p$$ means faster decay), and $$c$$ is a small offset that keeps the kernel finite at $$t = 0$$, so that each event raises the intensity by $$k / c^p$$. The condition $$p > 1$$ ensures that the total influence of a single event is finite.

This kernel is "heavy-tailed", so its influence fades at such a slow rate that events from long ago can still trigger new ones, whereas an exponential kernel's influence is effectively gone after a few multiples of $$1/\beta$$.

**A kernel can also be matched to a specific process.** $$\phi$$ is not necessarily a simple formula - it can be a flexible step function estimated from the data, or a shape chosen from prior knowledge. The second option is particularly important in epidemiology. 

The time between one person becoming infected and them infecting someone else is called the *generation interval*, and for many diseases its distribution has been measured and is well described by a gamma or lognormal curve. Using this curve as the kernel integrates the disease's biology directly into model - improving its efficacy.


**Quick summary:** use the exponential when speed and simplicity matter, the power law when influence persists over long periods, and a bespoke kernel when you possess data regarding the delays.

## 2.4. The branching ratio / reproduction number 

The branching ratio quantifies the number of further events which one event can trigger. As previously defined, the expected number of events produced by a rate was the area under the rate curve. An event adds $$\phi$$ to the rate, so the expected number of events it triggers is the area under the kernel: 

$$
n = \int_0^\infty \phi(t)\,dt.
$$

