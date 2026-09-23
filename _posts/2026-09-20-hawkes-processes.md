---
title: "Hawkes processes in epidemiology"
summary: A walkthrough of self-exciting point processes, inspired by Juliette Unwin's paper on malaria modelling.
category: hawkes-processes
math: true
---
Work in progress. 

**Note:** I am a biologist by background, and by no means a mathematician. Please bear with my elementary LaTeX skills and let me know if you spot any mistakes.

Suppose you work in a country that has almost eliminated malaria, and a new case is reported. It might come from a traveller who was infected abroad and brought the parasite home with them. If nothing follows, it was a dead end. On the other hand, it may be the first visible link in a local chain, where a mosquito bites an infected person, subsequently biting someone else, and transmission soon becomes widespread. Notably, local transmission doesn't yield an immediate increase in the likelihood of further cases - which the kernel informs (Section 1.4).

These two situations require completely different responses, yet on the day the case is reported they look identical. Often the only evidence is a list of cases, the dates they occurred, and some travel histories - partial labels of which cases are imported, which Hawkes models can exploit.

However, what you can do is identify a temporal pattern. We can assume that imported cases are not random, due to holidays, pilgrimages, school terms etc. Locally transmitted cases arrive in clusters, because each one makes further cases more likely for a while. Distinguishing these two patterns, solely relying on the timing of events, is exactly the type of problem that Hawkes process are well suited to. 

This post introduces Hawkes processes from scratch, building up to their applications in epidemiology. 

# 1. Prerequisites

## 1.1. Point processes / counting processes

A point process can be described as a random set of event times $$t_1, t_2, t_3, ...$$ on a timeline. These can be found everywhere, for example reported case times in epidemiology, earthquakes, and neuron firings. The counting process $$N(t)$$ is the running tally up to $$t$$, which is essentially a staircase that starts at 0 and increases by 1 at each sequential event ([Laub *et al.,* 2015](https://arxiv.org/pdf/1507.02822)). Formally, a counting process satisfies $$N(0) = 0$$, takes non-negative integer values, is non-decreasing, and is right-continuous, meaning that at an event time $$t_i$$, the count already includes that event. We also assume that the process is simple, where no two events occur at the exact same time, so each step of the staircase has a height of 1. 

Think of point processes as a list of times, and counting processes as the staircase. Everything before time $$t$$ is known as the *history*, which is wrriten as $$\mathcal{H}(t) = \{t_i:t_i<t\}$$, with a strict inequality. The intensity at $$t$$ must not depend on whether an event happens at $$t$$ itself.

Returning to the example of malaria, the events are the dates on which cases are reported. The point process is the list of those dates, and $$N(t)$$ is the total number of cases reported by time $$t$$.

## 1.2. Homogeneous Poisson process

The simplest point process, and the natural null model, is one in which events occur completely at random at a constant rate $$\lambda$$. It is defined by two properties. First, in a small interval of width $$h$$,

$$
\begin{aligned}
&P(\text{one event in } (t, t+h]) = \lambda h + o(h), \\
&P(\text{two or more events in } (t, t+h]) = o(h),
\end{aligned}
$$

where $$o(h)$$ denotes terms that become negligible relative to $$h$$ as $$h \to 0$$. The second condition says events occur one at a time. Second, the numbers of events in non-overlapping intervals are independent.

Because the rate is constant and intervals are independent, the gaps between events are independent and exponentially distributed with rate $$\lambda$$ (mean $$1/\lambda$$). The exponential is the only continuous distribution that is *memoryless*: however long you have already waited, the remaining wait has the same distribution. The count in a window of length $$T$$ is $$N(T) \sim \text{Poisson}(\lambda T)$$, so its variance equals its mean. This gives a simple first diagnostic: clustered data, whatever the cause, are *overdispersed*, with variance greater than the mean.

Two further properties will be useful later:
- **Superposition**: combining independent Poisson processes gives a Poisson process whose rate is the sum of their rates.
- **Thinning**: keeping each event independently with probability $$p$$ gives a Poisson process with rate $$p\lambda$$. Superposition underpins the branching view of Hawkes processes, and thinning is the basis of Ogata's algorithm for simulating them.

The Poisson process is widely used across STEM, from queuing theory to reliability engineering. Its key limitation is that events cannot influence one another. In an epidemic, however, each case can cause further cases. Capturing this *self-excitation* is the core motivation for Hawkes processes.

## 1.3. Inhomogeneous Poisson process

In this case, the constant rate $$\lambda$$ is replaced by a function of time, $$\lambda(t)$$, which is non-negative and deterministic (fixed in advance, rather than random). Counts in non-overlapping intervals remain independent, but the background tendency for events rises and falls. For example, the number of flu cases varies with the seasons. The number of events between $$a$$ and $$b$$ is Poisson distributed with mean

$$
\int_a^b \lambda(t)\,dt.
$$

Crucially, the rate varies because of external forces (seasons, weather), not because previous events change it. This distinction matters because a time-varying rate produces clusters of events, and so does self-excitation, so the two are easily confused in data. Fitting a self-exciting model with a constant background to seasonally driven data will attribute the seasonal peaks to transmission, inflating the apparent strength of self-excitation ([Filimonov & Sornette, 2015](https://www-tandfonline-com.bris.idm.oclc.org/doi/full/10.1080/14697688.2015.1032544)). Conversely, fitting an inhomogeneous Poisson process to self-exciting data will attribute transmission to the background, underestimating it.

If the rate is itself random, for example driven by unobserved fluctuations in mosquito abundance, the result is a *Cox* (or doubly stochastic) process. This is a third route to clustering, and arguably the hardest to separate from self-excitation. In practice, these mechanisms also interact: in malaria, a wet season increases not only the background rate but also how much onward transmission each case generates. We return to this when discussing time-varying baselines in Hawkes models.

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
n = \int_0^\infty \phi(s)\,ds.
$$

The sequential events triggered by another are described as the "offspring", in which each event has, on average, $$n$$ offspring. Each sequential event has its own offspring, so one spontaneous event is followed by $$\approx$$ $$n$$ events in the first generation, $$n^2$$ in the second, $$n^3$$ in the hitds, etc. The value of $$n$$ determines the subsequent response:

**When $$n < 1$$** each generation is smaller than the previous, so each chain of events will eventually die out. The expected size of the entire cluster (including original event) is $$1 + n + n^2 + \cdots = 1/1(1-n)$$ ([Laub *et al.*, 2025](https://arxiv.org/pdf/1507.02822)) For example, when $$n = 0.8$$ each spontaneous event leads to a cluster of five events, on average. The process has a steady state that it typically returns to.

**When $$n = 1$$** each event will replicate itself exactly (most of the time). In this case, individual chains will still die out, but they can reach enormous sizes before doing so. Even the slightest change in $$n$$ will tip it into one regime (e.g., settling) or the other (exploding). For example, in epidemiology, if $$R = 1$$, a disease is neither disappearing nor spreading. 

**When $$n = 1$$** each generation of events is larger than the last, and the number of events grows without limit. 

As stated, this directly links to the example of the reproduction number $$R$$ in epidemiology. To make this correspondence, the kernel can be split into two 

$$
\phi(t) = R \, g(t), 
$$

where $$g(t)$$ is the generation interval distribution. $$g$$ is a probability distribution so its area is 1, meaning $$n = R$$. Therefore, applying the determinants above $$R < 1$$, $$R = 1$$, and $$R > 1$$. $$R$$ varying overtime produces the time varying reproduction number $$R_t$$. 


## 2.5. Immigrant-offspring representation 

