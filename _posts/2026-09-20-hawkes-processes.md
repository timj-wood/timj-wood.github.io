---
title: "Hawkes processes in epidemiology"
summary: A walkthrough of self-exciting point processes, inspired by Juliette Unwin's paper on malaria modelling.
category: hawkes-processes
math: true
---
Work in progress. 

**Note:** I am a biologist by background, and by no means a mathematician. Please bear with my elementary LaTeX skills and let me know if you spot any mistakes.

Suppose you work in a country that has almost eliminated malaria, and a new case is reported. It might be a traveller infected abroad, a dead end if nothing follows. Or it might be the first visible link in a local chain, where a mosquito bites the infected person and later bites someone else. These situations need very different responses, yet on the day the case is reported they look identical. Often the only evidence is a list of case dates and incomplete travel histories, so the timing of cases has to do much of the work.

There are two temporal patterns. Imported cases arrive at a rate set by external factors such as holidays, pilgrimages and school terms, and are not caused by earlier local cases. Local cases are caused by earlier cases, so each one raises the chance of further cases for a while afterwards. In malaria this rise is delayed by weeks, because the parasite must develop in the mosquito and then incubate in the new host. Both patterns produce clusters, which is what makes them hard to separate.

Hawkes processes are designed for exactly this problem. They treat events as a mix of background events (here, importations) and triggered events (local transmission), where every case, imported or not, can have "offspring". The average number of offspring per case, the branching ratio, plays the role of the reproduction number, and whether it sits below one is the key question for elimination. [Unwin et al. (2021)](https://www.researchgate.net/publication/350572708_Using_Hawkes_Processes_to_model_imported_and_local_malaria_cases_in_near-elimination_settings) used this approach to separate imported from locally acquired malaria in near-elimination settings.

This review introduces Hawkes processes from scratch, building up to their applications in epidemiology.

# 1. Prerequisites

## 1.1. Point processes / counting processes

A point process is a random set of event times $t_1, t_2, t_3, \ldots$ on a timeline. Examples are everywhere, such as reported cases in epidemiology, earthquakes, and neuron firings. The associated counting process $N(t)$ is the running tally of events up to time $t$, a staircase that starts at 0 and steps up by 1 at each event ([Laub *et al.,* 2015](https://arxiv.org/pdf/1507.02822)). i.e., the point process is a list of times, and the counting process is a staircase.

Formally, a counting process satisfies $N(0) = 0$, takes non-negative integer values, is non-decreasing, and is right-continuous, meaning that at an event time $t_i$ the count already includes that event. We also assume only finitely many events occur in any bounded interval, and that the process is *simple*, meaning that no two events occur at exactly the same time, so every step of the staircase has height 1.

Everything that has happened before time $t$ is called the *history*, written $\mathcal{H}(t) = \{t_i : t_i < t\}$. The inequality is strict on purpose. When we define the intensity of a process in Section 1.4, it must depend only on the past, not on whether an event happens at $t$ itself.

In the case of malaria, the events are the dates on which cases occur, the point process is the list of those dates, and $N(t)$ is the total number of cases by time $t$. However, there are two practical limitations to this:
 
 - Surveillance data are usually recorded to the day, so several cases often share a date, which breaks the simplicity assumption. Common fixes include spreading tied events randomly within the day, or using discrete-time versions of the models. 
 - The recorded date is usually the date of diagnosis or reporting, not of infection, so the observed timeline is a shifted and blurred version of the transmission timeline. Both issues resurface when we fit Hawkes models to real data.

## 1.2. Homogeneous Poisson process

The simplest point process, and the natural null model, is one in which events occur completely at random at a constant rate $\lambda$. It is defined by two properties.

First, events occur at a constant rate and one at a time. In a small interval of width $h$, the probability of exactly one event is

$$
P(\text{one event in } (t, t+h]) = \lambda h + o(h),
$$

where $o(h)$ denotes terms that become negligible relative to $h$ as $h \to 0$. The probability of more than one event is negligible:

$$
P(\text{two or more events in } (t, t+h]) = o(h),
$$

so the probability of no event is $1 - \lambda h + o(h)$.

Second, the numbers of events in non-overlapping intervals are independent.

Together, these imply that the gaps between events are independent and exponentially distributed with rate $\lambda$ (mean $1/\lambda$). The exponential is the only continuous distribution that is *memoryless*: however long you have already waited, the remaining wait has the same distribution. The count in a window of length $T$ is $N(T) \sim \text{Poisson}(\lambda T)$, so its variance equals its mean.

This gives a simple first diagnostic. If you split the timeline into windows and the counts vary more than their mean (*overdispersion*), the data are not a homogeneous Poisson process. Two caveats apply: the result depends on the window length you choose, and overdispersion only tells you that events cluster, not *why*. As Section 1.3 shows, several different mechanisms produce clustering.

Two further properties will be useful later:

- **Superposition**: combining independent Poisson processes gives a Poisson process whose rate is the sum of their rates. This underpins the branching view of Hawkes processes, in which the full process is built from a background process plus the offspring processes triggered by each event.
- **Thinning**: keeping each event independently with probability $p$ gives a Poisson process with rate $p\lambda$. More generally, keeping an event at time $t$ with probability $p(t)$ gives an inhomogeneous Poisson process (Section 1.3) with rate $p(t)\lambda$. This general form is the basis of Ogata's algorithm for simulating Hawkes processes: propose candidate events from a homogeneous process with rate $M$, an upper bound on the intensity, then keep each candidate at time $t$ with probability $\lambda^*(t)/M$, where $\lambda^*(t)$ is the intensity defined in Section 1.4.

The Poisson process is widely used across STEM, from queueing theory to reliability engineering. Its key limitation is that events cannot influence one another. In an epidemic, however, each case can cause further cases. Capturing this *self-excitation* is the core motivation for Hawkes processes.

## 1.3. Inhomogeneous Poisson process

In this case, the constant rate $\lambda$ is replaced by a function of time, $\lambda(t)$, which is non-negative and deterministic (fixed in advance, rather than random). Counts in non-overlapping intervals remain independent, but the background tendency for events rises and falls. For example, the number of flu cases varies with the seasons. The number of events in the interval $(a, b]$ is Poisson distributed,

$$
N(a, b] \sim \text{Poisson}\big(\Lambda(a, b)\big),
$$

with mean

$$
\Lambda(a, b) = \int_a^b \lambda(t)\,dt.
$$

This integrated rate reappears in Section 1.4 as the *compensator*; for an inhomogeneous Poisson process, the compensator $\Lambda(t)$ is simply $\Lambda(0, t)$.

Crucially, the rate varies because of external forces (the seasons, the weather, travel patterns), not because previous events change it. This distinction matters because a time-varying rate and self-excitation both produce clusters of events, so the two are easily confused in data. Fitting a self-exciting model with a constant background to data with an unmodelled time-varying rate will attribute the peaks to self-excitation, inflating its apparent strength. [Filimonov & Sornette (2015)](https://doi.org/10.1080/14697688.2015.1032544) demonstrated this bias in financial data. The reverse problem also occurs: if a Hawkes model is given a very flexible background, the background can absorb genuine transmission clusters, and the strength of self-excitation is underestimated.

If the rate is itself random, for example driven by unobserved fluctuations in mosquito abundance, the result is a *Cox* (or doubly stochastic) process. This is a third route to clustering, and arguably the hardest to separate from self-excitation.

In practice, these mechanisms interact. In the near-elimination setting from the introduction, the background rate of importations is driven mainly by travel, whereas the wet season acts mainly on transmission, increasing how many onward cases each case generates. Seasonality may therefore belong in the self-exciting part of the model rather than the background. We return to this when discussing time-varying baselines and kernels in Hawkes models.

## 1.4. The conditional intensity function

The most important concept in this review is the **conditional intensity**: the instantaneous expected rate of events at time $t$, given the entire history of the process up to that point.

$$
\lambda^*(t) = \lim_{h \downarrow 0} \frac{\mathbb{E}\big[\,N(t+h) - N(t) \;\big|\; \mathcal{H}(t)\,\big]}{h}
$$

The asterisk is shorthand for "conditional on the history". Breaking down the formula:

- $N(t+h) - N(t)$ is the number of events in the short window $(t, t+h]$.
- $\mathbb{E}[\,\cdots \mid \mathcal{H}(t)\,]$ is the expected number of those events, given the history $\mathcal{H}(t)$ of all events before $t$.
- Dividing by $h$ turns this expected count into a rate, and $\lim_{h \downarrow 0}$ shrinks the window to zero, giving the rate at "the next instant".

In practical terms, $\lambda^*(t)\,dt$ is approximately the probability of an event in the next instant, given everything that has happened so far.

The conditional intensity unifies the processes seen so far. For a homogeneous Poisson process, the history is irrelevant and $\lambda^{\ast}(t) = \lambda$. For an inhomogeneous Poisson process, the history is still irrelevant, but the rate varies with time: $\lambda^{\ast}(t) = \lambda(t)$. For a Hawkes process, the history matters:

$$
\lambda^*(t) = \mu(t) + \sum_{t_i < t} \phi(t - t_i).
$$

Here $\mu(t)$ is the **background rate**, generating events that arrive independently of the past (importations, in our malaria example), and $\phi$ is the **kernel**, describing how much each past event at $t_i$ raises the intensity at a time $t - t_i$ later. Each event therefore adds its own contribution to the intensity, and events generated through the kernel are the *triggered* events, or *offspring*, of earlier ones. The shape of the kernel sets the timing. An exponentially decaying kernel raises the intensity immediately after each event, whereas a kernel that starts near zero and peaks weeks later captures the delay in malaria transmission described in the introduction. The total area under the kernel,

$$
\eta = \int_0^\infty \phi(u)\,du,
$$

is the **branching ratio** $\eta$: the expected number of offspring per event, which plays the role of the reproduction number.

### Why it matters

For a simple point process, the conditional intensity completely specifies the process, and it gives the likelihood of an observed set of event times $t_1, \ldots, t_n$ on $[0, T]$ directly ([Rasmussen, 2018](https://arxiv.org/abs/1806.00221)):

$$
\log L = \sum_{i=1}^{n} \log \lambda^*(t_i) \;-\; \int_0^T \lambda^*(s)\,ds.
$$

The first term rewards the model for assigning a high intensity at the times events actually occurred. The second term is the log-probability of seeing no events in the gaps between them, including the final stretch from the last event $t_n$ to the end of observation $T$. It penalises the model for predicting events that didn't happen.

Maximising this likelihood is a common way to fit Hawkes models to data. Bayesian methods, which place priors on the parameters and sample from the posterior, and expectation–maximisation (EM) algorithms, which exploit the branching structure by estimating which events triggered which, are also widely used.

One practical complication is **edge effects**. The Hawkes intensity depends on all past events, including any that occurred before observation began at time 0. If surveillance starts partway through ongoing transmission, the earliest local cases have no recorded "parents" and may be wrongly attributed to the background.

The integral in the second term of the likelihood is known as the **compensator**:

$$
\Lambda(t) = \int_0^t \lambda^*(s)\,ds.
$$

It is the cumulative intensity up to time $t$. The difference $N(t) - \Lambda(t)$ is a *martingale*: given the history, its expected future change is zero, so on average the compensator "compensates" for the events that occur.

The compensator also provides a goodness-of-fit check, via the **time-rescaling theorem**. If the model is correct, the transformed event times $\Lambda(t_1), \Lambda(t_2), \ldots$ form a homogeneous Poisson process with rate 1. Equivalently, the rescaled gaps $\Lambda(t_i) - \Lambda(t_{i-1})$ are independent and exponentially distributed with mean 1, which can be checked with a Q–Q plot or a Kolmogorov–Smirnov test. This form of residual analysis was introduced for self-exciting models by [Ogata (1988)](https://doi.org/10.1080/01621459.1988.10478560). When the parameters have been estimated from the same data, the test is only approximate.

# 2. What is a Hawkes processes?

## 2.1. History and motivation

Hawkes processes are named after the British statistician Alan G. Hawkes, who introduced them in a pair of papers ([Hawkes, 1971a](https://academic-oup-com.bris.idm.oclc.org/biomet/article/58/1/83/224809); [Hawkes, 1971b](https://academic-oup-com.bris.idm.oclc.org/jrsssb/article/33/3/438/7027167)). His idea was to write down a point process in which events are "self-exciting", where each event temporarily raises the rate of future events. Soon after, [Hawkes and Oakes (1974)](https://www-cambridge-org.bris.idm.oclc.org/core/journals/journal-of-applied-probability/article/abs/cluster-process-representation-of-a-selfexciting-process/E836A3D07D808068E2F9F3E7E366B081) demonstrated that the same process can be viewed as a family tree, in which some events arrive spontaneously and each event goes on to produce "offspring" events of its own. This concept of branching has turned out to be the most intuitive way to think about the model.

The model soon became famous for its applications in seismology. An earthquake triggers aftershocks, which can trigger aftershocks of their own, and this idea was built into the Epidemic-Type Aftershock Sequence (ETAS) model ([Ogata, 1988](https://www-tandfonline-com.bris.idm.oclc.org/doi/abs/10.1080/01621459.1988.10478560?casa_token=W-srgC9krLYAAAAA:3HduSWc9el-inU3aIGRQcXjM6qyHG7mSpMB156f02R4JE8XDi_soQjGcxtBbyL8VQozYUJWXgp6IxQ); [Ogata, 1998](https://link-springer-com.bris.idm.oclc.org/article/10.1023/A:1003403601725)). To describe earthquakes, seisomologists borrowed the language of contagion from epidemiology. Decades later, epidemiologists would adopt this seismological model to describe contagion. 

## 2.2. Core definition 

A Hawkes process is defined by its conditional intensity. In the simplest case, with a single stream of events, it is

$$
\lambda^*(t) = \mu + \sum_{t_i < t} \phi(t - t_i).
$$

$\mu > 0$ is the background rate: the rate at which events occur spontaneously, regardless of what has occured previously. For example, in malaria, this would be the reported cases. 

$\phi(\cdot) \geq 0$ is the triggering kernel: the extra rate presently contributed by a past event. It's input, $t-t_i$, is the time elapsed since the event at $t_i$. Typically, a kernel will start high and decay, so an event's influence is strongest immediately after its occurence, and then fades. In the case of malaria, this is local transmission. 

The sum runs over every event before time $t$. Each past event adds a decaying spike to the rate, which stack on top of the background. 

This is the definition of **self-exciting**, where an event raises $\lambda^*(t)$, increasing the likelihood of another event occuring, thus raising $\lambda^*(t)$ again. This feedback is what produces clusters of events in time, and unlike the inhomogeneous Poisson process (Section 1.3), the clustering is generated by the actual events occuring, not external forces. 

## 2.3. Formalising kernels 

The kernel $\phi$ describes how an event's influence plays out overtime, and different choices of kernel will produce different model outcomes. In the context of this review, there are three worth stating. 

**The exponential kernel** is the original and most common kernel:

$$
\phi(t) = \alpha e^{-\beta t}, \qquad \alpha, \beta > 0
$$

Each event instantly raises the intensity by $\alpha$, where the increase then decays at rate $\beta$, so an event's influence lasts for $\approx 1/\beta$ units of time. The popularity of this kernal lies in its practicality - the exponential's lack of memory means that the model can be quickly fitted to data. 

**The power-law kernel** decays at a much slower rate:

$$
\phi(t) = \frac{k}{(c + t)^p}, \qquad k, c > 0, \; p > 1.
$$

Here $k$ sets the overall strength of triggering, $p$ controls how quickly an event's influence fades (larger $p$ means faster decay), and $c$ is a small offset that keeps the kernel finite at $t = 0$, so that each event raises the intensity by $k / c^p$. The condition $p > 1$ ensures that the total influence of a single event is finite.

This kernel is "heavy-tailed", so its influence fades at such a slow rate that events from long ago can still trigger new ones, whereas an exponential kernel's influence is effectively gone after a few multiples of $1/\beta$.

**A kernel can also be matched to a specific process.** $\phi$ is not necessarily a simple formula - it can be a flexible step function estimated from the data, or a shape chosen from prior knowledge. The second option is particularly important in epidemiology. 

The time between one person becoming infected and them infecting someone else is called the *generation interval*, and for many diseases its distribution has been measured and is well described by a gamma or lognormal curve. Using this curve as the kernel integrates the disease's biology directly into model - improving its efficacy.


**Quick summary:** use the exponential when speed and simplicity matter, the power law when influence persists over long periods, and a bespoke kernel when you possess data regarding the delays.

## 2.4. The branching ratio / reproduction number 

The branching ratio quantifies the number of further events which one event can trigger. As previously defined, the expected number of events produced by a rate was the area under the rate curve. An event adds $\phi$ to the rate, so the expected number of events it triggers is the area under the kernel: 

$$
n = \int_0^\infty \phi(s)\,ds.
$$

The sequential events triggered by another are described as the "offspring", in which each event has, on average, $n$ offspring. Each sequential event has its own offspring, so one spontaneous event is followed by $\approx n$ events in the first generation, $n^2$ in the second, $n^3$ in the hitds, etc. The value of $n$ determines the subsequent response:

**When $n < 1$** each generation is smaller than the previous, so each chain of events will eventually die out. The expected size of the entire cluster (including original event) is $1 + n + n^2 + \cdots = 1/1(1-n)$ ([Laub *et al.*, 2025](https://arxiv.org/pdf/1507.02822)) For example, when $n = 0.8$ each spontaneous event leads to a cluster of five events, on average. The process has a steady state that it typically returns to.

**When $n = 1$** each event will replicate itself exactly (most of the time). In this case, individual chains will still die out, but they can reach enormous sizes before doing so. Even the slightest change in $n$ will tip it into one regime (e.g., settling) or the other (exploding). For example, in epidemiology, if $R = 1$, a disease is neither disappearing nor spreading. 

**When $n = 1$** each generation of events is larger than the last, and the number of events grows without limit. 

As stated, this directly links to the example of the reproduction number $R$ in epidemiology. To make this correspondence, the kernel can be split into two 

$$
\phi(t) = R \, g(t), 
$$

where $g(t)$ is the generation interval distribution. $g$ is a probability distribution so its area is 1, meaning $n = R$. Therefore, applying the determinants above $R < 1$, $R = 1$, and $R > 1$. $R$ varying overtime produces the time varying reproduction number $R_t$. 


## 2.5. Immigrant-offspring representation
