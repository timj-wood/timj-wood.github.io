---
title: "Hawkes processes in epidemiology (Part I: theory)"
category: stats
tags: hawkes
---
Let's suppose you work in a country that has almost eliminated malaria, and a new case is reported. It could be from a traveller who was infected abroad, in which case it's a dead end. Or it could be the first visible link in a local chain: a mosquito bites the infected person, later bites someone else, and so on. These situations need very different responses, but on the day the case is reported they are indiscernible.

[Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830) tackled this with *Hawkes processes*, separating imported from locally acquired malaria in Yunnan (China) and Eswatini using only the timing of cases. Travel histories were demoted to a check on the answer. I wanted to understand how that's possible, so here we are.

This post is the theory: what a Hawkes process is and why it suits an epidemic, with as little maths as I can get away with. In [part 2]({% post_url 2026-10-14-hawkes-processes-part-2 %}) we simulate a malaria outbreak and fit the model to it in R. (I'm a biologist by background, not a mathematician, so if you spot a mistake please let me know!)

## Events on a timeline

A *point process* is a random list of times at which something happens: reported disease cases, earthquakes, neuron firings. We don't care what the events are, only when they occur. At any moment, the list of events that have already happened is the **history**, and the question a point process model answers is: given the history, how likely is another event in the next little while?

## The Poisson process: events that ignore each other

The simplest model is the *Poisson process*, where events happen at a constant rate $\lambda$ (say, cases per day) and nothing that has happened so far changes that rate. Three things follow from this:

- The gaps between events are exponentially distributed, with an average gap of $1/\lambda$. The process is **memoryless**: however long it has been since the last event, the expected wait for the next one is the same.
- The number of events in a window of length $T$ is $\text{Poisson}(\lambda T)$, so the variance of the count equals its mean. If your weekly case counts vary *more* than this (overdispersion), that's an early hint that events cluster, though it doesn't tell you why.
- The rate can change over time, as long as the change is fixed in advance. Importations might rise around the holidays when more people travel, giving a rate $\lambda(t)$ that depends on the date. This is an *inhomogeneous* Poisson process, and the expected number of events in any window is just the area under the rate curve over that window.

The catch is in the first sentence: the events ignore each other. That's exactly wrong for an epidemic, where every case makes further cases more likely. Here's a Poisson process next to a Hawkes process with the same background rate of 0.1 events per day. In the Hawkes process each event bumps the rate up, the bump fades, and events arrive in bursts:

![A Poisson process and a Hawkes process with the same background rate](/images/hawkes/poisson_vs_hawkes.png)

## Letting the past matter

The way around memorylessness is to let the rate depend on what has already happened. The **conditional intensity**, written $\lambda^{\ast}(t)$, is the rate of events at time $t$ given the history up to $t$. The asterisk is shorthand for "given the history". If you multiply it by a tiny slice of time $dt$, you get roughly the probability of an event in that slice:

$$
P(\text{an event between } t \text{ and } t + dt \mid \text{history}) \approx \lambda^{\ast}(t)\,dt.
$$

It's a rate rather than a probability, so it can exceed 1. For a Poisson process, $\lambda^{\ast}(t)$ is just $\lambda$, or $\lambda(t)$ if the rate varies over time: the history makes no difference. The whole point of a Hawkes process is that it does.

## The Hawkes process

A Hawkes process is defined by its conditional intensity:

$$
\lambda^{\ast}(t) = \mu(t) + \sum_{t_i < t} \phi(t - t_i),
$$

where the sum runs over every event $t_i$ that has happened before $t$. There are only two pieces:

- $\mu(t) \geq 0$ is the **background rate**: events that arrive regardless of the past. For us, imported cases.
- $\phi(u) \geq 0$ is the **kernel**: how much a past event raises the rate $u$ days later. For us, local transmission. It's zero for $u \leq 0$, because a case can't cause cases before it happens.

So at any moment, the intensity is the background plus one bump for every previous case, each bump being however much that case is still "pushing" given how long ago it was. That's what **self-exciting** means: an event raises $\lambda^{\ast}(t)$, making more events likely, which raise $\lambda^{\ast}(t)$ again, and so on. In an inhomogeneous Poisson process, clustering comes from outside (the seasons, the holidays); here the events generate it themselves.

A bit of history: Alan Hawkes introduced these in [1971](https://doi.org/10.1093/biomet/58.1.83), and [Hawkes & Oakes (1974)](https://doi.org/10.2307/3212693) showed that the same process can be read as a family tree, which we'll lean on below. The model made its name in seismology, where aftershocks trigger aftershocks, through [Ogata's ETAS model](https://doi.org/10.1080/01621459.1988.10478560). So seismologists borrowed the language of contagion for earthquakes, and epidemiologists later borrowed their model back.

## Choosing a kernel

It helps to split the kernel into *how many* and *when*:

$$
\phi(u) = \eta \, g(u).
$$

Here $\eta$ is a single number, the total size of the bump, and $g(u)$ is a probability distribution for the delay between a case and the cases it causes. Three kernels are worth knowing about.

**Exponential.** $\phi(u) = \kappa e^{-\beta u}$. Each event immediately raises the rate by $\kappa$, and the bump fades away over roughly $1/\beta$ days. It's the most common choice because it's cheap to compute: the intensity at any moment can be updated from its value at the previous event, so fitting scales with the number of events rather than its square.

**Power law.** Something like $\phi(u) = k/(c+u)^p$. This decays much more slowly, so old events still matter long after they happen. It's the standard for aftershocks, where a big earthquake can keep triggering events for years.

**Built from the biology.** Both kernels above peak the instant an event happens, which is wrong for malaria. A newly infected person has to become infectious to mosquitoes, the parasite has to develop inside the mosquito, and then it has to incubate in the next person ([Huber *et al.*, 2016](https://doi.org/10.1186/s12936-016-1537-6)). So $g$ should be near zero for the first couple of weeks, then rise and fall. Unwin *et al.* used a kernel that is exactly zero for 15 days and then follows a hump-shaped curve (a Rayleigh distribution, chosen for convenience). If the event times are symptom onset dates, $g$ is essentially the *serial interval*: the distribution of the delay between one case's symptoms and the next's.

## The branching ratio

The area under the kernel,

$$
\eta = \int_0^\infty \phi(u)\,du,
$$

is called the **branching ratio**, and it has a direct meaning: the average number of cases that each case directly causes. In other words, it's a reproduction number. The paper calls it the *case reproduction number*, $R_c$, the reproduction number with current interventions in place.

Now the family tree. One imported case has on average $\eta$ children, $\eta^2$ grandchildren, $\eta^3$ great-grandchildren, and so on:

- **If $\eta < 1$**, each generation is smaller than the last and every chain dies out. The expected chain size, including the first case, is $1 + \eta + \eta^2 + \cdots = 1/(1-\eta)$. With $\eta = 0.8$, each imported case leads to 5 cases on average.
- **If $\eta = 1$**, each case replaces itself on average. Chains still die out eventually, but they can get enormous first.
- **If $\eta > 1$**, each generation is bigger than the last, and chains can grow without limit.

Two consequences for malaria. First, with a constant background $\mu$ and $\eta < 1$, the process settles down to an average rate of $\mu/(1-\eta)$ cases per day. So $\eta < 1$ doesn't mean zero cases: as long as importations continue, so do local cases ([Routledge *et al.*, 2018](https://doi.org/10.1038/s41467-018-04577-y)). Second, in that settled state, a randomly chosen case is local with probability $\eta$ and imported with probability $1 - \eta$. The branching ratio is also the expected *fraction* of cases that are local, which is exactly the number a malaria programme wants.

## Imported or local, in principle

This was the most enticing part of the paper, and the idea is simple once you have the model. At the moment case $i$ occurs, the intensity is a sum of named contributions:

$$
\lambda^{\ast}(t_i) = \underbrace{\mu(t_i)}_{\text{importation}} + \underbrace{\phi(t_i - t_1)}_{\text{case 1}} + \underbrace{\phi(t_i - t_2)}_{\text{case 2}} + \cdots
$$

Each term is a candidate explanation for case $i$, and its size is how strongly it was pushing at that moment. So the probability that case $i$ was imported is the background's share of the total,

$$
P(\text{case } i \text{ imported}) = \frac{\mu(t_i)}{\lambda^{\ast}(t_i)},
$$

and the probability that case $j$ infected case $i$ is *its* share, $\phi(t_i - t_j)/\lambda^{\ast}(t_i)$. A case in the off-season with nothing in the preceding five weeks gets a probability of importation near 1. A case at the peak of the season, five weeks after a cluster, is ambiguous, and the model says so.

That's the whole trick: fit $\mu(t)$ and $\phi(u)$ to the case times, then read off each case's probability of being imported. Whether it actually works, and what the case times can and can't tell you, is the subject of [part 2]({% post_url 2026-10-14-hawkes-processes-part-2 %}), where we simulate an outbreak with a known answer and see what the model recovers.

## Further reading

- [Laub, Taimre & Pollett (2015)](https://arxiv.org/abs/1507.02822): a free, readable introduction to the maths, including the family-tree view. (They use $\lambda$ for the background and $\mu$ for the kernel, the reverse of this post.)
- [Rasmussen (2018)](https://arxiv.org/abs/1806.00221): short lecture notes on the conditional intensity, the likelihood and simulation.
- [Reinhart (2018)](https://doi.org/10.1214/17-STS629): a broader review, including fitting methods and applications beyond seismology.
- [Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830): the paper this post builds towards (open access), and [epihawkes](https://github.com/mrc-ide/epihawkes), the R package released with it.

If you spot any mistakes or have questions, feel free to email me. Thanks for reading!
