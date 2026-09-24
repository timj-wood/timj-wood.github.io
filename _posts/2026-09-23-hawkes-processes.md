---
title: "Hawkes processes in epidemiology"
summary: A beginner friendly walkthrough of self-exciting point processes, working towards Juliette Unwin's method for separating imported and local malaria cases.
category: hawkes-processes
math: true
---
**Note:** I am an MSc Bioinformatics student, and a biologist by background; by no means a mathematician. If you spot any mistakes then please let me know!

Let's suppose you work in a country that has almost eliminated malaria, and a new case is reported. This could be from a traveller, who was infected abroad. If there are no subsequent infections, it's a dead end. However, the case might be the first visible link in a local chain where a mosquito bites the infected person and later bites someone else, and so on. These situations obviouslt require different different responses, but the day the initial case is reported they cannot be distinguished.

This is the exact problem which [Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830) aimed to address using *Hawkes processes*. To do so, they separated imported from locally acquired malaria in Yunnan (China) and Eswatini, only using the timing of cases, demoting travel histories to a simple check to their answers. Subsequently, I wanted to understand how this is possible, so here we are.

Aims: 

1. Cover enough background on random events (point processes) to understand Hawkes processes.
2. Define Hawkes processes.
3. Simulate a malaria outbreak.
4. Discuss what can occur when you apply real data.

# (Some) Background

## Events in time: point processes

To start off: a point process is a random list of event times $t_1, t_2, t_3, \ldots$ on a timeline. Point processes are observed in everywhere, for example reported disease cases, earthquakes, and neuron firings. The counting process $N(t)$ is the running tally of events up to and including time $t$, essentially a staircase that starts at 0 and sequentially increases by 1 at each event. Throughout, we'll assume that no two events can occur at the same moment (final section will cover why, in epidemiology, this is discredited). 

Everything which occurred before time $t$ is the history, $\mathcal{H}(t) = \\{t_i : t_i < t\\}$. The strict inequality is required, because whatever we predict for time $t$ may use events from the past, but not the event we are trying to predit. 

## The Poisson process

The most simplistic point process model is the homogeneous Poisson process, where events occur randomly at a constant rate $\lambda$. In any small interval of width $h$, the change of one event is roughly $\lambda h$, the chance of two or more is negligable, and the next occurrence in one interval is independent of every other interval. 

Three key points: 

- The gaps between events are exponentially distributed with mean $1/\lambda$. The process is *memoryless*: however long since the last event, the wait for the next one is the same. 
- The number of events in a window of length $T$ is $\text{Poisson}(\lambda T)$, so the variance is equal to the mean. Counts which vary more than the mean (known as overdispersion) are an early indicator that events cluster, though the cause cannot be deciphered from this. 
- If the rate changes overtime but is fixed in advance, for example $\lambda(t)$ rising around the holidays when more people travel, this produces an *inhomogeneous* Poisson process. The expected number of events between $a$ and $b$ is then the area under the rate curve, $\int_a^b \lambda(t)\,dt$.

The key limitation of the Poisson process is that events essentially ignore each other. During an epidemic, where each case causes more transmission, this is not the applicable.

<figure>
  <img src="/assets/posts/hawkes/poisson_vs_hawkes.png" alt="A Poisson process and a Hawkes process with the same background rate">
  <figcaption>Fig 1. A Poisson process and a Hawkes process with the same background rate. The Hawkes intensity jumps at each event and decays, and its events come in bursts.</figcaption>
</figure>

Fig 1 provides an illustration of the Poisson process, and alludes to the Hawkes process. Both processes have the same background rate of 0.1 events per day, but each event in a Hawkes process model causes an increase in frequency which slowly fades. In sooth, events arrive in bursts. 

## The conditional intensity 

The way around *memorylessness* is to let the rate depend on what has previously occurred. The **conditional intensity** is the expected rate of events at time $t$, given the history: 

$$
\lambda^{\ast}(t) = \lim_{h \downarrow 0} \frac{\mathbb{E}\big[\,N(t+h) - N(t) \;\big|\; \mathcal{H}(t)\,\big]}{h}.
$$

To give a brief overview of what this means: count the events in a short window after $t$, take the expected value given everything that has occurred so far, divide by the window width to turn a count into a rate, and shrink the window to zero. The asterisk is shorthand for "given the history". For small $dt$, $\lambda^{\ast}(t)\,dt$ is $\approx$ the probability of an event in the next instant. Since this function is a rate, it can exceed 1.

In the case of a Poisson process, where events are independently occurring: $\lambda^{\ast}(t) = \lambda$, or $\lambda(t)$ for inhomogeneous processes. For a Hawkes process, the history **is** accounted for. 

# The Hawkes process

## History

Hawkes processes are named after the British statistician Alan Hawkes, who introduced point processes in which each event temporarily raises the rate of future events ([Hawkes, 1971](https://doi.org/10.1093/biomet/58.1.83)). Later, [Hawkes & Oakes (1974)](https://doi.org/10.2307/3212693) demonstrated that the same process can be interpreted as a family tree, where some events arrive spontaneously, and every event can have "offspring" of its own. This view turned out to be the most intuitive way of approaching Hawkes process models. 

Initially, the model became renowned in seismology, where earthquakes trigger aftershocks that trigger aftershocks of their own, and so on. [Ogata, 1988](https://doi.org/10.1080/01621459.1988.10478560) integrated this into the Epidemic-Type Aftershock Sequence (ETAS) model. Interestingly, seismologists borrowed the language of contagion to describe earthquakes, and decades later epidemiologists would borrow a model from the seismologists to describe contagion. 

## Definition

A Hawkes process is defined by its conditional intensity:

$$
\lambda^{\ast}(t) = \mu(t) + \sum_{t_i < t} \phi(t - t_i).
$$

There are only two components:

- $\mu(t) \geq 0$ is the **background rate**: events arrive regardless of the past. In the example of epidemiology, these are imported cases. 
- $\phi(u) \geq 0$ is the **kernel**: how much a past event raises the rate $u$ days later. Only past events count, so so $\phi(u) = 0$ for $u \leq 0$. In the example, this is local transmission. 


Every previous case increases the rate, which stacks on top of the background rate. This is the precise definition of *self-exciting*, where an event raises $\lambda^{\ast}(t)$, making further events more likely, which raise $\lambda^{\ast}(t)$ again, and so on. Unlike inhomogeneous Poisson processes, where external forces were responsible for clustering, these clusters are generated by the events themselves.

## Choosing a kernel 

The kernel controls how an event's influence effects the model. Essentially, it helps to split the quantity and timing of a model: 

$$
\phi(u) = \eta \, g(u),
$$

where $\eta$ is a number and $g$ is a probability density describing the delay between a case and the cases it causes. For the context of this post, there are three different types of kernels worth noting.

**Exponential.** $\phi(u) = \alpha e^{-\beta u}$, where $\alpha, \beta > 0$. Each event immediately raises the rate by $\alpha$, and this increase fades at rate $\beta$, lasting roughly $1/\beta$. This is the most common kernel because it is quick to work with: the intensity at any moment can be updated from it value at the previous event, so fitting scales with the number of events rather than its square. Its total area is $\eta = \alpha/\beta$.

**Power law.** $\phi(u) = k/(c+u)^p$, with $k, c > 0$ and $p > 1$. This kernel decays much slower, so events which occurred long before can still have an effect. $c$ keeps the kernel finite at $u = 0$, and $p > 1$ is needed for the total area, $\eta = k / \big((p-1)c^{p-1}\big)$, to be finite. This kernel is standard for aftershocks. 

**Constructed with the biology.** The kernels above peak the instant an event occurs, which is not the case for malaria. A newly infected person must first become infectious to mosquitoes, the parasite must then develop inside the mosquito, and it must then incubate in the next person, and so on ([Huber *et al.*, 2016](https://doi.org/10.1186/s12936-016-1537-6)). Therefore, in malaria, $g$ should be close to zero for the first couple of weeks, then rise and fall. Unwin *et al.* (2021) used a kernel that is exactly zero for the first 15 days and then follows a Rayleigh (hump shaped) curve. If the event times are symptom to onset dates, $g$ is essentially the *serial interval* distribution: the time from one person's symptons to the next person's.

## The branching ratio 

As previously established, $\eta$ is the area under the kernel, 

$$
\eta = \int_0^\infty \phi(u)\,du,
$$

and it is also the **branching ratio**: the expected number of cases each case directly causes. In the epidemiology example, it is a reproduction number. Unwin *et al.* (2021) referred to this as the *case reproduction number* $R_c$, the reproduction number with whatever interventions are in place.

Let's go back to the family tree visualisation for a second. One imported case has on average $\eta$ children, $\eta^2$ grandchildren, $\eta^3$ great grandchildren, and so on: 

- **If $\eta < 1$**, each generation is smaller than the last and every chain dies out. The expected size of the whole chain, including the fist case, is $1 + \eta + \eta^2 + \cdots = 1/(1-\eta)$. With $\eta = 0.8$, each imported case leads to a chain of 5 cases on average. 
- **If $\eta = 1$**, each case replaces itself on average. Chains still die out eventually, but they can grow enormous before the ydo - this is a balanced process.
- **If $\eta > 1$**, each generation is larger than the last, and chains can grow without limitation. 

In the context of malaria, there are two consequences to this. First, with a constant background rate $\mu$ and $\eta < 1$, the process settles to an average rate of $\mu/(1-\eta)$ ([Laub *et al.*, 2015](https://arxiv.org/abs/1507.02822)). So $\eta < 1$ doesn't mean that there are zero cases, because as long as importations continue, cases continue ([Routledge *et al.*, 2018](https://doi.org/10.1038/s41467-018-04577-y)). Second, in a settled state, a randomly chosen case is locally acquired with probability $\eta$ and imported with probability $1 - \eta$. The branching ratio is therefore also the expected fraction of cases that are local. 



# Further reading

- [Laub, Taimre & Pollett (2015), *Hawkes processes*](https://arxiv.org/abs/1507.02822): a free, readable introduction to the maths, including the family-tree view. (Note that they use $\lambda$ for the background and $\mu$ for the kernel, the reverse of this post.)
- [Rasmussen (2018), *Temporal point processes and the conditional intensity function*](https://arxiv.org/abs/1806.00221): short lecture notes on the conditional intensity, the likelihood and simulation.
- [Reinhart (2018), *A review of self-exciting spatio-temporal point processes and their applications*](https://doi.org/10.1214/17-STS629): a broader review, including fitting methods and applications beyond seismology.
- [Unwin *et al.* (2021), *Using Hawkes Processes to model imported and local malaria cases in near-elimination settings*](https://doi.org/10.1371/journal.pcbi.1008830): the paper this post builds towards (open access).
