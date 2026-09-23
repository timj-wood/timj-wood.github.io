---
title: "Hawkes processes in epidemiology"
summary: A beginner-friendly walkthrough of self-exciting point processes, working towards Juliette Unwin's method for separating imported and local malaria cases.
category: hawkes-processes
math: true
---

Suppose you work in a country that has almost eliminated malaria, and a new case is reported. It might be a traveller infected abroad, a dead end if nothing follows. Or it might be the first visible link in a local chain, where a mosquito bites the infected person and later bites someone else. These situations need very different responses, yet on the day the case is reported they can look identical. Often the only evidence is a list of case dates and some incomplete travel histories.

I'm an MSc bioinformatics student, and a biologist by background; by no means a mathematician. I came across [Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830), who aimed to address this exact problem using something called a *Hawkes process*: they separated imported from locally acquired malaria in Yunnan (China) and Eswatini, only using the timing of cases, and kept the travel histories back to check their answer. I wanted to understand how that is possible, and this post is me working through it from scratch.

The plan is to provide:

1. Enough background on random events for understanding. 
2. A definition for a Hawkes process itself.
3. A simulated malaria outbreak.
4. Information of what changes when the data are real.

If you spot any mistakes (there are likely many), please let me know.

# Some background

## Events in time

A *point process* is a random list of event times $t_1, t_2, t_3, \ldots$ on a timeline, these are observed in many phenomena such as reported cases, earthquake aftershocks, and neurons firing. Its *counting process* $N(t)$ is the running tally of events up to and including time $t$, a staircase that starts at 0 and steps up by 1 at each event. Throughout, I'll assume no two events happen at exactly the same moment (the last section will cover why real surveillance data break this assumption).

Everything that happened before time $t$ is the *history*, written $\mathcal{H}(t) = \{t_i : t_i < t\}$. The strict inequality is required, because whatever we predict for time $t$ may use the past, but not the event we are trying to predict.

## The Poisson process

The simplest point process model is the homogeneous Poisson process, where events happen completely at random at a constant rate $\lambda$. In any tiny interval of width $h$, the chance of one event is about $\lambda h$, the chance of two or more is negligible, and what happens in one interval is independent of every other interval.

Three consequences are worth knowing:

- The gaps between events are exponentially distributed with mean $1/\lambda$. The process has no memory: however long since the last event, the wait for the next is the same.
- The number of events in a window of length $T$ is $\text{Poisson}(\lambda T)$, so its variance equals its mean. Counts that vary *more* than their mean (overdispersion) are a first hint that events cluster, although not a hint about why.
- If the rate changes over time but is fixed in advance, say $\lambda(t)$ rising around holidays when more people travel, we get an *inhomogeneous* Poisson process. The expected number of events between $a$ and $b$ is then the area under the rate curve, $\int_a^b \lambda(t)\,dt$.

The key limitation is in the section title: events *ignore each other*. In an epidemic, that's exactly wrong, because each case can cause more cases.

<figure>
  <img src="/assets/posts/hawkes/poisson_vs_hawkes.png" alt="A Poisson process and a Hawkes process with the same background rate">
  <figcaption>Fig 1. A Poisson process and a Hawkes process with the same background rate. The Hawkes intensity jumps at each event and decays, and its events come in bursts.</figcaption>
</figure>

The figure above previews where we're heading. Both processes have the same background rate of 0.1 events per day, but in the Hawkes process each event gives the rate a jump that then fades away, so events arrive in bursts.

## The conditional intensity

The way out is to let the rate depend on what has already happened. The **conditional intensity** is the expected rate of events at time $t$, given the history:

$$
\lambda^{\ast}(t) = \lim_{h \downarrow 0} \frac{\mathbb{E}\big[\,N(t+h) - N(t) \;\big|\; \mathcal{H}(t)\,\big]}{h}.
$$

Read it from the inside out: count the events in a short window after $t$, take the expected value given everything so far, divide by the window width to turn a count into a rate, and shrink the window to zero. The asterisk is shorthand for "given the history". For small $dt$, $\lambda^{\ast}(t)\,dt$ is roughly the probability of an event in the next instant. (It is a rate, such as cases per day, not a probability, so it can exceed 1.)

For a Poisson process, the history doesn't matter: $\lambda^{\ast}(t) = \lambda$, or $\lambda(t)$ in the inhomogeneous case. For a Hawkes process, it does.

# The Hawkes process

## Where it comes from

Hawkes processes are named after the statistician Alan Hawkes, who introduced point processes in which each event temporarily raises the rate of future events ([Hawkes, 1971](https://doi.org/10.1093/biomet/58.1.83)). A few years later, [Hawkes & Oakes (1974)](https://doi.org/10.2307/3212693) showed that the same process can be viewed as a family tree: some events arrive spontaneously, and every event can have "offspring" of its own. That family-tree view turns out to be the most intuitive way to think about the model.

The model became famous in seismology, where earthquakes trigger aftershocks that trigger aftershocks of their own. Ogata built this into the Epidemic-Type Aftershock Sequence (ETAS) model ([Ogata, 1988](https://doi.org/10.1080/01621459.1988.10478560)). Seismologists borrowed the language of contagion to describe earthquakes, and decades later epidemiologists borrowed the seismologists' model to describe contagion.

## The definition

A Hawkes process is defined by its conditional intensity:

$$
\lambda^{\ast}(t) = \mu(t) + \sum_{t_i < t} \phi(t - t_i).
$$

There are only two components:

- $\mu(t) \geq 0$ is the **background rate**: events that arrive regardless of the past. In our example, these are imported cases.
- $\phi(u) \geq 0$ is the **kernel**: how much a past event raises the rate $u$ days later. Only past events count, so $\phi(u) = 0$ for $u \leq 0$. In our example, this is local transmission.

Every past case adds its own bump to the rate, and the bumps stack on top of the background. That is what *self-exciting* means: an event raises $\lambda^{\ast}(t)$, making further events more likely, which raise $\lambda^{\ast}(t)$ again. Unlike the holiday-driven clusters of an inhomogeneous Poisson process, these clusters are generated by the events themselves.

## Choosing a kernel

The kernel controls how an event's influence plays out over time. It helps to split it into "how many" and "when":

$$
\phi(u) = \eta \, g(u),
$$

where $\eta$ is a number (more on it below) and $g$ is a probability density describing the delay between a case and the cases it causes. Three choices come up again and again.

**Exponential.** $\phi(u) = \alpha e^{-\beta u}$, with $\alpha, \beta > 0$. Each event immediately raises the rate by $\alpha$, and the bump fades at rate $\beta$, lasting roughly $1/\beta$. It is the most common kernel because it is fast to work with: the intensity at any moment can be updated from its value at the previous event, so fitting scales with the number of events rather than its square. Its total area is $\eta = \alpha/\beta$.

**Power law.** $\phi(u) = k/(c+u)^p$, with $k, c > 0$ and $p > 1$. This decays much more slowly, so events from long ago can still matter. $c$ keeps the kernel finite at $u = 0$, and $p > 1$ is needed for the total area, $\eta = k / \big((p-1)c^{p-1}\big)$, to be finite. It is standard for aftershocks.

**Built from biology.** Both kernels above peak the instant an event happens, which is wrong for malaria. A newly infected person must first become infectious to mosquitoes, the parasite must develop inside the mosquito, and it must then incubate in the next person ([Huber *et al.*, 2016](https://doi.org/10.1186/s12936-016-1537-6)). So in malaria, $g$ should be close to zero for the first couple of weeks, then rise and fall. Unwin *et al.* used a kernel that is exactly zero for the first 15 days and then follows a hump-shaped (Rayleigh) curve. If the event times are symptom-onset dates, $g$ is essentially the *serial interval* distribution: the time from one person's symptoms to the next person's.

## The branching ratio

The number $\eta$ is the area under the kernel,

$$
\eta = \int_0^\infty \phi(u)\,du,
$$

and it is the **branching ratio**: the expected number of cases each case directly causes. In epidemiological terms, it is a reproduction number. Unwin *et al.* call it the *case reproduction number* $R_c$, the reproduction number with whatever interventions are in place.

Think of the family tree again. One imported case has on average $\eta$ children, $\eta^2$ grandchildren, $\eta^3$ great-grandchildren, and so on:

- **If $\eta < 1$**, each generation is smaller than the last and every chain dies out. The expected size of a whole chain, including the first case, is $1 + \eta + \eta^2 + \cdots = 1/(1-\eta)$. With $\eta = 0.8$, each imported case leads to a chain of 5 cases on average.
- **If $\eta = 1$**, each case replaces itself on average. Chains still die out eventually, but they can grow very large first; the process is balanced on a knife-edge.
- **If $\eta > 1$**, each generation is larger than the last, and chains can grow without limit.

Two consequences matter for the malaria question. First, with a constant background rate $\mu$ and $\eta < 1$, the process settles to an average rate of $\mu/(1-\eta)$ ([Laub *et al.*, 2015](https://arxiv.org/abs/1507.02822)). So $\eta < 1$ does *not* mean zero cases: as long as importations continue, cases continue ([Routledge *et al.*, 2018](https://doi.org/10.1038/s41467-018-04577-y)). Second, in that settled state a randomly chosen case is locally acquired with probability $\eta$ and imported with probability $1 - \eta$. The branching ratio is therefore also the expected *fraction* of cases that are local.

# Simulating a malaria outbreak

The best way I found to understand all this is to simulate an outbreak where I know the truth, and then see whether the model can recover it. All the code for this section is in a single Python script, [hawkes_malaria.py](/assets/posts/hawkes/hawkes_malaria.py), which uses only NumPy, SciPy and Matplotlib.

## Setting it up

For illustration (these are made-up numbers, not estimates), I simulated five years of cases in which:

- imported cases arrive at a seasonal background rate, $\mu(t) = m\,\big(1 + a\cos(2\pi(t - t_{\text{peak}})/365)\big)$, averaging $m = 0.1$ per day (about 3 per month), with a strong annual cycle ($a = 0.8$) peaking around day 200 of each year;
- each case causes on average $\eta = 0.6$ local cases;
- the delay kernel $g$ is zero for the first 15 days and then follows a Rayleigh curve with spread $\sigma = 20$ days, so onward cases are most likely about five weeks later.

Because $\eta = 0.6$, we expect about 60% of cases to be local, and every chain eventually dies out.

## Simulating it

There are two natural ways to simulate a Hawkes process.

**Family tree.** First draw the imported cases from the background rate. Then, for each case, draw its number of children from a Poisson distribution with mean $\eta$, and give each child a delay drawn from $g$. Repeat for the children, the grandchildren, and so on, until a generation has no children. This is the Hawkes & Oakes picture turned directly into an algorithm, and it has a bonus: we know exactly which cases are imported and who infected whom.

**Thinning.** [Ogata (1981)](https://doi.org/10.1109/TIT.1981.1056305) gave the standard general-purpose method. Propose candidate events at a rate $m$ that is an upper bound on $\lambda^{\ast}(s)$ over a short window ahead, given the events so far, and accept a candidate at time $s$ with probability $\lambda^{\ast}(s)/m$. After each step, update the history and the bound. For kernels that only decay, the current intensity is a valid bound until the next event, but for a delayed kernel like ours it is not, which is why Unwin *et al.* had to design their own version.

I used the family-tree method for the outbreak (and thinning for the exponential-kernel example in the first figure). The run shown here produced 426 cases, 174 imported and 252 local, so 59% local, close to the expected 60%.

<figure>
  <img src="/assets/posts/hawkes/simulated_outbreak.png" alt="Simulated outbreak split into imported and local cases">
  <figcaption>Fig 2. The simulated outbreak. The total intensity is split into the seasonal background (imported cases) and the triggered part (local transmission); below, each case is marked by its true origin.</figcaption>
</figure>

Notice how the local cases lag behind each seasonal peak of importations, and how the triggered part of the intensity is often larger than the background itself.

## Fitting it

Now pretend we only have the case dates. The conditional intensity gives us the likelihood of the observed times $t_1, \ldots, t_n$ on $[0, T]$ directly ([Rasmussen, 2018](https://arxiv.org/abs/1806.00221)):

$$
\log L = \sum_{i=1}^{n} \log \lambda^{\ast}(t_i) \;-\; \int_0^T \lambda^{\ast}(s)\,ds.
$$

The first term rewards the model for putting a high rate where cases actually happened. The second is the log-probability of seeing no cases in all the gaps, including the stretch after the last case, and penalises the model for predicting cases that didn't happen.

Maximising this over the parameters is how Unwin *et al.* fitted their model. The likelihood can have several peaks, so it is worth repeating the fit from different starting points. The alternatives are Bayesian methods and the EM algorithm, which comes next.

## Telling imported from local

This is the payoff. Once we have a fitted model, look at the intensity at the moment a case $i$ occurred. It is a sum of pieces: the background, plus one bump from each earlier case. Each piece's share of the total is the probability that it "caused" case $i$:

$$
P(\text{case } i \text{ is imported}) = \frac{\mu(t_i)}{\lambda^{\ast}(t_i)},
$$

$$
P(\text{case } i \text{ was caused by case } j) = \frac{\phi(t_i - t_j)}{\lambda^{\ast}(t_i)}.
$$

These probabilities add up to 1 for each case. They are exactly what the EM algorithm estimates on each pass: treat the unknown family tree as missing data, guess it given the current parameters, re-fit the parameters given the guess, and repeat. In seismology, labelling events this way is called *stochastic declustering* ([Zhuang, Ogata & Vere-Jones, 2002](https://doi.org/10.1198/016214502760046925)). For us, it is an imported-versus-local classifier that uses timing alone.

In the toy outbreak we know the true labels, so we can check how well it works, just as Unwin *et al.* checked theirs against travel histories.

Fitting the model by maximum likelihood (with the 15-day delay fixed, as Unwin *et al.* did) recovered most parameters well:

| Parameter | True | Fitted |
|---|---|---|
| Average import rate $m$ (per day) | 0.100 | 0.102 |
| Seasonal strength $a$ | 0.80 | 0.75 |
| Seasonal peak (day of year) | 200 | 203 |
| Branching ratio $\eta$ | 0.60 | 0.57 |
| Delay spread $\sigma$ (days) | 20 | 14 |

The branching ratio and the import pattern come back close to the truth. The spread of the delay is less well pinned down, which makes sense: with only a few hundred cases, many different delay shapes explain the data almost equally well.

<figure>
  <img src="/assets/posts/hawkes/classification.png" alt="Estimated probability of being imported, by true origin">
  <figcaption>Fig 3. Each case's estimated probability of being imported, split by its true origin.</figcaption>
</figure>

The classification results were the most interesting part for me. Adding up the probabilities gives an expected 185 imported cases, against a true 174, so the model gets the overall split roughly right. Individual cases are much harder: labelling each case as imported when its probability is above 0.5 gets 68% right, and the two groups overlap a lot in the figure. That's because an imported case arriving during a burst of local transmission looks, in timing alone, just like a local case. Timing tells you a lot about *how many* cases are imported, but much less about *which* ones, and this is exactly where partial travel histories become valuable.

## Checking the fit

The integral in the likelihood, $\Lambda(t) = \int_0^t \lambda^{\ast}(s)\,ds$, gives a neat goodness-of-fit check. If the model is right, the gaps $\Lambda(t_i) - \Lambda(t_{i-1})$ should look like independent draws from an exponential distribution with mean 1: stretching time by the fitted intensity turns the data into a plain Poisson process. You can check this with a Q–Q plot or a Kolmogorov–Smirnov test, which is the *residual analysis* Ogata (1988) introduced for earthquake models and Unwin *et al.* used for malaria. When the parameters were fitted to the same data, the test is lenient, so passing it is encouraging rather than conclusive.

<figure class="small">
  <img src="/assets/posts/hawkes/time_rescaling_qq.png" alt="Q–Q plot of rescaled gaps against an exponential distribution">
  <figcaption>Fig 4. Q–Q plot of the rescaled gaps against an exponential distribution with mean 1. The points lie close to the diagonal.</figcaption>
</figure>

For the toy outbreak, the points lie close to the diagonal and the Kolmogorov–Smirnov test finds no evidence against the model ($p = 0.70$), as it should, since the model is correct by construction.

# What changes with real data

A simulation is the best case: the model is correct by construction. Real surveillance data raise problems at every step, and most of the hard work in Unwin *et al.*'s paper is in handling them.

**Tied dates.** Cases are recorded to the day, so several often share a date, which breaks the one-event-at-a-time assumption. Unwin *et al.* spread tied cases randomly within the day ("jittering"). That is reasonable when the kernel changes little within a day, as with malaria's multi-week delay, but with a kernel that peaks immediately, same-day cases would wrongly appear to cause each other. The alternative is to model daily counts directly with a discrete-time version of the model.

**Which date?** The recorded date is usually symptom onset, diagnosis or reporting, not infection. The observed timeline is shifted relative to transmission and, because delays vary between people, blurred. With onset dates, the kernel describes the serial interval rather than the time between infections ([Huber *et al.*, 2016](https://doi.org/10.1186/s12936-016-1537-6)).

**Other routes to clustering.** Holiday travel produces clusters of imported cases that look, in a list of dates, much like a burst of transmission. If the background is assumed constant when it really varies, the model blames the peaks on transmission: [Filimonov & Sornette (2015)](https://doi.org/10.1080/14697688.2015.1032544) showed that a Hawkes fit to simulated data with shifting rates, and no self-excitation at all, can report a branching ratio close to 1. The reverse can happen too: a very flexible background can soak up genuine transmission clusters. Hidden drivers cause similar trouble. An unobserved reservoir of asymptomatic infections produces reported cases that seem to have no parent, and look imported. Even with the right model, labelling individual events as background or triggered is harder than it looks ([Sornette & Utkin, 2009](https://doi.org/10.1103/PhysRevE.79.061110)).

**Seasonality.** Seasonality can act on importations (seasonal travel, and seasonal transmission where travellers come from) and on local transmission (more mosquitoes in the wet season means more onward cases, and even changes their timing). Unwin *et al.* put it in the background rate, as an annual cycle plus a long-term trend. Putting it in the kernel instead, as a reproduction number that varies through the year, is also defensible. Which is right depends on the setting.

**Edge effects.** The intensity depends on all past cases, including any before surveillance began, but the likelihood assumes there were none. So if data start partway through ongoing transmission, the earliest local cases have no recorded parents and look imported. With a multi-week kernel, this affects the first weeks to months of data. A common fix is to use an initial period as a burn-in: its cases feed into the intensity but not the likelihood ([Reinhart, 2018](https://doi.org/10.1214/17-STS629)). At the other end, cases near the end of the data have children that haven't happened yet, so naive counts of offspring per case are too low.

**Under-reporting and the reproduction number.** $\eta$ counts *reported* offspring, so missed cases pull it down. Unwin *et al.* explored this by randomly deleting simulated cases. Reading $\eta$ as a reproduction number also assumes a large pool of susceptible people, which is reasonable at near-elimination but not in a large epidemic.

# Further reading

- [Laub, Taimre & Pollett (2015), *Hawkes processes*](https://arxiv.org/abs/1507.02822): a free, readable introduction to the maths, including the family-tree view. (Note that they use $\lambda$ for the background and $\mu$ for the kernel, the reverse of this post.)
- [Rasmussen (2018), *Temporal point processes and the conditional intensity function*](https://arxiv.org/abs/1806.00221): short lecture notes on the conditional intensity, the likelihood and simulation.
- [Reinhart (2018), *A review of self-exciting spatio-temporal point processes and their applications*](https://doi.org/10.1214/17-STS629): a broader review, including fitting methods and applications beyond seismology.
- [Unwin *et al.* (2021), *Using Hawkes Processes to model imported and local malaria cases in near-elimination settings*](https://doi.org/10.1371/journal.pcbi.1008830): the paper this post builds towards (open access).
