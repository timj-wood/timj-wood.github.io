---
title: "Hawkes processes in epidemiology"
summary: A beginner friendly walkthrough of self-exciting point processes, working towards Juliette Unwin's method for separating imported and local malaria cases.
category: hawkes
math: true
---
Let's suppose you work in a country that has almost eliminated malaria, and a new case is reported. It could be from a traveller who was infected abroad, in which case it's a dead end. Or it could be the first visible link in a local chain: a mosquito bites the infected person, later bites someone else, and so on. These situations need very different responses, but on the day the case is reported they look identical.

[Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830) tackled this with *Hawkes processes*, separating imported from locally acquired malaria in Yunnan (China) and Eswatini using only the timing of cases. Travel histories were demoted to a check on the answer. I wanted to understand how that's possible, so here we are.

We'll cover just enough about point processes to define the Hawkes process, then simulate a malaria outbreak and see what the model can recover from it (I'm a biologist by background, not a mathematician, so if you spot a mistake please let me know!).

## Point processes

A point process is a random list of event times $t_1, t_2, t_3, \ldots$ on a timeline: reported disease cases, earthquakes, neuron firings. The counting process $N(t)$ is the running tally of events up to and including time $t$, a staircase that starts at 0 and steps up by 1 at each event. We'll assume no two events happen at exactly the same moment.

Everything before time $t$ is the history, $\mathcal{H}(t) = \\{t_i : t_i < t\\}$. The inequality is strict because a prediction for time $t$ can use the past, but not the event we're trying to predict.

## The Poisson process

The simplest model is the homogeneous Poisson process, where events occur at a constant rate $\lambda$. In any small interval of width $h$, the chance of one event is roughly $\lambda h$, the chance of two or more is negligible, and what happens in one interval is independent of every other. This gives us three useful facts:

- The gaps between events are exponentially distributed with mean $1/\lambda$. The process is *memoryless*: however long since the last event, the wait for the next one is the same.
- The number of events in a window of length $T$ is $\text{Poisson}(\lambda T)$, so the variance equals the mean. Counts that vary more than this (overdispersion) are an early hint that events cluster, though not of why.
- If the rate changes over time but is fixed in advance, for example $\lambda(t)$ rising around the holidays when more people travel, we get an *inhomogeneous* Poisson process. The expected number of events between $a$ and $b$ is the area under the rate curve, $\int_a^b \lambda(t)\,dt$.

The catch is that events ignore each other. That's exactly wrong for an epidemic, where each case causes more transmission. Here's a Poisson process next to a Hawkes process with the same background rate of 0.1 events per day. In the Hawkes process each event bumps the rate up, the bump fades, and events arrive in bursts:

![A Poisson process and a Hawkes process with the same background rate](/images/hawkes/poisson_vs_hawkes.png)

## The conditional intensity

The way around memorylessness is to let the rate depend on what has already happened. The **conditional intensity** is the expected rate of events at time $t$, given the history:

$$
\lambda^{\ast}(t) = \lim_{h \downarrow 0} \frac{\mathbb{E}\big[\,N(t+h) - N(t) \;\big|\; \mathcal{H}(t)\,\big]}{h}.
$$

In words: count the events in a short window after $t$, take the expectation given everything so far, divide by the window width to turn a count into a rate, and shrink the window to zero. The asterisk is shorthand for "given the history". For small $dt$, $\lambda^{\ast}(t)\,dt$ is roughly the probability of an event in the next instant (it's a rate, so it can exceed 1).

For a Poisson process, $\lambda^{\ast}(t)$ is just $\lambda$, or $\lambda(t)$ in the inhomogeneous case. The history doesn't come into it. For a Hawkes process, it does.

## The Hawkes process

A Hawkes process is defined by its conditional intensity:

$$
\lambda^{\ast}(t) = \mu(t) + \sum_{t_i < t} \phi(t - t_i).
$$

There are only two pieces:

- $\mu(t) \geq 0$ is the **background rate**: events that arrive regardless of the past. For us, imported cases.
- $\phi(u) \geq 0$ is the **kernel**: how much a past event raises the rate $u$ days later, with $\phi(u) = 0$ for $u \leq 0$. For us, local transmission.

Every previous case adds its own bump on top of the background. That's what *self-exciting* means: an event raises $\lambda^{\ast}(t)$, making more events likely, which raise $\lambda^{\ast}(t)$ again, and so on. In an inhomogeneous Poisson process clustering comes from outside; here the events generate it themselves.

(A bit of history: Alan Hawkes introduced these in [1971](https://doi.org/10.1093/biomet/58.1.83), and [Hawkes & Oakes (1974)](https://doi.org/10.2307/3212693) showed the same process can be read as a family tree, which we'll lean on below. The model made its name in seismology, where aftershocks trigger aftershocks, through [Ogata's ETAS model](https://doi.org/10.1080/01621459.1988.10478560). So seismologists borrowed the language of contagion for earthquakes, and epidemiologists later borrowed their model back.)

## Choosing a kernel

It helps to split the kernel into *how many* and *when*:

$$
\phi(u) = \eta \, g(u),
$$

where $\eta$ is a number and $g$ is a probability density for the delay between a case and the cases it causes. Three kernels are worth knowing about here.

**Exponential.** $\phi(u) = \kappa e^{-\beta u}$, with $\kappa, \beta > 0$. Each event immediately raises the rate by $\kappa$, and the bump fades at rate $\beta$, lasting roughly $1/\beta$. It's the most common choice because it's cheap: the intensity at any moment can be updated from its value at the previous event, so fitting scales with the number of events rather than its square. Its total area is $\eta = \kappa/\beta$.

**Power law.** $\phi(u) = k/(c+u)^p$, with $k, c > 0$ and $p > 1$. This decays much more slowly, so old events still matter. $c$ keeps the kernel finite at $u = 0$, and $p > 1$ keeps the total area, $\eta = k / \big((p-1)c^{p-1}\big)$, finite. This is the standard for aftershocks.

**Built from the biology.** Both kernels above peak the instant an event happens, which is wrong for malaria. A newly infected person has to become infectious to mosquitoes, the parasite has to develop inside the mosquito, and then it has to incubate in the next person ([Huber *et al.*, 2016](https://doi.org/10.1186/s12936-016-1537-6)). So $g$ should be near zero for the first couple of weeks, then rise and fall. Unwin *et al.* used a kernel that is exactly zero for 15 days and then follows a hump-shaped Rayleigh curve. If the event times are symptom onset dates, $g$ is essentially the *serial interval* distribution.

## The branching ratio

The area under the kernel,

$$
\eta = \int_0^\infty \phi(u)\,du,
$$

is the **branching ratio**: the expected number of cases each case directly causes. In other words, a reproduction number (the paper calls it the *case reproduction number* $R_c$, the reproduction number with current interventions in place).

Back to the family tree. One imported case has on average $\eta$ children, $\eta^2$ grandchildren, $\eta^3$ great grandchildren, and so on:

- **If $\eta < 1$**, each generation is smaller than the last and every chain dies out. The expected chain size, including the first case, is $1 + \eta + \eta^2 + \cdots = 1/(1-\eta)$. With $\eta = 0.8$, each imported case leads to 5 cases on average.
- **If $\eta = 1$**, each case replaces itself on average. Chains still die out eventually, but they can get enormous first.
- **If $\eta > 1$**, each generation is bigger than the last, and chains can grow without limit.

Two consequences for malaria. First, with a constant background $\mu$ and $\eta < 1$, the process settles to an average rate of $\mu/(1-\eta)$ ([Laub *et al.*, 2015](https://arxiv.org/abs/1507.02822)), so $\eta < 1$ doesn't mean zero cases: as long as importations continue, so do cases ([Routledge *et al.*, 2018](https://doi.org/10.1038/s41467-018-04577-y)). Second, once settled, a randomly chosen case is local with probability $\eta$ and imported with probability $1 - \eta$. The branching ratio is also the expected fraction of cases that are local.

## Simulating a malaria outbreak

The best method which I have found in understanding this is to simulate an outbreak where the true transmission tree is known, then fit the model and see what it recovers. I used the simulator, likelihood, gradient and intensity plot from Unwin, Mishra and Bhatt's excellent [epihawkes](https://github.com/mrc-ide/epihawkes) package, and my scripts are at [hawkes-malaria](https://github.com/timj-wood/hawkes-malaria). The package isn't on CRAN, so install it from GitHub:

```r
remotes::install_github("mrc-ide/epihawkes")
library(epihawkes)
```

The setup (numbers made up for illustration) is five years of cases where:

- Imported cases arrive at a seasonal background rate, $\mu(t) = m\,\big(1 + a\cos(2\pi(t - t_{\text{peak}})/365.25)\big)$, averaging $m = 0.1$ per day ($\approx 3$ per month), with a strong annual cycle ($a = 0.8$) peaking around day 200.
- Each case causes, on average, $\eta = 0.6$ local cases.
- The delay kernel $g$ is zero for 15 days, then follows a Rayleigh curve with spread $\sigma = 20$ days, so onward cases are most likely about five weeks later.

Since $\eta = 0.6$, about 60% of cases should be local, and every chain should eventually die out.

`epihawkes` uses a different parameterisation, so first we translate. Its Rayleigh kernel, `ray_kernel()`, is

$$
\phi(u) = \alpha\,(u - d)\,e^{-\frac{1}{2}\delta (u - d)^2} \quad \text{for } u > d,
$$

where $d$ is the 15-day delay. This is $\eta\, g(u)$ with $\sigma^2 = 1/\delta$ and $\eta = \alpha/\delta$. Its seasonal background, `mu_sinusoidal_constant()`, is $A + M\cos(\omega t) + N\sin(\omega t)$ with $\omega = 2\pi/365.25$, which is my $\mu(t)$ with $A = m$ and $(M, N) = m a\,(\cos \omega t_{\text{peak}}, \sin \omega t_{\text{peak}})$:

```r
T_max <- 5 * 365
m <- 0.10; a <- 0.80; t_peak <- 200; eta <- 0.60; sigma <- 20; delay <- 15
w <- 2 * pi / 365.25

true <- list(alpha = eta / sigma^2, delta = 1 / sigma^2, delay = delay,
             A = m, M = m * a * cos(w * t_peak), N = m * a * sin(w * t_peak))
```

I'll keep talking about $\eta$ and $\sigma$ since they mean something directly (a reproduction number and a spread in days); $\alpha$ and $\delta$ are just what the optimiser sees.

`hawkes_simulation()` uses Ogata's thinning algorithm: walk forward in time, bound the intensity over the near future, propose the next event at that bound, and accept it with probability $\lambda^{\ast}(t)/\text{bound}$. It works for any kernel and background because it only ever evaluates $\lambda^{\ast}(t)$.

```r
events_thin <- hawkes_simulation(events = c(0), T_max = T_max, kernel = ray_kernel,
                                 parameters = true, mu_fn = mu_sinusoidal_constant,
                                 mu_fn_diff = mu_diff_sinusoidal_constant, print_level = 0)
length(events_thin)   # 438
```

Thinning gives case times but doesn't track which cases were imported, and to test the model I need the true answer. So from here on I simulate via the family tree instead: draw imported cases from $\mu(t)$, give each case $\text{Poisson}(\eta)$ children at delays of $15 + \text{Rayleigh}(\sigma)$ days, and repeat for each new generation until one comes up empty, recording every case's parent. That's `simulate_tree()` in the script (about 15 lines). The two simulators produce the same process; one run of each gave 438 and 466 cases.

```r
tree     <- simulate_tree(T_max, true)
events   <- tree$times          # used from here on
imported <- tree$parent == 0    # the truth, for marking later
```

This run had 178 imported and 288 local cases, a local fraction of 0.62 against the 0.60 expected. Here's the intensity from the package's `plot_intensity_imported()`, with each case sitting on the curve, coloured imported (red) or local (blue, which the package calls "within country"):

![Conditional intensity with imported and local cases marked](/images/hawkes/simulated_outbreak_r.png)

## Fitting the model

Now pretend we only have the case times. The parameters are fitted by maximum likelihood, and the log-likelihood of a point process is quite a nice equation:

$$
\log L = \sum_i \log \lambda^{\ast}(t_i) - \int_0^T \lambda^{\ast}(t)\,dt.
$$

The first term rewards high intensity at the moments cases actually occurred. The second, the *compensator*, charges for intensity everywhere else, so the model can't win by making $\lambda^{\ast}$ enormous. Put the rate where the cases are, and nowhere else. The integral splits into the area under the background (`mu_int_sinusoidal_constant()`) plus, for each case, $\eta$ times the fraction of its kernel that falls before the end of the window (`int_ray()`).

That's `neg_log_likelihood()`, and `ray_derivatives()` gives the analytic gradient (derived in the paper's S2 Text). Following the paper I minimise with BFGS from ten random starts and keep the best:

```r
fit_once <- function(events, start) {
  optim(par = start, fn = neg_log_likelihood, gr = ray_derivatives, method = "BFGS",
        events = events, delay = delay, kernel = ray_kernel,
        mu_fn = mu_sinusoidal_constant, mu_diff_fn = mu_diff_sinusoidal_constant,
        mu_int_fn = mu_int_sinusoidal_constant, print_level = 0)
}

fits <- lapply(1:10, function(k) fit_once(events, random_start()))
best <- fits[[which.min(sapply(fits, `[[`, "value"))]]
fit  <- as.list(best$par); fit$delay <- delay   # optim drops the fixed delay
```

where `random_start()` draws each parameter uniformly from a sensible range (`fit_epihawkes()` in the script does all of this in one go).

Most starts landed on the same answer, which is reassuring about the optimiser. The answer itself was less reassuring:

| Parameter | True | Estimate |
|---|---|---|
| $m$ | 0.10 | 0.165 |
| $a$ | 0.80 | 0.58 |
| $t_{\text{peak}}$ | 200 | 220 |
| $\eta$ | 0.60 | 0.36 |
| $\sigma$ | 20 | 22.5 |

The background is too high and $\eta$ too low, so naturally I assumed I'd done something wrong. But the analytic gradient matches a numerical one to six figures, and restarting the optimiser *at the true parameters* walks away to the same estimate, with a better likelihood than the truth. The fit is fine; the data are ambiguous. With a handful of cases a week, "more importation, less transmission" and "less importation, more transmission" produce similar case timings, and this particular five years happened to look like the former.

To see how much of this is chance, I repeated the whole thing on twenty simulated datasets (dashed lines are the true values):

![Parameter estimates across twenty simulated datasets](/images/hawkes/replicates_r.png)

The estimates scatter around the truth, with $\eta$ anywhere from 0.1 to 0.7 and the background pulled up whenever $\eta$ is pulled down. The two trade off against each other, and 450 cases isn't enough to pin both down tightly. The paper dealt with this by comparing fits across background terms, and it's worth remembering whenever you read an $R_c$ from this kind of model.

## Imported or local?

This was the most enticing part of the paper, and it's quite simple once the model is fitted. At the instant case $i$ occurs, the intensity is a sum of named contributions:

$$
\lambda^{\ast}(t_i) = \underbrace{\mu(t_i)}_{\text{importation}} + \underbrace{\phi(t_i - t_1)}_{\text{case 1}} + \underbrace{\phi(t_i - t_2)}_{\text{case 2}} + \cdots
$$

The probability that case $i$ was imported is the background's share of the total, $\mu(t_i)/\lambda^{\ast}(t_i)$, and the probability that case $j$ infected case $i$ is its share, $\phi(t_i - t_j)/\lambda^{\ast}(t_i)$. With the package's functions it's a few lines:

```r
lags    <- outer(events, events, "-")                 # lags[i, j] = t_i - t_j
contrib <- ray_kernel(lags, fit) * (lags > fit$delay) # earlier cases' pull on each case
mu_i    <- mu_sinusoidal_constant(events, fit)
p_imported <- mu_i / (mu_i + rowSums(contrib))
```

This is a probability, and should be treated as one. A case in the off-season with nothing in the preceding five weeks gets $P(\text{imported})$ near 1. A case at the peak of the season, five weeks after a cluster, is ambiguous, and the model says so. Here's how the estimates split against what actually happened:

![Estimated probability of importation, split by the truth](/images/hawkes/classification_r.png)

The two histograms overlap heavily: thresholding at 0.5 labels fewer than half the cases correctly, so you may as well guess. Case timing tells you how much transmission there is, in aggregate and over years, but very little about where any individual case came from. That's why the paper kept the travel histories: to check the population-level answer, which the model nails quite nicely.

It's worth playing with the parameters to see what moves this. A sharper kernel (say $\sigma = 5$) improves individual classification a lot, because "five weeks after a case" becomes a very precise statement. Pushing $\eta$ towards 0.9 makes chains longer, excitation dominates the intensity, and the background becomes hard to see.

## Further reading

- [Laub, Taimre & Pollett (2015)](https://arxiv.org/abs/1507.02822): a free, readable introduction to the maths, including the family-tree view. (They use $\lambda$ for the background and $\mu$ for the kernel, the reverse of this post.)
- [Rasmussen (2018)](https://arxiv.org/abs/1806.00221): short lecture notes on the conditional intensity, the likelihood and simulation.
- [Reinhart (2018)](https://doi.org/10.1214/17-STS629): a broader review, including fitting methods and applications beyond seismology.
- [Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830): the paper this post builds towards (open access), and [epihawkes](https://github.com/mrc-ide/epihawkes), the R package released with it.

If you spot any mistakes or have questions, feel free to email me. Thanks for reading!
