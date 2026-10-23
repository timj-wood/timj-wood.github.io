---
title: "Hawkes processes in epidemiology (Part II: code)"
redirect_from: /posts/hawkes-processes-part-2/
category: stats
tags: hawkes
---
In [part 1]({% post_url 2026-10-13-hawkes-processes-part-1 %}) we met the Hawkes process and the idea behind [Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830): fit a model to the timing of malaria cases, then read off how much transmission there is and which cases were probably imported. This post puts it to the test. The best way I found to understand the method was to simulate an outbreak where the true transmission tree is known, fit the model to the case times alone, and see what it recovers.

If you haven't read part 1, the short version is this. The rate of new cases at time $t$, given everything that has happened so far, is

$$
\lambda^{\ast}(t) = \mu(t) + \sum_{t_i < t} \eta\, g(t - t_i),
$$

where $\mu(t)$ is the **background rate** of imported cases, $\eta$ is the **branching ratio** (the average number of local cases each case causes, a reproduction number), and $g$ is the distribution of the **delay** between a case and the cases it causes. The probability that a case was imported is the background's share of the intensity at the moment it occurred.

I used the simulator, likelihood, gradient and intensity plot from Unwin, Mishra and Bhatt's excellent [epihawkes](https://github.com/mrc-ide/epihawkes) package, and my scripts are at [hawkes-malaria](https://github.com/timj-wood/hawkes-malaria). The package isn't on CRAN, so install it from GitHub:

```r
remotes::install_github("mrc-ide/epihawkes")
library(epihawkes)
```

## Setting up the outbreak

The setup (numbers made up for illustration) is five years of cases where:

- Imported cases arrive at a seasonal background rate, $\mu(t) = m\,\big(1 + a\cos(2\pi(t - t_{\text{peak}})/365.25)\big)$, averaging $m = 0.1$ per day ($\approx 3$ per month), with a strong annual cycle ($a = 0.8$) peaking around day 200.
- Each case causes, on average, $\eta = 0.6$ local cases.
- The delay $g$ is zero for 15 days, then follows a Rayleigh curve with spread $\sigma = 20$ days, so onward cases are most likely about five weeks later.

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

## Simulating

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

The estimates scatter around the truth, with $\eta$ anywhere from 0.1 to 0.7 and the background pulled up whenever $\eta$ is pulled down. The two trade off against each other, and 450 cases isn't enough to pin both down tightly. The paper dealt with this by comparing fits across different background terms, and it's worth remembering whenever you read an $R_c$ from this kind of model.

## Imported or local?

From part 1, the probability that case $i$ was imported is the background's share of the intensity when it happened, $\mu(t_i)/\lambda^{\ast}(t_i)$, and the probability that case $j$ infected it is the share of $j$'s kernel. With the package's functions it's a few lines:

```r
lags    <- outer(events, events, "-")                 # lags[i, j] = t_i - t_j
contrib <- ray_kernel(lags, fit) * (lags > fit$delay) # earlier cases' pull on each case
mu_i    <- mu_sinusoidal_constant(events, fit)
p_imported <- mu_i / (mu_i + rowSums(contrib))
```

This is a probability, and should be treated as one. Here's how the estimates split against what actually happened:

![Estimated probability of importation, split by the truth](/images/hawkes/classification_r.png)

The two histograms overlap heavily: thresholding at 0.5 labels fewer than half the cases correctly, so you may as well guess. Case timing tells you how much transmission there is, in aggregate and over years, but very little about where any individual case came from. That's why the paper kept the travel histories: to check the population-level answer, which the model nails quite nicely.

It's worth playing with the parameters to see what moves this. A sharper kernel (say $\sigma = 5$) improves individual classification a lot, because "five weeks after a case" becomes a very precise statement. Pushing $\eta$ towards 0.9 makes chains longer, excitation dominates the intensity, and the background becomes hard to see.

## Wrapping up

Two things stuck with me. First, the maths of the Hawkes process is clean enough that a biologist can follow it, and the fitting is a few dozen lines of R on top of `epihawkes`. Second, the model's honest answer is a population-level one: with this many cases, $\eta$ and the background are hard to separate, and individual cases are mostly unclassifiable from timing alone. The clever part of Unwin *et al.* isn't that the model labels cases, it's that it estimates how much local transmission is left, and that the travel histories agree.

The scripts for everything above are at [hawkes-malaria](https://github.com/timj-wood/hawkes-malaria). If you spot any mistakes or have questions, feel free to email me. Thanks for reading!
