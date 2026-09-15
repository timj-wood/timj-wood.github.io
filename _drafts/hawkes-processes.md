---
layout: post
title: "Events that cause events: Hawkes processes and malaria elimination"
date: 2026-09-20
summary: "A literature review of self-exciting point processes, centred on their use in modelling malaria transmission in near-elimination settings."
---

*[~100 words. A case appears in a district that is almost malaria-free. Was it imported, or seeded by a case last week? The question is about timing — that is what this post is about.]*

## Why a constant rate isn't enough

*[Poisson process in one paragraph: events arrive at a steady rate with no memory.]*

$$
\lambda(t) = \lambda \quad \text{(constant)}
$$

*[Where it fails: aftershocks, outbreaks, neurons — events cluster. State the idea in one sentence: an event temporarily raises the chance of another.]*

## The Hawkes process

*[Intensity in words first: baseline rate plus a bump from every past event, which fades.]*

$$
\lambda(t \mid \mathcal{H}_t) = \mu + \sum_{t_i < t} \phi(t - t_i)
$$

where $$\mu$$ is the baseline rate, $$\mathcal{H}_t$$ is the history of events before $$t$$, and $$\phi$$ is the *kernel* — the shape of the memory.

*[The exponential kernel: the simplest and most common choice.]*

$$
\phi(s) = \alpha\, e^{-\beta s}, \qquad s \ge 0
$$

$$\alpha$$ is the jump in the rate immediately after an event; $$\beta$$ is how fast that jump decays. Put together:

$$
\lambda(t) = \mu + \alpha \sum_{t_i < t} e^{-\beta (t - t_i)}
$$

<figure>
  <img src="{{ '/assets/img/posts/hawkes-processes/fig1-intensity.png' | relative_url }}" alt="Intensity of a Hawkes process over time, with event times marked as ticks beneath the curve.">
  <figcaption>Figure 1. [Intensity for a simulated exponential-kernel Hawkes process, event times as ticks. Compare against a Poisson process with the same mean rate.]</figcaption>
</figure>

### The branching ratio

*[Each event's total contribution is the area under its kernel — the expected number of "offspring" events it triggers.]*

$$
n = \int_0^\infty \phi(s)\, \mathrm{d}s = \frac{\alpha}{\beta}
$$

*[Stability: if $$n < 1$$ bursts die out; if $$n \ge 1$$ the process explodes. Explicit $$R_0$$ analogy here.]*

$$
\text{stable} \iff n = \frac{\alpha}{\beta} < 1
$$

*[Long-run average rate, when stable — each baseline event heads a whole family tree:]*

$$
\bar{\lambda} = \frac{\mu}{1 - n}
$$

*[One line on other kernels — power-law for earthquakes, gamma or shifted kernels for disease serial intervals.]*

## A short history

*[Hawkes 1971 — what problem he was solving.]*

*[Ogata 1981 — simulating one by thinning. Ogata 1988 — ETAS and why seismology adopted it first.]*

*[Spread to finance, social media, epidemiology.]*

## Fitting a Hawkes process to data

*[Maximum likelihood, conceptually: pick $$\mu, \alpha, \beta$$ that make the observed timings most probable.]*

$$
\ell(\mu, \alpha, \beta) = \sum_{i=1}^{N} \log \lambda(t_i) \;-\; \int_0^T \lambda(t)\, \mathrm{d}t
$$

The first term rewards high intensity *at* the events; the second penalises high intensity everywhere else.

*[Why the exponential kernel makes this cheap: the sum inside $$\lambda(t_i)$$ can be updated recursively rather than recomputed — one sentence, no derivation.]*

$$
A_i = e^{-\beta (t_i - t_{i-1})}\,(1 + A_{i-1}), \qquad A_1 = 0, \qquad \lambda(t_i) = \mu + \alpha A_i
$$

*[Goodness of fit: the time-rescaling idea in words. If the model is right, transforming each event time by the integrated intensity should give a unit-rate Poisson process.]*

$$
\tau_i = \int_0^{t_i} \lambda(t)\, \mathrm{d}t \qquad \Rightarrow \qquad \tau_i - \tau_{i-1} \sim \text{Exp}(1)
$$

*[Software pointer: `tick`, `HawkesPyLib`.]*

## Case study: malaria in near-elimination settings

*[Setting: Unwin et al. 2021 — China and Eswatini. Why near-elimination is the interesting regime: few cases, and the question is whether they are imported or local.]*

| Model term | Meaning for malaria |
|---|---|
| $$\mu$$, baseline | *[Imported cases arriving from outside the region]* |
| Self-excitation, $$\sum \phi(t - t_i)$$ | *[Local onward transmission via mosquitoes]* |
| Kernel $$\phi$$ | *[Serial interval — time from one case to the cases it causes, including the mosquito stage]* |
| Branching ratio $$n$$ | *[Effective reproduction number $$R_{\text{eff}}$$]* |

*[What "semi-mechanistic" means here — where the biology enters the kernel.]*

$$
\phi(s) = \; ? \qquad \text{[copy the kernel form from the paper]}
$$

*[Main findings, in your own words.]*

*[What this can do that a compartmental (SIR-type) model can't — and vice versa.]*

*[Anything published since 2021 — check.]*

## The same idea elsewhere

*[COVID-19 — discrete-time self-exciting model on daily counts. One paragraph.]*

$$
\lambda_t = \mu_t + \sum_{s < t} \phi(t - s)\, y_s
$$

*[Earthquakes — ETAS. One paragraph. Note the power-law kernel and magnitude-dependent productivity.]*

$$
\lambda(t) = \mu + \sum_{t_i < t} K\, e^{\alpha (M_i - M_0)}\, (t - t_i + c)^{-p}
$$

*[Neurons — two sentences.]*

## Limitations and open directions

*[Kernel choice; non-stationary baselines $$\mu(t)$$; multivariate versions.]*

$$
\lambda_i(t) = \mu_i + \sum_{j=1}^{D} \sum_{t_{j,k} < t} \phi_{ij}(t - t_{j,k})
$$

*[Likelihood becomes intractable for richer versions — simulation-based inference as the way round it (SB-ETAS). One sentence connecting to your own MSc direction.]*

## Further reading

1. Laub, Taimre & Pollett (2015). *Hawkes Processes.* [arXiv:1507.02822](https://arxiv.org/abs/1507.02822). *[one-line note]*
2. Rizoiu, Lee, Mishra & Xie (2017). *A Tutorial on Hawkes Processes for Events in Social Media.* [arXiv:1708.06401](https://arxiv.org/abs/1708.06401). *[note]*
3. Hawkes (1971). *Spectra of some self-exciting and mutually exciting point processes.* Biometrika 58(1):83–90.
4. Ogata (1981). *On Lewis' simulation method for point processes.* IEEE Trans. Inf. Theory 27(1):23–31.
5. Ogata (1988). *Statistical models for earthquake occurrences and residual analysis for point processes.* JASA 83(401):9–27.
6. Unwin et al. (2021). *Using Hawkes Processes to model imported and local malaria cases in near-elimination settings.* [PLOS Comput Biol 17(4):e1008830](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008830). *[note]*
7. *[COVID-19 discrete-time self-exciting paper, PLOS ONE 2021]*
8. Stockman, Lawson & Werner (2024). *SB-ETAS.* [arXiv:2404.16590](https://arxiv.org/abs/2404.16590). *[note]*
