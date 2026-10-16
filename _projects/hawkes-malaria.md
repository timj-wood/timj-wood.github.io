---
title: Hawkes processes for imported and local malaria
summary: Simulating and fitting a Hawkes process to separate imported and locally acquired malaria cases.
period: Oct 2026
order: 1
repo: https://github.com/timj-wood/hawkes-malaria
blog: /posts/hawkes-processes/
figure: /assets/posts/hawkes/simulated_outbreak_r.png
figure_alt: Simulated malaria outbreak with imported and locally acquired cases
---
An R walkthrough of the method from Unwin *et al.* (2021), which uses Hawkes processes to separate imported from locally acquired malaria cases using only their timing. Using the [epihawkes](https://github.com/mrc-ide/epihawkes) package, I simulated five years of a seasonal outbreak, fitted the model, and repeated the fit on twenty replicate datasets. The model recovers overall transmission reasonably well, but case timing alone says little about whether any individual case was imported.
