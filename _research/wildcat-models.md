---
title: Scottish wildcat demographic models
summary: Fitting 2D dadi models to a folded joint site frequency spectrum for wildcat and domestic cat.
period: Jun - Sept 2026
order: 1
with: Mark Beaumont, Dennis Prangle and Grace Yan
repo: https://github.com/timj-wood/wildcat_models
report: /assets/docs/dadi_report.pdf
blog: /posts/composite-likelihood/
figure: /assets/img/corr.png
figure_alt: Parameter correlation matrix from the Godambe-corrected covariance
---
Summer internship in Bristol's School of Mathematics. I fitted three two-population demographic models in dadi to the joint site frequency spectrum of 16 Scottish wildcats and 6 domestic cats, compared them with the composite-likelihood AIC, and corrected their confidence intervals using the Godambe information matrix. The uncorrected intervals overstated the certainty considerably, and only four parameters could be estimated to within a factor of two.

The preferred model places the wildcat–domestic cat split at roughly 485,000 years ago, with a domestic size change at roughly 10,600 years ago. It is intended as a baseline for later simulation-based inference work.
