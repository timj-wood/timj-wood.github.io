---
title: Scottish wildcat demographic models
summary: Fitting 2D dadi models to a folded joint site frequency spectrum for wildcat and domestic cat.
period: Jun - Sept 2026
order: 1
with: Mark Beaumont, Dennis Prangle and Grace Yan
repo: https://github.com/timj-wood/wildcat_models
pdf: /assets/docs/dadi_report.pdf
math: true
---
I fitted three two-population demographic models to a folded joint site frequency spectrum built from 16 wild-caught Scottish wildcats and 6 domestic cats, using 22.5 Mb of sequence on chromosomes A1 and A2 and dadi 2.4.4. The models share one topology — an ancestral population splits, then the branches exchange migrants — and differ in how much structure they impose on population size and gene flow.

Because dadi maximises a composite likelihood, I compared models with the composite-likelihood AIC and built confidence intervals from the Godambe information matrix, estimated from 100 block-bootstrap replicates. The effective parameter count was 293 against a nominal 12, which is how far the uncorrected likelihood overstates its own certainty. The numerical details are written up separately.

A secondary-contact model fitted 806 CLAIC units worse than a model with size changes at domestication and the loss of the British land bridge. The best-scoring model finished on a parameter bound, so the simpler one is reported: split at roughly 485,000 years ago, domestic size change at roughly 10,600 years ago. Only four parameters are known to within a factor of two, and the correlation matrix shows the deeper parameters lie along a single ridge in the likelihood — the interval widths are not independent. The wildcat size change lands at 673 years ago, matching no known event, so that label should be dropped.

The analysis is intended as a baseline for the simulation-based inference methods. 