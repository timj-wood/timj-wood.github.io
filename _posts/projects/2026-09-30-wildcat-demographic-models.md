---
title: "Scottish wildcat demographic models"
short: "Wildcat demographic models"
summary: "Composite-likelihood demographic inference for the Scottish wildcat and the domestic cat from a folded joint site frequency spectrum."
affiliation: University of Bristol
when: jun – sept 2026
with: Mark Beaumont, Dennis Prangle and Grace Yan
links:
  - { label: repo, url: "https://github.com/timj-wood/wildcat_models" }
  - { label: blog, url: /posts/composite-likelihood.html }
---

## Approach

A summer internship with the School of Mathematics. Three two-population demographic models of the split between Scottish wildcats and domestic cats were fitted with [dadi](https://dadi.readthedocs.io/en/latest/) to a folded joint site frequency spectrum built from 16 wild-caught wildcats and 6 domestic cats, using 22.5 Mb of callable sequence on chromosomes A1 and A2. Because dadi maximises a composite likelihood, its model-selection criteria and confidence intervals are over-optimistic, so models were compared with the composite-likelihood AIC and intervals were taken from the Godambe information matrix, both estimated from 100 block-bootstrap replicates.

## Findings

A model in which the two populations were isolated and then came back into contact fitted substantially worse than two models placing changes in population size and migration at domestication and at the loss of the British land bridge. Of those two, the model with continuous size change ended with a parameter on its bound, so the model with instantaneous size changes is the primary result. It places the split between the *silvestris* and *lybica* lineages at roughly 485,000 years ago and the size change in the domestic lineage at roughly 10,600 years ago. Four parameters are estimated with narrow intervals; the others only to within a factor of four or more, and the correlation matrix shows the deeper parameters lying along a single ridge in the likelihood, so the intervals cannot be read independently. One estimate, a size change in the wildcat lineage 673 years ago, matches no known event.

The code for this project is in the accompanying [repository](https://github.com/timj-wood/wildcat_models). The composite-likelihood corrections, and the mistakes made applying them, are the subject of a [blog]({{ '/posts/composite-likelihood.html' | relative_url }}) post.
