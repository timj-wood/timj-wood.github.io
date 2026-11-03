---
title: "A Hawkes process for imported and local malaria cases"
short: "Malaria Hawkes model"
summary: "A self-exciting point process fitted to five years of simulated malaria cases, to see how well case timing alone recovers the level of local transmission and the origin of each case."
when: oct 2026
links:
  - { label: repo, url: "https://github.com/timj-wood/hawkes-malaria" }
  - { label: write-up, url: /posts/hawkes-processes-part-2.html }
  - { label: the theory, url: /posts/hawkes-processes-part-1.html }
---

## Approach

[Unwin *et al.* (2021)](https://doi.org/10.1371/journal.pcbi.1008830) separated imported from locally acquired malaria cases using only their timing, by fitting a Hawkes process whose background rate stands for importation and whose kernel stands for onward transmission. To test what case timing can and cannot recover, five years of cases were simulated from a known model with a seasonal background, a branching ratio of 0.6 and a transmission delay of about 15 days, so that the true origin of every case was known. The model was then refitted to the case times alone in R, using the [epihawkes](https://github.com/mrc-ide/epihawkes) package, and the exercise was repeated on twenty replicate datasets.

## Findings

The fit was sound but the estimates were not. The branching ratio came out at 0.36 against a true value of 0.6, and across the replicates it ranged from 0.1 to 0.7, trading off against the background rate throughout. Classifying individual cases as imported or local from their fitted probabilities did no better than chance. With a few hundred cases, timing constrains the overall level of transmission far better than it identifies where any one case came from, which is why the original study kept travel histories as a check on the answer.

The code is in the [repository](https://github.com/timj-wood/hawkes-malaria). How the outbreak was simulated and fitted, with the figures, is in the [write-up]({{ '/posts/hawkes-processes-part-2.html' | relative_url }}), and the theory behind the model is in [an explainer]({{ '/posts/hawkes-processes-part-1.html' | relative_url }}).
