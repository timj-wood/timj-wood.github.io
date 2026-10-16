---
title: Research
permalink: /research/
redirect_from:
  - /projects/wildcats/
  - /projects/wildcat-introgression/
---
{%- comment -%}
  One <section class="entry"> per piece of work: figure on the left, text on the right.
  Copy a block to add a new one; drop the <img> and the "has-figure" class if there’s no figure.
{%- endcomment -%}

<section class="entry" id="wildcat-models">
  <h2>Scottish wildcat demographic models</h2>
  <div class="entry-body has-figure">
    <img class="entry-figure" src="{{ '/images/corr.png' | relative_url }}" alt="Parameter correlation matrix from the Godambe-corrected covariance">
    <div class="entry-text">
      <p class="entry-meta">Jun - Sept 2026 · with Mark Beaumont, Dennis Prangle and Grace Yan</p>
      <p>Summer internship in Bristol’s School of Mathematics. I fitted three two-population demographic models in dadi to the joint site frequency spectrum of 16 Scottish wildcats and 6 domestic cats, compared them with the composite-likelihood AIC, and corrected their confidence intervals using the Godambe information matrix. The uncorrected intervals overstated the certainty considerably, and only four parameters could be estimated to within a factor of two.</p>
      <p>The preferred model places the wildcat–domestic cat split at roughly 485,000 years ago, with a domestic size change at roughly 10,600 years ago. It is intended as a baseline for later simulation-based inference work.</p>
      <p class="entry-links"><a href="https://github.com/timj-wood/wildcat_models">Repo</a>, <a href="{{ '/docs/dadi_report.pdf' | relative_url }}">Report</a>, <a href="{{ '/posts/composite-likelihood/' | relative_url }}">Blog</a></p>
    </div>
  </div>
</section>

<section class="entry" id="wildcat-introgression">
  <h2>Introgression and selection in Scottish wildcats</h2>
  <div class="entry-body has-figure">
    <img class="entry-figure" src="{{ '/images/manhattan.png' | relative_url }}" alt="Manhattan plot of pcadapt outlier p-values across the cat genome, with candidate genes labelled">
    <div class="entry-text">
      <p class="entry-meta">Sept - Dec 2025 · with Mark Beaumont</p>
      <p>BSc dissertation. I used pcadapt to scan genome-wide SNP data from 108 cats for regions shaping population structure in the European wildcat, a species threatened by hybridisation with domestic cats. Seven candidate genes with strong outlier signals pointed to GABA transport, taurine transport and calcium signalling, suggesting that behavioural and metabolic pathways may be under selection. PCA showed patterns of introgression, local adaptation and regional divergence, which could inform conservation efforts but would need further testing.</p>
      <p class="entry-links"><a href="{{ '/docs/bsc_diss.pdf' | relative_url }}">Report</a></p>
    </div>
  </div>
</section>
