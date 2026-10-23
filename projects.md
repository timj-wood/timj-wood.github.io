---
title: Projects
permalink: /projects.html
redirect_from: /projects/
hide_title: true
---
{%- comment -%}
  One <section class="entry"> per project: title, then the figure at full width (linking to the repo),
  then the date, a one- or two-sentence summary and links. The full write-up lives in the blog post.
  Copy a block to add a new one.
{%- endcomment -%}
<p class="page-intro">Side projects, mostly things I wanted to understand properly by building them myself. For academic work, head over to <a href="{{ '/research.html' | relative_url }}">research</a>.</p>

<section class="entry" id="splice-site-nn">
  <h2>Can a mini neural network learn where a gene is spliced?</h2>
  <div class="entry-body stacked">
    <a class="entry-figure-link" href="https://github.com/timj-wood/splice-site-nn" target="_blank" rel="noopener"><img class="entry-figure" src="{{ '/images/splice/neuron_weights.png' | relative_url }}" alt="Weights of a single neuron trained to detect splice junctions, laid out by sequence position and base, with the strongest weights on AG before the junction and GT after it"></a>
    <p class="entry-meta">Oct 2026</p>
    <p class="entry-summary">My first (mini) neural network, built from a single neuron with hand-written gradient descent up to a small multi-layer perceptron in PyTorch, on 3,190 DNA sequences. The model rediscovers the GT…AG intron consensus by itself and, usefully, does not outperform logistic regression.</p>
    <p class="entry-links"><a href="https://github.com/timj-wood/splice-site-nn">Repo</a>, Blog (Coming soon)</p>
  </div>
</section>
