---
title: Projects
permalink: /projects/
hide_title: true
---
{%- comment -%}
  One <section class="entry"> per project: figure on the left (linking to the repo), text on the right
  with the title in bold at the start. Copy a block to add a new one.
{%- endcomment -%}
<p class="page-intro">Side projects, mostly things I wanted to understand properly by building them myself. For academic work, head over to <a href="{{ '/research/' | relative_url }}">research</a>.</p>

<section class="entry" id="hawkes-malaria">
  <div class="entry-body has-figure">
    <a class="entry-figure-link" href="https://github.com/timj-wood/hawkes-malaria" target="_blank" rel="noopener"><img class="entry-figure" src="{{ '/images/hawkes/simulated_outbreak_r.png' | relative_url }}" alt="Simulated malaria outbreak with imported and locally acquired cases"></a>
    <div class="entry-text">
      <p><strong class="entry-title">Hawkes processes for imported and local malaria.</strong> An R walkthrough of the method from Unwin <em>et al.</em> (2021), which uses Hawkes processes to separate imported from locally acquired malaria cases using only their timing. Using the <a href="https://github.com/mrc-ide/epihawkes">epihawkes</a> package, I simulated five years of a seasonal outbreak, fitted the model, and repeated the fit on twenty replicate datasets. The model recovers overall transmission reasonably well, but case timing alone says little about whether any individual case was imported.</p>
      <p class="entry-links"><a href="https://github.com/timj-wood/hawkes-malaria">Repo</a>, <a href="{{ '/posts/hawkes-processes/' | relative_url }}">Blog</a></p>
    </div>
  </div>
</section>
