---
title: Projects
permalink: /projects/
hide_title: true
---
{%- comment -%}
  One <section class="entry"> per project: title, then the figure at full width (linking to the repo),
  then the date, a one- or two-sentence summary and links. The full write-up lives in the blog post.
  Copy a block to add a new one.
{%- endcomment -%}
<p class="page-intro">Side projects, mostly things I wanted to understand properly by building them myself. For academic work, head over to <a href="{{ '/research/' | relative_url }}">research</a>.</p>

<section class="entry" id="hawkes-malaria">
  <h2>Hawkes processes for imported and local malaria</h2>
  <div class="entry-body stacked">
    <a class="entry-figure-link" href="https://github.com/timj-wood/hawkes-malaria" target="_blank" rel="noopener"><img class="entry-figure" src="{{ '/images/hawkes/simulated_outbreak_r.png' | relative_url }}" alt="Simulated malaria outbreak with imported and locally acquired cases"></a>
    <p class="entry-meta">Oct 2026</p>
    <p class="entry-summary">Simulating and fitting a Hawkes process in R to separate imported from locally acquired malaria cases, following Unwin <em>et al.</em> (2021).</p>
    <p class="entry-links"><a href="https://github.com/timj-wood/hawkes-malaria">Repo</a>, <a href="{{ '/posts/hawkes-processes/' | relative_url }}">Blog</a></p>
  </div>
</section>
