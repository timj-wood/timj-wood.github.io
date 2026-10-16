---
title: Projects
permalink: /projects/
hide_title: true
---
<p class="page-intro">Side projects, mostly things I wanted to understand properly by building them myself. For academic work, head over to <a href="{{ '/research/' | relative_url }}">research</a>.</p>

{% if site.projects.size > 0 %}
{% include entries.html items=site.projects inline=true %}
{% else %}
<p class="entry-empty">Nothing here yet.</p>
{% endif %}
