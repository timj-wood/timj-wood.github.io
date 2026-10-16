---
title: Blog
permalink: /blog/
---
{% assign by_year = site.posts | group_by_exp: "post", "post.date | date: '%Y'" %}
{% for year in by_year %}
<h2 class="year">{{ year.name }}</h2>
<ol class="feed">
  {% for post in year.items %}
  <li>
    <a href="{{ post.url | relative_url }}">{{ post.title }}</a>
    <span class="when">{{ post.date | date: "%-d %B" }}{% if post.category %} · {{ post.category }}{% endif %}</span>
    {% if post.summary %}<p>{{ post.summary }}</p>{% endif %}
  </li>
  {% endfor %}
</ol>
{% endfor %}
