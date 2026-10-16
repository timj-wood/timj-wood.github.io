---
title: Blog
permalink: /blog/
hide_title: true
---
<div class="archive">
  <p class="archive-switch">
    <a href="#by-date" class="switch-date">By date</a>
    <span aria-hidden="true">|</span>
    <a href="#by-tag" class="switch-tag">By tag</a>
  </p>

  <div id="by-date" class="archive-view view-date">
    {% assign by_year = site.posts | group_by_exp: "post", "post.date | date: '%Y'" %}
    {% for year in by_year %}
    <h2 class="year">{{ year.name }}</h2>
    <ul class="archive-list">
      {% for post in year.items %}
      <li>
        <a href="{{ post.url | relative_url }}">{{ post.title }}</a>
        <time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %-d" }}</time>
      </li>
      {% endfor %}
    </ul>
    {% endfor %}
  </div>

  <div id="by-tag" class="archive-view view-tag">
    {% assign tags = site.categories | sort %}
    {% for tag in tags %}
    <h2 class="year">{{ tag[0] }}</h2>
    <ul class="archive-list">
      {% for post in tag[1] %}
      <li>
        <a href="{{ post.url | relative_url }}">{{ post.title }}</a>
        <time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %-d, %Y" }}</time>
      </li>
      {% endfor %}
    </ul>
    {% endfor %}
  </div>
</div>
