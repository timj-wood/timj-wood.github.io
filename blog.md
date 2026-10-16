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
    {% assign all_tags = "" | split: "" %}
    {% for post in site.posts %}{% assign all_tags = all_tags | concat: post.categories | concat: post.tags %}{% endfor %}
    {% assign all_tags = all_tags | uniq | sort_natural %}
    <div class="tag-chips">
      {% for tag in all_tags %}
      {% assign count = 0 %}
      {% for post in site.posts %}{% if post.categories contains tag or post.tags contains tag %}{% assign count = count | plus: 1 %}{% endif %}{% endfor %}
      <button type="button" class="tag-chip" data-tag="{{ tag | escape }}" aria-pressed="false">{{ tag }} ({{ count }})</button>
      {% endfor %}
    </div>
    <p class="tag-instruction" hidden>Click a tag above to see related posts.</p>
    {% for year in by_year %}
    <div class="tag-year">
      <h2 class="year">{{ year.name }}</h2>
      <ul class="archive-list">
        {% for post in year.items %}
        {% assign post_tags = post.categories | concat: post.tags | uniq %}
        <li data-tags="|{{ post_tags | join: '|' | escape }}|">
          <a href="{{ post.url | relative_url }}">{{ post.title }}</a>
          <time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %-d" }}</time>
        </li>
        {% endfor %}
      </ul>
    </div>
    {% endfor %}
  </div>
</div>

<script>
  (function () {
    var view = document.getElementById("by-tag");
    var chips = view.querySelectorAll(".tag-chip");
    var hint = view.querySelector(".tag-instruction");
    function select(tag) {
      hint.hidden = tag !== null;
      chips.forEach(function (c) { c.setAttribute("aria-pressed", c.dataset.tag === tag ? "true" : "false"); });
      view.querySelectorAll(".tag-year").forEach(function (group) {
        var shown = 0;
        group.querySelectorAll("li").forEach(function (li) {
          var match = tag !== null && li.dataset.tags.indexOf("|" + tag + "|") !== -1;
          li.hidden = !match;
          if (match) shown++;
        });
        group.hidden = shown === 0;
      });
    }
    chips.forEach(function (c) { c.addEventListener("click", function () { select(c.dataset.tag); }); });
    select(null);
  })();
</script>
