# timj-wood.github.io

Personal site, built with Jekyll (GitHub Pages builds it automatically on push).

## Writing a new post

1. Create `_posts/YYYY-MM-DD-slug.md`. The slug becomes the URL: `/posts/slug/`.
2. Start it with front matter:

   ```
   ---
   layout: post
   title: "Post title"
   date: YYYY-MM-DD
   summary: "One line shown on the homepage."
   ---
   ```

3. Write Markdown below. Maths uses `$$ ... $$` — on its own line for display, inline for inline.
4. Put images in `assets/img/posts/slug/` and reference them as
   `{{ '/assets/img/posts/slug/name.png' | relative_url }}`.
5. Push. The homepage Writing list updates itself.

Posts dated in the future are not published until that date. To work on something
unfinished, keep it in `_drafts/` (no date needed in the filename) and move it to
`_posts/` with a date when it's ready.

## Layout

```
_config.yml        site settings, permalink pattern
_layouts/          default.html (page chrome) and post.html (article wrapper + MathJax)
_posts/            published writing
_drafts/           unpublished writing
index.html         homepage
assets/            style.css, brand/, docs/, img/
```

## Previewing locally (optional)

```
gem install bundler jekyll
bundle init && bundle add jekyll
bundle exec jekyll serve --drafts
```
