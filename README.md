# Woodworks

Jekyll on GitHub Pages.

## Structure

```
index.md              home page: bio + projects + writing feed  (layout: home)
about.md              /about/
_posts/               YYYY-MM-DD-slug.md → /posts/slug/
_drafts/              unpublished posts (no date needed)
_projects/            one file per project → /projects/name/
_layouts/             default, home, page, post, project
_includes/            head (fonts, SEO tags, MathJax), footer
assets/css/style.css  the one stylesheet
assets/docs/          cv.pdf, report PDFs
assets/img/           images
```

## Front matter

Posts: `title`, `summary` (shown on home + as standfirst), `category` (free text), `math: true` to load MathJax.
Projects: `title`, `summary`, `period`, `order` (sort), `with`, `repo`, `pdf`.

## Local preview

```
bundle install
bundle exec jekyll serve --drafts --livereload
```

## Deploy

Push to `main` of `timj-wood/timj-wood.github.io`. GitHub Pages builds it. Custom domain: add a `CNAME` file containing the domain and point DNS at GitHub.

## Migrating from the old site

1. Copy the body of each old post into `_posts/YYYY-MM-DD-<old-slug>.md` — URLs are preserved as `/posts/<slug>/`.
2. Move `cv.pdf` and the report PDFs into `assets/docs/`.
3. Delete the old `index.html`, `posts/` and `assets/scan.js`.
