# Woodworks

Jekyll on GitHub Pages.

## Structure

```
index.md              landing page = About: photo, name, role, contact address, bio  (layout: home)
about.md              redirects old /about/ links to /
blog.md               /blog/ — all posts grouped by year
research.md           /research/ — research entries, shown in full on one page
projects.md           /projects/ — other projects, same format
_data/navigation.yml  header links and the links at the end of every post
_posts/               YYYY-MM-DD-slug.md → /posts/slug/
_drafts/              unpublished posts (no date needed)
_research/            one file per research entry (no page of its own)
_projects/            one file per project (no page of its own)
_layouts/             default, home, page, post
_includes/            head (fonts, SEO tags, MathJax), footer, entries (research/projects listing), end-matter ("written on" line + related posts under each post)
assets/css/style.css  the one stylesheet
assets/docs/          cv.pdf, report PDFs
assets/img/           images
```

## Front matter

Posts: `title`, `summary` (shown on home + as standfirst), `category` (free text), `math: true` to load MathJax.
About (index.md): `portrait`, `role`.
Posts in the same `category` are listed under each other as "More on …".
Research / projects: `title`, `period`, `order` (sort), `with`, and optional `repo`, `report` (PDF path), `blog` (post URL), `figure` + `figure_alt` (figure shown on the left). The body text is shown in full.

## Local preview

```
bundle install
bundle exec jekyll serve --drafts --livereload
```

## Deploy

Push to `main` of `timj-wood/timj-wood.github.io`. GitHub Pages builds it. Custom domain: add a `CNAME` file containing the domain and point DNS at GitHub.
