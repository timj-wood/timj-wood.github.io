# Woodworks

Tim Wood's site. Jekyll on GitHub Pages, with the folder layout of [stmorse.github.io](https://github.com/stmorse/stmorse.github.io).

## Structure

```
index.html         About page (landing page): photo, name, role, contact, bio
blog.md            /blog/ — all posts grouped by year
research.md        /research/ — research entries, written straight into the page
projects.md        /projects/ — side projects, same format
404.html           not-found page
_config.yml        title, author links, plugins
_layouts/          default, page, post
_includes/         head, header, menu (header links), footer, social-icons, related-posts
_posts/            YYYY-MM-DD-slug.md → /posts/slug/
_drafts/           unpublished posts
css/main.css       the one stylesheet
images/            images; post figures in images/<post>/
docs/              CV, report PDFs, supporting scripts
```

## Front matter

Posts: `title`, `summary`, `category` (posts in the same category are listed first under "Related posts"), `math: true` to load MathJax.

## Adding research or a project

Copy one of the `<section class="entry">` blocks in `research.md` or `projects.md`.

## Local preview

```
bundle install
bundle exec jekyll serve --drafts --livereload
```
