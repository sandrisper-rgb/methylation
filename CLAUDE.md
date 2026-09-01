# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Quarto website for an exploratory biomarker study, deployed to GitHub Pages by
`.github/workflows/publish.yml` on every push to `main`. Most of the ~60 `.qmd` files are
prose; only six execute code.

## Setup — do this first

```bash
python3 -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements-quarto.txt
python -m ipykernel install --user --name project-python \
  --display-name "Python (methylation project)"
```

The `project-python` kernel is **not optional and not just for the Python posts**.
`posts/index.qmd` is a listing page, so Quarto indexes every post on any render. Without
the kernel, *every* render fails — prose-only and R-only included — with a misleading
`ERROR: Jupyter kernel 'project-python' not found`.

On Fedora: install Quarto from the `.rpm` at quarto.org (there is no dnf repo), and
`sudo dnf install R`. The venv above is what avoids the `externally-managed-environment`
(PEP 668) error that a bare `pip install` hits on Fedora 38+. These Fedora specifics are
untested from this machine.

## Commands

```bash
quarto preview                        # live reload on port 4200
quarto render                         # full site to _site/
quarto render posts/some-post/index.qmd   # one document
```

There is no test suite (`tests/` and `R/` are empty placeholders). Correctness is
asserted with `stopifnot()` inside chunks — see `workflows/example-analysis.qmd:36` — so a
failed assertion aborts the render and a green build is the check.

## The freeze contract

`_quarto.yml` sets `freeze: auto` and `_freeze/` is **committed**. CI runs `setup-r` but
installs **no R packages**, so the R documents can only build from the freeze cache.

Quarto's freeze hash covers the whole source file, so editing a *sentence* in an R
document invalidates it and CI will try to re-execute code it cannot run. After touching
any R document — prose changes included — re-render it locally and commit the updated
`_freeze/` output alongside.

The R documents are `posts/limma-versus-dss-for-methylation.qmd` and
`workflows/example-analysis.qmd`. Rendering them locally needs `dplyr`, `ggplot2`,
`knitr`, and Bioconductor's `limma`, `bsseq`, `DSS`.
The three Python posts declare `jupyter: project-python` in their front matter.

## Authoring conventions

- Posts are either a flat `posts/name.qmd` or a `posts/name/index.qmd` directory when
  they carry figures or assets.
- Citations use `references.bib` and `apa.csl` at the repo root; posts in subdirectories
  reference them as `../../references.bib`.
- Raw HTML nested inside a `::: {.class}` fenced div must be unindented with no interior
  blank lines. Pandoc ends a raw-HTML block at a blank line, so a blank line followed by
  indented markup is silently escaped into a visible code block. `<script>` and `<style>`
  are exempt.

## Data governance

Never commit participant-level clinical or genomic data, identifiers, consent forms, or
linkage keys.
