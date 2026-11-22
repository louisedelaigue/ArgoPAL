# Development

## Install for development

    git clone https://github.com/louisedelaigue/ArgoPAL.git
    cd ArgoPAL
    pip install -e ".[ais,docs]"
    pip install pytest

## Tests

```{include} ../README.md
---
start-after: "## Tests"
relative-docs: "docs/"
relative-images: true
---
```

## Documentation

This site is built with Sphinx from `docs/`, on Read the Docs at every
push.  To build it locally:

    pip install -e ".[docs]"
    sphinx-build -b html docs docs/_build/html

The tutorial page is a copy of `examples/argopal_tutorial.ipynb`, shown
with its saved outputs.  After changing the tutorial, rebuild and re-run
it (`python examples/make_tutorial.py`, then run it from `examples/`).
