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

The tutorial pages are copies of the notebooks in `examples/`, shown with
their saved outputs.  Each notebook is built by a script,
`examples/make_*.py`: after changing it, run the script, then run the
notebook from `examples/` so that its outputs are saved.
