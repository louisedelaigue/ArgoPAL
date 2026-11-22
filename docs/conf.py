"""Sphinx configuration for the ArgoPAL documentation (Read the Docs)."""

import os
import shutil
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import argopal


# The tutorial lives in examples/; a copy is rendered with its saved
# outputs (it downloads a float from the GDAC, so it is not re-run here).
shutil.copyfile(
    os.path.join(ROOT, "examples", "argopal_tutorial.ipynb"),
    os.path.join(HERE, "tutorial.ipynb"),
)

project = "ArgoPAL"
author = "Louise Delaigue"
copyright = "2026, Louise Delaigue"
release = argopal.__version__
version = release

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Markdown (MyST) and notebooks
myst_enable_extensions = ["colon_fence", "dollarmath"]
myst_heading_anchors = 3
nb_execution_mode = "off"
# README sections are included into pages under their own heading
suppress_warnings = ["myst.header"]

# API reference from the NumPy-style docstrings; `name` renders as code
default_role = "literal"
autosummary_generate = True  # one page per function, in docs/api/generated
autodoc_default_options = {"members": False}
autodoc_typehints = "none"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# Theme
html_theme = "furo"
html_title = f"ArgoPAL {release}"
html_logo = "img/argopal_logo.png"
html_static_path = []
html_theme_options = {
    "source_repository": "https://github.com/louisedelaigue/ArgoPAL",
    "source_branch": "main",
    "source_directory": "docs/",
}
copybutton_prompt_text = r">>> |\$ "
copybutton_prompt_is_regexp = True
