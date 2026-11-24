"""Sphinx configuration for the ArgoPAL documentation (Read the Docs)."""

import os
import shutil
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import argopal


# The tutorials live in examples/; copies are rendered with their saved
# outputs (they download a float from the GDAC, so they are not re-run).
TUTORIALS = [
    "argopal_tutorial.ipynb",
    "pal_data_tutorial.ipynb",
    "multi_float_tutorial.ipynb",
    "depth_correction_tutorial.ipynb",
]
os.makedirs(os.path.join(HERE, "tutorials"), exist_ok=True)
for name in TUTORIALS:
    shutil.copyfile(
        os.path.join(ROOT, "examples", name),
        os.path.join(HERE, "tutorials", name),
    )

project = "ArgoPAL"
author = "Louise Delaigue"
copyright = "2026, Louise Delaigue"
release = argopal.__version__
version = release

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
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
autodoc_default_options = {"members": False}
autodoc_typehints = "none"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# Theme
html_theme = "furo"
html_title = f"ArgoPAL {release}"
html_logo = "img/argopal_logo.png"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "source_repository": "https://github.com/louisedelaigue/ArgoPAL",
    "source_branch": "main",
    "source_directory": "docs/",
}
copybutton_prompt_text = r">>> |\$ "
copybutton_prompt_is_regexp = True
