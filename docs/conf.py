# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Make the package importable for autodoc (repo root is one level up).
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
project = "gwfolding"
author = "Anirban Ain, Prathamesh Dalvi, Sanjit Mitra, Erik Floden"
copyright = "2026, Anirban Ain, Prathamesh Dalvi, Sanjit Mitra"
release = "0.1.0"
version = "0.1.0"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
]

# The package imports these lazily; mock them so docs build without the
# IGWN frame library (Fr) or heavy scientific deps being installed.
autodoc_mock_imports = ["Fr", "numpy", "scipy"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- HTML output -------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_logo = "../assets/logo.png"
html_theme_options = {
    "logo_only": False,
    "navigation_depth": 3,
}
