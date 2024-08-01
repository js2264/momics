#!/usr/bin/env python

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath("../../src/momics/"))

from datetime import datetime
from importlib.metadata import metadata

info = metadata("momics")
project = info["Name"]
author = "Jacques Serizay"
copyright = f"2023-{datetime.now():%Y}, {author}."
copyright = f"{datetime.now():%Y}, {author}."
release = info["Version"]
version = release.rsplit(".", maxsplit=1)[0]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "recommonmark",
]

autosummary_generate = True
numpydoc_show_class_members = True
napoleon_use_rtype = True
autodoc_typehints = "description"
autodoc_class_signature = "separated"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
pygments_style = "sphinx"
todo_include_todos = False
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "furo"
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
htmlhelp_basename = "momicsdoc"
html_last_updated_fmt = "%b %d, %Y"
