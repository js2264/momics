#!/usr/bin/env python

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

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

numpydoc_show_class_members = False
napoleon_use_rtype = False
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

html_theme = "furo"
html_static_path = ["_static"]
htmlhelp_basename = "momicsdoc"

# -- Options for manual page output ---------------------------------------
man_pages = [(master_doc, "momics", "momics Documentation", [author], 1)]

# -- Options for Texinfo output -------------------------------------------
texinfo_documents = [
    (
        master_doc,
        "momics",
        "momics Documentation",
        author,
        "momics",
        "One line description of project.",
        "Miscellaneous",
    ),
]
