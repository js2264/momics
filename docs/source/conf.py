#!/usr/bin/env python

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath("../../src/"))

from datetime import datetime
from importlib.metadata import metadata

info = metadata("momics")
project = info["Name"]
author = "Jacques Serizay"
copyright = f"2023-{datetime.now():%Y}, {author}"
release = info["Version"]
version = release.rsplit(".", maxsplit=1)[0]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "autoapi.extension",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx_click.ext",
    "myst_parser",
]

myst_enable_extensions = ["colon_fence", "substitution"]
myst_heading_anchors = 2
source_suffix = [".rst", ".md"]

numpydoc_show_class_members = True
napoleon_use_rtype = True
autodoc_typehints = "description"
autodoc_class_signature = "separated"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# pygments_style = "one-dark"
todo_include_todos = False
master_doc = "index"

html_static_path = ["_static"]
html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"
html_css_files = ["css/custom.css"]
html_sourcelink_suffix = ""
html_last_updated_fmt = ""
htmlhelp_basename = "momicsdoc"
html_last_updated_fmt = "%b %d, %Y"
html_title = "momics"

html_theme = "sphinx_rtd_theme"
# html_theme = "furo"
# html_theme = "pydata_sphinx_theme"
html_theme_options = {
    # "use_edit_page_button": True,
    # "external_links": "",
    "logo_only": False,
    # "github_url": "https://github.com/js2264/momics/",
    # "show_prev_next": False,
    # "search_bar_text": "Search the docs ...",
    "navigation_with_keys": False,
    "collapse_navigation": False,
    "navigation_depth": 4,
    # "show_nav_level": 2,
    # "show_toc_level": 2,
    # "footer_start": ["copyright"],
    # "footer_center": ["sphinx-version"],
    # "content_footer_items": ["last-updated"],
    # "rightsidebar": "false",
    # "relbarbgcolor": "black",
}
html_context = {
    "github_url": "https://github.com/js2264/momics/",
    "github_user": "js2264",
    "github_repo": "momics",
    "github_version": "devel",
    "doc_path": "",
}
html_static_path = ["_static"]

# -- Options for autoapi -------------------------------------------------------
autoapi_type = "python"
autoapi_dirs = ["../../src/momics"]
autoapi_ignore = ["*/cli*"]
autoapi_exclude = ["*/cli*"]
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
    "imported-members",
]
autoapi_keep_files = False
autoapi_root = "api"
autoapi_member_order = "groupwise"
