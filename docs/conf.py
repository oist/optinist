# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime
sys.path.insert(0, os.path.abspath('../'))
from optinist.version import VERSION, VERSION_SHORT

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#


# -- Project information -----------------------------------------------------

project = 'OptiNiSt'
copyright = f"{datetime.today().year}, OIST"
author = ''
version = VERSION_SHORT
release = VERSION

# -- readthedocs -------------------------------------------------------------
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# -- General configuration ---------------------------------------------------
master_doc = 'index'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    'sphinx.ext.napoleon',
    'sphinx.ext.coverage',
    'sphinxcontrib.apidoc',
    'sphinx_autodoc_typehints',
    'myst_parser',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.autosummary',
    'sphinx.ext.extlinks',
    'sphinx.ext.autodoc.typehints'
]

# Tell myst-parser to assign header anchors for h1-h3.
myst_heading_anchors = 4

suppress_warnings = ["myst.header"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'tests', '*test*']

source_suffix = [".rst", ".md"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_logo = '_static/optinist.png'
html_favicon = '_static/favicon.ico'

# disable document page source link
html_show_sourcelink = False

autosummary_generate = True

html_theme_options = {
    'canonical_url': '',
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'top',
    'style_external_links': False,
    'style_nav_header_background': '#C3EBE1',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}