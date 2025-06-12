# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import importlib.metadata
import sys
from os.path import abspath, dirname, join

sys.path.insert(0, abspath("."))

# Fetch the `__version__`
proj_folder = dirname(dirname(__file__))
bsb_init_file = join(proj_folder, "pyproject.toml")

# -- Project information -----------------------------------------------------

project = "DBBS Cerebellar models"
copyright = "2025, DBBS University of Pavia"
author = "DBBS University of Pavia"
release = importlib.metadata.version("cerebellar_models")


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx_extensions",
]
autodoc_type_aliases = {"JobPool": "bsb.services.pool.JobPool"}
autodoc_mock_imports = ["nest"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://scipy.github.io/devdocs/", None),
    "bsb": ("https://bsb.readthedocs.io/en/latest", None),
    "bsb-core": ("https://bsb-core.readthedocs.io/en/latest", None),
    "neo": ("https://neo.readthedocs.io/en/latest/", None),
}

autoclass_content = "both"
autodoc_typehints = "both"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

numfig = True  # allow figure numbering

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

html_context = {
    "maintainer": "Dimitri Rodarie",
    "project_pretty_name": "DBBS Cerebellar Models",
    "projects": {"DBBS Scaffold": "https://github.com/dbbs/cerebellar-models"},
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
