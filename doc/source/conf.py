# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import sys, os
sys.path.insert(0, os.path.abspath('../../fast_response/'))
#import fast_response

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FastResponseAnalysis'
copyright = '2023, Alex Pizzuto, Jessie Thwaites'
author = 'Alex Pizzuto, Jessie Thwaites'
release = '1.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autodoc',
]

templates_path = ['_templates']
exclude_patterns = []

autodoc_mock_imports = ['skylab', 'icecube']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'classic'
html_theme_options = {
    "sidebarwidth": 310,
}
html_static_path = ['_static']
