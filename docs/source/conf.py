import os
import sys

sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------

project = 'phenocor'
copyright = '2026, Jing Chen'
author = 'Jing Chen'
release = '0.2.0'

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",    
    "sphinx.ext.napoleon",    
    "nbsphinx",           
]

templates_path = ['_templates']
exclude_patterns = []

# -- HTML output -------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

nbsphinx_show_input_prompt = False
nbsphinx_execute = "auto"

autosummary_imported_members = False
add_module_names = False

autosummary_generate = False
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "private-members": False
}