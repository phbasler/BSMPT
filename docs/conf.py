# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BSMPT - Beyond the Standard Model Phase Transitions'
copyright = '2024, Philipp Basler, Lisa Biermann, Margarete Mühlleitner, Jonas Müller, Rui Santos and João Viana'
author = 'Philipp Basler, Lisa Biermann, Margarete Mühlleitner, Jonas Müller, Rui Santos and João Viana'
release = '3.0.7'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['breathe', 'sphinx.ext.graphviz', 'sphinx.ext.autodoc', 'sphinx.ext.viewcode',
              'sphinx.ext.napoleon', 'myst_parser', 'sphinx_github_changelog', 'sphinx_rtd_dark_mode']

# Sphinx Github Changelog Token
sphinx_github_changelog_token = os.getenv('SPHINX_GITHUB_CHANGELOG_TOKEN')

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# breathe_debug_trace_directives = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'classic'
html_theme = 'sphinx_rtd_theme'
default_dark_mode = True
html_static_path = ['_static']

# Custom CSS, word wrap
html_css_files = [
    'custom.css',
]

# Logo
html_logo = 'logo.png'

# LaTeX
latex_engine = 'xelatex'
latex_elements = {
    'fontpkg': r'''
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}
''',
    'preamble': r'''
\usepackage[titles]{tocloft}
\cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
''',
    'fncychap': r'\usepackage[Bjornstrup]{fncychap}',
    'printindex': r'\footnotesize\raggedright\printindex',
}
latex_show_urls = 'footnote'

breathe_projects = {
    "BSMPT": "/Users/joaofranciscoviana/mega/BSMPT/build/macos-armv8-release/xml/"}
breathe_default_project = "BSMPT"
