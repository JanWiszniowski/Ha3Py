# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in params values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


project = 'Ha3Py'
copyright = '2024, Jan Wiszniowski, Andrzej Kijko'
author = 'Jan Wiszniowski, Andrzej Kijko'


def current_version():
    changelog = open('../../CHANGELOG.md', 'r')
    changelog_items = changelog.readlines()
    version = "0.0.0"
    version_date = "2025-01-01"
    for item in changelog_items:
        if item.startswith("##"):
            pos = item.split();
            version = pos[1][1:]
            version_date = pos[3]
    return version, version_date


release = '0.0.1'
release_date = '2024-12-06'
release, release_date = current_version()

# -- General params ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import sys
import os
from datetime import datetime
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.join(os.path.abspath('../..'), 'ha3py'))

# from _version import get_versions

# __version__ = get_versions()['version']
# __release_date__ = get_versions()['date']

# The full version, including alpha/_beta/rc tags.
# release = __version__

# Release date in the format "Month DD, YYYY"
# release_date = datetime.strptime(
#     __release_date__, '%Y-%m-%dT%H:%M:%S%z'
# ).strftime('%b %d, %Y')
# rst_epilog = f'\n.. |release date| replace:: {release_date}'
release_date = str(datetime.now().date())

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
    'sphinx.ext.graphviz',
    'sphinx_favicon',
    'sphinx_rtd_theme',
    'sphinxcontrib.bibtex',
    'sphinx_mdinclude',
    'matplotlib.sphinxext.plot_directive',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
bibtex_bibfiles = ['refs.bib']
bibtex_reference_style = 'author_year'
bibtex_default_style = 'unsrt'
numfig = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme_options = {
    # 'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': 'blue',
    'flyout_display': 'hidden',
    'version_selector': True,
    'language_selector': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}
html_theme = "sphinx_rtd_theme"
# html_theme = "alabaster"
# html_theme = "classic"
html_static_path = ['_static']
html_logo = '../imgs/logo_hapy.png'

# -- Options for LaTeX output -------------------------------------------------
latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',
}
