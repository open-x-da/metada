import os
import sys

# Add both source and build directories to Python path
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath('../../build/tests'))

project = 'METADA'
copyright = '2024'
author = 'Xin Zhang'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme' 