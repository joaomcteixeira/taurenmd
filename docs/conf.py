# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinxarg.ext',
    'sphinx.ext.autosectionlabel',
]
source_suffix = '.rst'
master_doc = 'index'
project = 'taurenmd'
year = '2019'
author = 'Joao MC Teixeira'
copyright = '{0}, {1}'.format(year, author)
version = release = '0.7.2'

pygments_style = 'trac'
templates_path = ['.']
extlinks = {
    'issue': ('https://github.com/joaomcteixeira/taurenmd/issues/%s', '#'),
    'pr': ('https://github.com/joaomcteixeira/taurenmd/pull/%s', 'PR #'),
    'mda': ('https://www.mdanalysis.org', 'MDAnalysis'),
    'mdtraj': ('http://mdtraj.org/1.9.3', 'MDTraj'),
    'openmm': ('https://openmmtools.readthedocs.io/en/0.18.1', 'OpenMM'),
}
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

linkcheck_ignore = [r'https://codecov.io/*']

if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = 'sphinx_rtd_theme'

html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {
    '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
}
html_short_title = '%s-%s' % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False
