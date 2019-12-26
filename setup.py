#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import, print_function

import io
import re
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages, setup


def read(*names, **kwargs):
    with io.open(
            join(dirname(__file__), *names),
            encoding=kwargs.get('encoding', 'utf8')
            ) as fh:
        return fh.read()


setup(
    name='taurenmd',
    version='0.7.2',
    license='GNU GPLv2',
    description='Command-line and library interface for analysis routines in Molecular Dynamics',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
        ),
    author='Joao MC Teixeira',
    author_email='joaomcteixeira@gmail.com',
    url='https://github.com/joaomcteixeira/taurenmd',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 1 - Planning',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Environment :: Console',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Utilities',
        ],
    project_urls={
        'Documentation': 'https://taurenmd.readthedocs.io/',
        'Changelog': 'https://taurenmd.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/joaomcteixeira/taurenmd/issues',
        },
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
        ],
    python_requires='>=3, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, !=3.5.*',
    install_requires=[
        'bioplottemplates',
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
        ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
        },
    entry_points={
        'console_scripts': [
            'taurenmd = taurenmd.cli:maincli',
            'taurenmd_angle = taurenmd.cli_angle:maincli',
            'taurenmd_dist = taurenmd.cli_distances:maincli',
            'taurenmd_fext = taurenmd.cli_fext:maincli',
            'taurenmd_imagemol = taurenmd.cli_imagemol:maincli',
            'taurenmd_noSol = tauremd.cli_noSol:maincli',
            'taurenmd_report = taurenmd.cli_report:maincli',
            'taurenmd_rot = taurenmd.cli_rot:maincli',
            'taurenmd_rmsf = taurenmd.cli_rmsf:maincli',
            'taurenmd_rmsd = taurenmd.cli_rmsd:maincli',
            'taurenmd_trajedit = taurenmd.cli_trajedit:maincli',
            ]
        },
    )
