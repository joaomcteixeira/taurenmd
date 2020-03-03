"""
Installation setup.

Read the full installation instructions at:

    https://taurenmd.readthedocs.io/en/latest/installation.html

To install taurenmd run:

    >>> python setup.py

To install taurenmd in developer mode run:

    >>> python setup.py develop

To install taurenmd without ANY of its dependencies:

    >>> python setup.py --no-deps
"""
from __future__ import absolute_import, print_function

import io
import os
import re
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages, setup


install_requires = [
    'bioplottemplates>0.1',
    'pyquaternion',
    ]

if os.getenv('READTHEDOCS'):
    install_requires = []

# support deps
supdeps = [
    'matplotlib>=3,<4',
    'numpy>=1,<2',
    ]

# Molecular Dynamics deps
mddeps = [
    'MDAnalysis>=0.2,<1',
    'mdtraj>=1,<2',
    ]

alldeps = mddeps + supdeps


def _read(*names, **kwargs):
    with io.open(
            join(dirname(__file__), *names),
            encoding=kwargs.get('encoding', 'utf8')
            ) as fh:
        return fh.read()


setup(
    name='taurenmd',
    version='0.8.9',
    license='GNU GPLv2',
    description='A command-line interface for analysis routines in Molecular Dynamics data.',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', _read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', _read('CHANGELOG.rst'))
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
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: Unix',
        'Operating System :: POSIX',
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
        'Molecular Dynamics',
        'Proteins',
        'DNA',
        'RNA',
        'Structural Biology',
        'Molecular Biology',
        'Biochemistry',
        ],
    python_requires='>=3.6, <3.8',
    install_requires=install_requires,
    extras_require={
        'all': alldeps,
        'md': mddeps,
        'sup': supdeps,
        },
    entry_points={
        'console_scripts': [
            'taurenmd = taurenmd.cli:maincli',
            # add bellow, by alphabetical, order your newly developed
            # client, following the nomenclature provided,
            # where NAME is the command line name of your subroutine.
            # 'tmdNAME = taurenmd.cli_NAME:maincli',
            'tmdpangle = taurenmd.cli_pangle:maincli',
            'tmddist = taurenmd.cli_distances:maincli',
            'tmdfext = taurenmd.cli_fext:maincli',
            'tmdimagemol = taurenmd.cli_imagemol:maincli',
            'tmdnosol = tauremd.cli_nosol:maincli',
            'tmdreport = taurenmd.cli_report:maincli',
            'tmdrot = taurenmd.cli_rot:maincli',
            'tmdrmsf = taurenmd.cli_rmsf:maincli',
            'tmdrmsd = taurenmd.cli_rmsd:maincli',
            'tmdtrajedit = taurenmd.cli_trajedit:maincli',
            ]
        },
    )
