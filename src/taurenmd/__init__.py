# THIS FILE IS EXCLUDED FROM FLAKE8 IN TOX.ini WORKFLOW
r"""
Welcome to
```
_________ _______           _______  _______  _        _______  ______  
\__   __/(  ___  )|\     /|(  ____ )(  ____ \( (    /|(       )(  __  \ 
   ) (   | (   ) || )   ( || (    )|| (    \/|  \  ( || () () || (  \  )
   | |   | (___) || |   | || (____)|| (__    |   \ | || || || || |   ) |
   | |   |  ___  || |   | ||     __)|  __)   | (\ \) || |(_)| || |   | |
   | |   | (   ) || |   | || (\ (   | (      | | \   || |   | || |   ) |
   | |   | )   ( || (___) || ) \ \__| (____/\| )  \  || )   ( || (__/  )
   )_(   |/     \|(_______)|/   \__/(_______/|/    )_)|/     \|(______/ 
                                                                       
```
**A command-line interface for Molecular Dynamics Analysis routines.**

version: {}
"""
# link to logo
# http://patorjk.com/software/taag/#p=display&h=0&f=Epic&t=taurenmd
import logging

from taurenmd.core import Path  # noqa: F401
from taurenmd.logger import DEBUGFILE, LOGFILE

__version__ = '0.8.9'

_BANNER = __doc__.format(__version__)

_DOCUMENTATION = """
Documentation
=============
This project is fully documented at: https://taurenmd.readthedocs.io/
Client documentation is provided with the ``-h`` option, for example:

    >>> taurenmd -h

    or

    >>> taurenmd dist -h

Citing
======
taurenmd uses other scientific libraries to handle and generate
Molecular Dynamics data. You SHOULD by all means cite also the other
libraries when using taurenmd. Please visit our documentation page
for complete details on how to cite properly:

    https://taurenmd.readthedocs.io/en/latest/citing.html
""".format(__version__)

__doc__ = _BANNER + _DOCUMENTATION

_INTERFACE_DESCRIPTION = __doc__ + "Usage\n=====\n\n"

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_db = logging.FileHandler(DEBUGFILE, mode='w')
_db.setLevel(logging.DEBUG)
_db.setFormatter(
    logging.Formatter(
        '%(filename)s:%(name)s:%(funcName)s:%(lineno)d: %(message)s'
        )
    )

_info = logging.FileHandler(LOGFILE, mode='w')
_info.setLevel(logging.INFO)
_info.setFormatter(logging.Formatter('%(message)s'))

_ch = logging.StreamHandler()
_ch.setLevel(logging.INFO)
_ch.setFormatter(logging.Formatter('%(message)s'))

log.addHandler(_db)
log.addHandler(_info)
log.addHandler(_ch)

references = set()
