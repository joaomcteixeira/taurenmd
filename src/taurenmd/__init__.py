r"""
Welcome to
_________ _______           _______  _______  _        _______  ______  
\__   __/(  ___  )|\     /|(  ____ )(  ____ \( (    /|(       )(  __  \ 
   ) (   | (   ) || )   ( || (    )|| (    \/|  \  ( || () () || (  \  )
   | |   | (___) || |   | || (____)|| (__    |   \ | || || || || |   ) |
   | |   |  ___  || |   | ||     __)|  __)   | (\ \) || |(_)| || |   | |
   | |   | (   ) || |   | || (\ (   | (      | | \   || |   | || |   ) |
   | |   | )   ( || (___) || ) \ \__| (____/\| )  \  || )   ( || (__/  )
   )_(   |/     \|(_______)|/   \__/(_______/|/    )_)|/     \|(______/ 
                                                                       
**A command-line interface for Molecular Dynamics Analysis routines.**
"""
# link to logo
# http://patorjk.com/software/taag/#p=display&h=0&f=Epic&t=taurenmd

_BANNER = __doc__

_DOCSTRING = """
version: 0.7.2

Documentation
=============
This project is fully documented at: https://taurenmd.readthedocs.io/

Citing
======
taurenmd uses other scientific libraries to handle and generate
Molecular Dynamics data. You SHOULD by all means cite also the other
libraries when using taurenmd. Please visit our documentation page
for complete details on how to cite properly:

    https://taurenmd.readthedocs.io/en/latest/citing.html
"""
import logging

from taurenmd.core import Path
from taurenmd.logger import CMDFILE, DEBUGFILE, LOGFILE

__doc__ += _DOCSTRING

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

__version__ = '0.7.2'

references = set()
