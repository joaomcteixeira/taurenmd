import logging
import os
from pathlib import Path as _Path

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_db = logging.FileHandler('.taurenmd.debug', mode='w')
_db.setLevel(logging.DEBUG)
_db.setFormatter(
    logging.Formatter(
        '%(filename)s:%(name)s:%(funcName)s:%(lineno)d: %(message)s'
        )
    )

_ch = logging.StreamHandler()
_ch.setLevel(logging.INFO)
_ch.setFormatter(logging.Formatter('%(message)s'))

log.addHandler(_db)
log.addHandler(_ch)


CMDFILE = '.taurenmd_cmd_register.log'


class Path(type(_Path())):
    """
    Define a common Path to string interface.

    Avoids using os.fspath around libs.
    """
    def str(self):
        """Return string version of Path."""
        return os.fspath(self)
    
    def myparents(self):
        return self.resolve().parents[0]


__version__ = '0.7.2'
