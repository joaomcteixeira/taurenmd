import logging


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

_db = logging.StreamHandler()
_db.setLevel(logging.DEBUG)
_db.setFormatter('* DEBUG [%(lineno)d]*: %(message)s')

_ch = logging.StreamHandler()
_ch.setLevel(logging.INFO)
_ch.setFormatter('%(message)s')

log.addHandler(_db)
log.addHandler(_ch)

__version__ = '0.0.0'
