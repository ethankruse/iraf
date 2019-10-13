"""
Implementation of the IRAF ccdred package and routines.
"""
_imagetypes = "object|zero|dark|flat|illum|fringe|other|" \
              "comp|none|unknown".split('|')

from .instrument_routines import *
from . import utils
from .combine_routines import *
from .ccdproc_routines import *
del instrument_routines
del combine_routines
del ccdproc_routines
