"""
Implementation of the IRAF ccdred package and routines.
"""
imagetypes = "object|zero|dark|flat|illum|fringe|other|" \
             "comp|none|unknown".split('|')

from .instruments import *
from .combine import *
from .ccdproc import *
