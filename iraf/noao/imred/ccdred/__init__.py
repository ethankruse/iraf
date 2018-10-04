# XXX: include none and unknown as options as well?
imagetypes = "object|zero|dark|flat|illum|fringe|other|comp".split('|')

from .instruments import *
from .combine import *
from .ccdproc import *
