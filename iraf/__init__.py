try:
    __IRAF_SETUP__
except NameError:
    __IRAF_SETUP__ = False

__version__ = "0.1"
__hdrstring__ = 'AIRAF in Python'

if not __IRAF_SETUP__:
    from . import utils
    from .images import *
    from .noao import *
    from .plot import *
    from . import sys
