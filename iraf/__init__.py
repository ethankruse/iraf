# Taken from numpy/astropy etc.
# This is only set if we're running from the setup.py file.
try:
    __IRAF_SETUP__
except NameError:
    __IRAF_SETUP__ = False

__version__ = "0.1"
# the package name that gets inserted into images created by the package
# in the header value 'origin'
__hdrstring__ = "AIRAF in Python"

if not __IRAF_SETUP__:
    from . import utils
    from .images import *
    from .noao import *
    from .plot import *
    from . import sys
