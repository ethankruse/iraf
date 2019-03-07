"""
AIRAF: Astronomical Image Reduction and Analysis Facility
======

TODO: write overview of package and docs.

"""

# Taken from numpy/astropy etc.
# This is only set if we're running from the setup.py file.
try:
    __IRAF_SETUP__
except NameError:
    __IRAF_SETUP__ = False

__version__ = "0.1"
# the package name that gets inserted into images created by the package
# in the header value 'origin'
__hdrstring__ = f"AIRAF v{__version__} for Python"

if not __IRAF_SETUP__:
    from . import utils
    from . import sys
    from .images import *
    del images
    from .noao import *
    del noao
    from . import plot
