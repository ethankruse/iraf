import os

uparam_dir = os.path.join(os.getcwd(), 'iraf_uparam')

from _cl import *

del os


# XXX: somehow give first level packages direct access, e.g. iraf.images, and
# have their package point to the CL one.

# make a 'packages' directory and put everything in there. Then import from
# there using a loadpackage() function.
