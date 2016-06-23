import os
uparam_dir = os.path.join(os.getcwd(), 'iraf_uparam')
del os

from _cl import cl

from images import *
from noao import *
from plot import *
