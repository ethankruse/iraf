import os
uparam_dir = os.path.join(os.getcwd(), 'iraf_uparam')
del os

from _cl import load_task, cl

from images import *
from noao import *
from plot import *
