import os

uparam_dir = os.path.join(os.getcwd(), '.iraf_uparam')

from _cl import *

# this should be iraf.ccdred.instrument
instrument = None
logfile = os.path.join(os.getcwd(), 'logfile')
ssfile = os.path.join(os.getcwd(), 'subsets')

del os
