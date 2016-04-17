import os
from cl import *

uparam_dir = os.path.join(os.getcwd(), '/.uparam')

# this should be iraf.ccdred.instrument
instrument = None
logfile = os.path.join(os.getcwd(), 'logfile')
ssfile = os.path.join(os.getcwd(), 'subsets')

del os
