from iraf import load_task
from imstat import *


def _imutil():
    return

load_task(_imutil, 'imutil')
load_task(_imstatistics, 'imstatistics')
del load_task
