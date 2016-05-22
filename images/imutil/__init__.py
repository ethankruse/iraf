from iraf import loadpackage
from imstat import *


def _imutil():
    return

imutil = loadpackage(_imutil, 'imutil')
imstatistics = loadpackage(_imstatistics, 'imstatistics')
del loadpackage
