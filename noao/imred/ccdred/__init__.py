from iraf import loadpackage
from combine import *


def _ccdred():
    return

ccdred = loadpackage(_ccdred, 'ccdred')
combine = loadpackage(_combine, 'combine')
