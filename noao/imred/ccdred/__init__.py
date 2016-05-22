from iraf import loadparams
from combine import *


def _ccdred():
    return

ccdred = loadparams(_ccdred, 'ccdred')
combine = loadparams(_combine, 'combine')
