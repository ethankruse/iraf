from iraf import load_task
from combine import *


def _ccdred():
    return

load_task(_ccdred, 'ccdred')
load_task(_combine, 'combine')
