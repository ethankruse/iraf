from iraf._cl import load_task
from combine import *


def _ccdred():
    return

ccdred = load_task(_ccdred, 'ccdred')
combine = load_task(_combine, 'combine')
del load_task
