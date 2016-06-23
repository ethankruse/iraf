from iraf._cl import load_task


def _imred():
    return

imred = load_task(_imred, 'imred')

from ccdred import *
del load_task
