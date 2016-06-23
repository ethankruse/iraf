from iraf._cl import load_task


def _noao():
    return

noao = load_task(_noao, 'noao')

from imred import *
del load_task
