from iraf import load_task


def _noao():
    return

noao = load_task(_noao, 'noao')

from imred import *
