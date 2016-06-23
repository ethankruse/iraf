from iraf._cl import load_task


def _images():
    return

images = load_task(_images, 'images')

from imutil import *
del load_task
