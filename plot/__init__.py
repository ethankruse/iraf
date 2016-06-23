from iraf._cl import load_task
from implot import *


def _plot():
    return

plot = load_task(_plot, 'plot')
implot = load_task(_implot, 'implot')
del load_task
