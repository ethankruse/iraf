from iraf.utils import file_handler
from iraf import Instrument
# import numpy as np
# import os
# import sys
# from iraf.sys import image_open, image_close

__all__ = ['ccdproc']


def ccdproc(images, output, *, ccdtype='object', noproc=False, fixpix=True,
            overscan=True, trim=True, zerocor=True, darkcor=True, flatcor=True,
            illumcor=False, fringecor=False, readcor=False, scancor=False,
            readaxis='line', fixfile=None, biassec=None, trimsec=None,
            zero=None, dark=None, flat=None, illum=None, fringe=None,
            minreplace=1., scantype='shortscan', nscan=1, interactive=False,
            function='legendre', order=1, sample='*', naverage=1, niterate=1,
            low_reject=3., high_reject=3., grow=0., instrument=None):
    inputs = file_handler(images)
    # XXX: is output required? What happens if you don't give it one?
    outputs = file_handler(output)

    if 0 < len(outputs) != len(inputs):
        raise Exception("Input and output lists do not match")

    # was given a string or something else, so set up the instrument object
    if not isinstance(instrument, Instrument):
        instrument = Instrument(instrument)

    # this allows interactive to be 4 valued (yes, no, always yes, always no)
    # to allow for not prompting for every image
    # set_interactive("", interactive)
