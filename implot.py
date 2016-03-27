from .imstat import loadparams, file_handler
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def implot(*args, **kwargs):
    params = loadparams(*args, **kwargs)
    image = file_handler(params['image'].value)
    line = params['line'].value
    wcs = params['wcs'].value
    step = params['step'].value
    coords = params['coords'].value
    device = params['device'].value

    logscale = False
    overplot = False
    lineplot = True
    erase = False
    xnticks = 5
    ynticks = 5

    linestep = 1
    linetype = 1
    color = 1
    p_navg = 1
    navg = 1

    # we couldn't find any images to plot
    if len(image) == 0:
        return

    plt.figure()

    # number of images in the list to look through
    nim = len(image)
    # which image we're examining now
    index = 0

    # while True:

    # open the image
    try:
        hdulist = fits.open(image[index])
    except IOError:
        print "Error reading image {0} ...".format(image[index])
        # XXX: go to next image

    im = hdulist[0].data
    if im is None:
        # XXX: call error (1, "image has no pixels")
        pass

    ncols = im.shape[0]
    if len(im.shape) > 1:
        nlines = im.shape[1]
    else:
        nlines = 1

    if line is None:
        line = max(1, min(nlines, (nlines + 1) / 2))

    if step is None or step < 1:
        step = max(1, nlines / 10)

    npix = max(ncols, nlines)

    if not overplot:
        plt.clf()



    return params
