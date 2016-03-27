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

    fig, ax = plt.subplots()

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
        ndim = 2
    else:
        nlines = 1
        ndim = 1

    if line is None:
        line = max(1, min(nlines, (nlines + 1) / 2))

    if step is None or step < 1:
        step = max(1, nlines / 10)

    npix = max(ncols, nlines)

    if not overplot:
        plt.clf()

    if lineplot:
        y1 = max(0, min(nlines - 1, line))
        y2 = max(1, min(nlines, line + navg))
        if ndim == 1:
            yplot = im
        else:
            yplot = im[:, y1:y2].mean(axis=1)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(ncols)
    else:
        x1 = max(0, min(ncols - 1, line))
        x2 = max(1, min(ncols, line + navg))
        if ndim == 1:
            yplot = im[x1:x2]
        else:
            yplot = im[x1:x2, :].mean(axis=0)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(nlines)

    # if (format[1] == '%')
    #    call strcpy(format, fmt, SZ_FNAME)

    plt.plot(xplot, yplot)

    if ndim > 1:
        if navg > 1:
            if lineplot:
                txt = "lines"
                i1 = y1
                i2 = y2
                sz = nlines
            else:
                txt = "columns"
                i1 = x1
                i2 = x2
                sz = ncols
            title = "Average of {0} {1:d} to {2:d} of {3} in\n{4}".format(txt, i1, i2, sz, image[index])
        else:
            if lineplot:
                txt = "Line"
                sz = nlines
            else:
                txt = "Column"
                sz = ncols
            title = "{0} {1:d} of {2} in\n{3}".format(txt, line, sz, image[index])
    else:
        title = image[index]

    plt.title(title)

    overplot = False

    """
    while True:
        def press(event):
            print('press', event.key)
            import sys
            sys.stdout.flush()
            if event.key == 'x':
                visible = xl.get_visible()
                xl.set_visible(not visible)
                fig.canvas.draw()

        fig.canvas.mpl_connect('key_press_event', press)

        ax.plot(np.random.rand(12), np.random.rand(12), 'go')
        xl = ax.set_xlabel('easy come, easy go')

    """

    hdulist.close()

    return params
