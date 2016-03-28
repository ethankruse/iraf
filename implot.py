from .imstat import loadparams, file_handler
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.widgets import Button, CheckButtons


im_set = {'fig': None, 'ax': None, 'lineplot': True,
          'line': None, 'im': None, 'ndim': None, 'navg': None, 'ncols': None,
          'nlines': None, 'image': None, 'index': None, 'cid': None}


def implot_redraw():
    im_set['fig'].clf()
    ax = im_set['fig'].add_subplot(111)
    im_set['ax'] = ax


def implot_plot():
    if im_set['lineplot']:
        y1 = max(0, min(im_set['nlines'] - 1, im_set['line']))
        y2 = max(1, min(im_set['nlines'], im_set['line'] + im_set['ndim']))
        if im_set['ndim'] == 1:
            yplot = im_set['im']
        else:
            yplot = im_set['im'][:, y1:y2].mean(axis=1)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(im_set['ncols'])
    else:
        x1 = max(0, min(im_set['ncols'] - 1, im_set['line']))
        x2 = max(1, min(im_set['ncols'], im_set['line'] + im_set['ndim']))
        if im_set['ndim'] == 1:
            yplot = im_set['im'][x1:x2]
        else:
            yplot = im_set['im'][x1:x2, :].mean(axis=0)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(im_set['nlines'])

    # if (format[1] == '%')
    #    call strcpy(format, fmt, SZ_FNAME)

    im_set['ax'].plot(xplot, yplot)

    if im_set['ndim'] > 1:
        if im_set['ndim'] > 1:
            if im_set['lineplot']:
                txt = "lines"
                i1 = y1
                i2 = y2
                sz = im_set['nlines']
            else:
                txt = "columns"
                i1 = x1
                i2 = x2
                sz = im_set['ncols']
            title = "Average of {0} {1:d} to {2:d} of {3} in\n{4}"
            title = title.format(txt, i1, i2, sz,
                                 im_set['image'][im_set['index']])
        else:
            if im_set['lineplot']:
                txt = "Line"
                sz = im_set['nlines']
            else:
                txt = "Column"
                sz = im_set['ncols']
            title = "{0} {1:d} of {2} in\n{3}".format(txt, im_set['line'], sz,
                                                      im_set['image'][
                                                          im_set['index']])
    else:
        title = im_set['image'][im_set['index']]

    im_set['ax'].set_title(title)


def implot_keypress(event):
    print 'press', event.key
    if event.key == 'x':
        im_set['ax'].plot(np.random.randn(50) * 500, np.random.randn(50) * 200)
    if event.key == 'q':
        print 'q'
        im_set['ax'].set_xlim(np.random.randn() * 50,
                              np.random.randn() * 50 + 1500)
    plt.draw()

    if event.key == ':':
        im_set['fig'].canvas.mpl_disconnect(im_set['cid'])
        inp = raw_input('What do you want?')
        print inp
        im_set['cid'] = im_set['fig'].canvas.mpl_connect('key_press_event',
                                               implot_keypress)


def test(event):
    im_set['ax'].set_xlim(np.random.randn()*50, np.random.randn()*50 + 1500)
    im_set['ax'].plot(np.random.randn(50) * 500, np.random.randn(50) * 200)
    plt.draw()


def implot(*args, **kwargs):
    params = loadparams(*args, **kwargs)
    im_set['image'] = file_handler(params['image'].value)
    im_set['line'] = params['line'].value
    wcs = params['wcs'].value
    step = params['step'].value
    coords = params['coords'].value
    device = params['device'].value

    logscale = False
    overplot = False
    erase = False
    xnticks = 5
    ynticks = 5

    linestep = 1
    linetype = 1
    color = 1
    p_navg = 1
    im_set['ndim'] = 1

    # we couldn't find any images to plot
    if len(im_set['image']) == 0:
        return

    plt.ion()
    fig, ax = plt.subplots()
    im_set['fig'] = fig
    im_set['ax'] = ax

    # number of images in the list to look through
    nim = len(im_set['image'])
    # which image we're examining now
    im_set['index'] = 0

    # open the image
    try:
        hdulist = fits.open(im_set['image'][im_set['index']])
    except IOError:
        print "Error reading image {0} ...".format(
            im_set['image'][im_set['index']])
        # XXX: go to next image

    im_set['im'] = hdulist[0].data
    if im_set['im'] is None:
        # XXX: call error (1, "image has no pixels")
        pass

    im_set['ncols'] = im_set['im'].shape[0]
    if len(im_set['im'].shape) > 1:
        im_set['nlines'] = im_set['im'].shape[1]
        im_set['ndim'] = 2
    else:
        im_set['nlines'] = 1
        im_set['ndim'] = 1

    if im_set['line'] is None:
        im_set['line'] = max(1,
                             min(im_set['nlines'], (im_set['nlines'] + 1) / 2))

    if step is None or step < 1:
        step = max(1, im_set['nlines'] / 10)

    if not overplot:
        implot_redraw()

    implot_plot()

    fig.subplots_adjust(right=0.8)

    rax = plt.axes([0.85, 0.4, 0.15, 0.15])
    bnext = Button(rax, 'Next')
    nax = plt.axes([0.85, 0.6, 0.15, 0.15])
    check = CheckButtons(nax, ('2 Hz', '4 Hz', '6 Hz'), (False, True, True))
    check.on_clicked(test)
    bnext.on_clicked(test)
    im_set['cid'] = fig.canvas.mpl_connect('key_press_event', implot_keypress)

    ax._button = bnext
    ax._radio = check

    overplot = False

    plt.show(block=False)

    hdulist.close()

    return params
