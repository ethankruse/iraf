from .imstat import loadparams, file_handler
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.widgets import Button, CheckButtons, Slider


im_set = {'fig': None, 'ax': None, 'sax': None, 'lineplot': True,
          'line': None, 'im': None, 'ndim': None, 'navg': None, 'ncols': None,
          'nlines': None, 'image': None, 'index': None, 'cid': None,
          'input': '', 'iax': None, 'overplot': False, 'sid': None,
          '_slider': None, 'stat': None, 'xdata': None, 'ydata': None}

bell = '\a'

def implot_plot():
    if not im_set['overplot']:
        im_set['ax'].cla()

    if im_set['lineplot']:
        y1 = max(0, min(im_set['nlines'] - 1, im_set['line']))
        y2 = max(1, min(im_set['nlines'], im_set['line'] + im_set['navg']))
        if im_set['ndim'] == 1:
            yplot = im_set['im']
        else:
            yplot = im_set['im'][:, y1:y2].mean(axis=1)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(im_set['ncols'])
    else:
        x1 = max(0, min(im_set['ncols'] - 1, im_set['line']))
        x2 = max(1, min(im_set['ncols'], im_set['line'] + im_set['navg']))
        if im_set['ndim'] == 1:
            yplot = im_set['im'][x1:x2]
        else:
            yplot = im_set['im'][x1:x2, :].mean(axis=0)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(im_set['nlines'])

    # if (format[1] == '%')
    #    call strcpy(format, fmt, SZ_FNAME)

    im_set['ax'].plot(xplot, yplot)
    im_set['xdata'] = xplot
    im_set['ydata'] = yplot

    if im_set['ndim'] > 1:
        if im_set['navg'] > 1:
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
            title = title.format(txt, i1, i2-1, sz,
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

    im_set['fig'].canvas.draw_idle()


def implot_plot_line():
    was_set = False
    # we're currently plotting columns
    if not im_set['lineplot']:
        was_set = im_set['overplot']
        # turn overplot off
        if im_set['overplot']:
            im_set['_check'].set_active(0)
            im_set['overplot'] = False
        im_set['lineplot'] = True
        # redo the slider
        implot_remove_slider()
        implot_make_slider()

    implot_plot()
    if was_set:
        im_set['_check'].set_active(0)
        im_set['overplot'] = True


def implot_plot_col():
    was_set = False
    # we're currently plotting lines
    if im_set['lineplot']:
        # turn overplot off
        was_set = im_set['overplot']
        if im_set['overplot']:
            im_set['_check'].set_active(0)
            im_set['overplot'] = False

        im_set['lineplot'] = False
        # redo the slider
        implot_remove_slider()
        implot_make_slider()

    implot_plot()
    if was_set:
        im_set['_check'].set_active(0)
        im_set['overplot'] = True


def implot_keypress(event):
    # pressed for the second time
    if im_set['stat'] is not None:
        event.key = 's'

    # user input is complete and was one of the extended options
    if event.key == 'enter' and len(im_set['input']):
        full_inp = im_set['input']
        im_set['input'] = ''
        im_set['iax'].set_text('')
        im_set['fig'].canvas.draw_idle()
        # nothing input
        if len(full_inp.strip()) == 1:
            return
        # print out for the user's history
        print full_inp
        # ignore the leading :
        full_inp = full_inp[1:]
        pieces = full_inp.strip().split()
        try:
            # change the number of lines/columns to average together
            if pieces[0] == 'a' and len(pieces) > 1:
                im_set['navg'] = int(pieces[1])
                if im_set['navg'] <= 0:
                    im_set['navg'] = 1
            elif pieces[0] == 'i' and len(pieces) > 1:
                infile = file_handler(pieces[1])
                if len(infile) > 0:
                    implot_open_image(infile[0])
            elif pieces[0] == 'w' and len(pieces) > 1:
                pass
                # XXX: Change wcs type.
            elif pieces[0] == 'f' and len(pieces) > 1:
                pass
                # XXX: Change label format.
            # plot a selected line/column or average between 2 lines/columns
            elif (pieces[0] == 'l' or pieces[0] == 'c') and len(pieces) > 1:
                if len(pieces) == 2:
                    if int(pieces[1]) >= 0:
                        im_set['line'] = int(pieces[1])
                        if pieces[0] == 'l':
                            implot_plot_line()
                        else:
                            implot_plot_col()
                    else:
                        print bell
                else:
                    if int(pieces[1]) >= 0 and int(pieces[2]) >= 0:
                        im_set['line'] = min(int(pieces[1]), int(pieces[2]))
                        im_set['navg'] = np.abs(int(pieces[1]) - int(pieces[2])) - 1
                        if pieces[0] == 'l':
                            implot_plot_line()
                        else:
                            implot_plot_col()
                    else:
                        print bell
            # toggle log y scale on or off
            elif pieces[0] == 'log+':
                im_set['ax'].set_yscale('log')
            elif pieces[0] == 'log-':
                im_set['ax'].set_yscale('linear')

            else:
                print bell


        except ValueError:
            print bell

        return

    # user is continuing their long input
    if len(im_set['input']):
        # ignore keypress events like 'up', 'down', 'shift'.
        if len(event.key) == 1:
            im_set['input'] += event.key

        im_set['iax'].set_text(im_set['input'])
        im_set['fig'].canvas.draw_idle()
        return
    # user is beginning a long input
    if event.key == ':':
        im_set['input'] += ':'
        im_set['iax'].set_text(im_set['input'])
        # XXX: make a global variable of 'previous inputs' that resets when
        # you hit the enter key. Then interpret that.
        # Also make use of the stdout manipulation to be printing what you're
        # typing on the command line.

    # wants to plot lines
    if event.key == 'l':
        # we're currently plotting columns
        if not im_set['lineplot']:
            # which column to plot
            fracthru = im_set['line'] * 1. / im_set['ncols']
            im_set['line'] = int(im_set['nlines'] * fracthru)

            implot_plot_line()

    # wants to plot columns
    if event.key == 'c':
        # we're currently plotting lines
        if im_set['lineplot']:
            if im_set['nlines'] == 1:
                print bell
                return

            # which line to plot
            fracthru = im_set['line'] * 1. / im_set['nlines']
            im_set['line'] = int(im_set['ncols'] * fracthru)

            implot_plot_col()

    # go to previous image in list
    if event.key == 'm':
        if im_set['index'] > 0:
            im_set['index'] -= 1
            implot_open_image()

    # go to next image in list
    if event.key == 'n':
        if im_set['index'] < len(im_set['image']) - 1:
            im_set['index'] += 1
            implot_open_image()

    if event.key == 'p':
        # XXX: implement impprofile.x
        """
        case 'p':
        # Profile analysis.
        x1 = x
        y1 = y
        call printf (again)
        if (clgcur ("gcur", x2, y2, wcs, key, command,
            SZ_FNAME) == EOF)
            next

        call imp_profile (gp, Memr[xnew], Memr[ynew], npix,
            x1, y1, x2, y2, sl, sline)
        call printf (Memc[sl_getstr(sl,sline)])
        """
        pass

    # print statistics on a region
    if event.key == 's':
        if (event.inaxes is not im_set['ax'] or event.xdata is None or
                    event.ydata is None):
            im_set['stat'] = None
            return

        if im_set['stat'] is None:
            im_set['stat'] = (event.xdata, event.ydata)
            print 'Again'
            return

        # have 2 successive presses, calculate the statistics
        x1 = min(im_set['stat'][0], event.xdata)
        x2 = max(im_set['stat'][0], event.xdata)

        region = np.where((im_set['xdata'] >= x1) & (im_set['xdata'] <= x2))[0]
        if len(region) == 0:
            ind1 = np.abs(im_set['xdata'] - x1).argmin()
            if ind1 != len(im_set['xdata']):
                region = np.array([ind1, ind1+1])
            else:
                region = np.array([ind1, ind1 - 1])
        mean = im_set['ydata'][region].mean()
        std = im_set['ydata'][region].std()
        sum = im_set['ydata'][region].sum()
        median = np.median(im_set['ydata'][region])
        outstr = "Median={0:g}, mean={1:g}, rms={2:g}, sum={3:g}, npix={4:d}"
        print outstr.format(median, mean, std, sum, len(region))
        im_set['stat'] = None

    if event.key == ' ':
        if event.xdata is not None:
            pixel = np.abs(event.xdata - im_set['xdata']).argmin()
            txt = 'Pixel = [{0}, {1}] Value = {2}'
            if im_set['lineplot']:
                txt = txt.format(im_set['line'], pixel, im_set['ydata'][pixel])
            else:
                txt = txt.format(pixel, im_set['line'], im_set['ydata'][pixel])
            print txt

    if event.key == '?':
        # XXX: implement pagefiles
        """
        # Print command summary.
        call
        gpagefile(gp, KEYSFILE, "implot cursor commands")
        """
        pass

    im_set['fig'].canvas.draw_idle()


def test(label):
    if label == 'Overplot':
        if not im_set['overplot']:
            im_set['overplot'] = True
        else:
            im_set['overplot'] = False


def implot_change_line(value):
    im_set['line'] = int(value)
    implot_plot()
    pass


def implot_open_image(infile=None):
    if infile is None:
        infile = im_set['image'][im_set['index']]

    # open the image
    try:
        hdulist = fits.open(infile)
    except IOError:
        print "Error reading image {0} ...".format(
            im_set['image'][im_set['index']])
        return
        # XXX: go to next image

    im_set['im'] = hdulist[0].data
    hdulist.close()

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
        im_set['line'] = max(0, min(im_set['nlines'],
                                    (im_set['nlines'] + 1) / 2) - 1)
        im_set['lineinit'] = im_set['line']
        im_set['colinit'] = max(0,
                                min(im_set['ncols'], (im_set['ncols'] + 1) / 2) - 1)
    # redo the slider
    implot_remove_slider()
    implot_make_slider()

    implot_plot()


def implot_remove_slider():
    if im_set['_slider'] is None:
        return

    im_set['_slider'].disconnect(im_set['sid'])
    im_set['fig'].delaxes(im_set['sax'])
    im_set['sax'] = None
    im_set['sid'] = None
    im_set['_slider'] = None


def implot_make_slider():
    if im_set['lineplot']:
        txt = 'Line'
    else:
        txt = 'Column'
    # slider can't adjust bounds, so make sure it has room for the largest
    # dimension, then just prevent the smaller dimension from having higher
    # values
    if im_set['lineplot']:
        smax = im_set['nlines'] - 1
    else:
        smax = im_set['ncols'] - 1
    # smax = max(im_set['ncols'], im_set['nlines']) - 1
    im_set['sax'] = plt.axes([0.83, 0.6, 0.15, 0.05], zorder=1)
    slider = Slider(im_set['sax'], txt, 0, smax, valinit=im_set['line'],
                    valfmt='%.0f')
    im_set['sid'] = slider.on_changed(implot_change_line)
    slider.vline.set_visible(False)
    im_set['sax'].legend(loc='top')
    # slider.poly.set_fc('r')
    slider.label.set_x(0.5)
    slider.label.set_y(1.4)
    slider.label.set_ha('center')
    slider.valtext.set_x(0.5)
    slider.valtext.set_y(-0.5)
    slider.valtext.set_ha('center')

    im_set['_slider'] = slider


def implot(*args, **kwargs):
    # XXX: where does this go?
    # Disable default Matplotlib shortcut keys:
    keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
    for key in keymaps:
        plt.rcParams[key] = ''

    params = loadparams(*args, **kwargs)
    im_set['image'] = file_handler(params['image'].value)
    im_set['line'] = params['line'].value
    wcs = params['wcs'].value
    step = params['step'].value
    coords = params['coords'].value
    device = params['device'].value

    logscale = False
    erase = False
    xnticks = 5
    ynticks = 5

    linestep = 1
    linetype = 1
    color = 1
    p_navg = 1
    im_set['navg'] = 1

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

    # open the image and plot it
    implot_open_image()

    im_set['cid'] = fig.canvas.mpl_connect('key_press_event', implot_keypress)

    fig.subplots_adjust(right=0.8)

    inptxt = fig.text(0.02, 0.02, '')
    im_set['iax'] = inptxt

    nax = plt.axes([0.85, 0.65, 0.2, 0.25], zorder=-5)
    check = CheckButtons(nax, ('Overplot',), (im_set['overplot'], ))
    check.on_clicked(test)

    im_set['_check'] = check
    nax.set_axis_off()


    plt.show(block=False)



    return nax
