from .imstat import loadparams, file_handler
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.widgets import Button, CheckButtons, Slider
import copy
import functools


im_set_init = {'fig': None, 'ax': None, 'sax': None, 'lineplot': True,
          'line': None, 'im': None, 'ndim': None, 'navg': None, 'ncols': None,
          'nlines': None, 'image': None, 'index': None, 'cid': None,
          'input': '', 'iax': None, 'overplot': False, 'sid': None,
          '_slider': None, 'stat': None, 'xdata': None, 'ydata': None}

bell = '\a'


def implot_plot(fig):
    if not fig.im_set['overplot']:
        fig.im_set['ax'].cla()

    if fig.im_set['lineplot']:
        y1 = max(0, min(fig.im_set['nlines'] - 1, fig.im_set['line']))
        y2 = max(1, min(fig.im_set['nlines'], fig.im_set['line'] + fig.im_set['navg']))
        if fig.im_set['ndim'] == 1:
            yplot = fig.im_set['im']
        else:
            yplot = fig.im_set['im'][:, y1:y2].mean(axis=1)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(fig.im_set['ncols'])
    else:
        x1 = max(0, min(fig.im_set['ncols'] - 1, fig.im_set['line']))
        x2 = max(1, min(fig.im_set['ncols'], fig.im_set['line'] + fig.im_set['navg']))
        if fig.im_set['ndim'] == 1:
            yplot = fig.im_set['im'][x1:x2]
        else:
            yplot = fig.im_set['im'][x1:x2, :].mean(axis=0)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(fig.im_set['nlines'])

    # if (format[1] == '%')
    #    call strcpy(format, fmt, SZ_FNAME)

    fig.im_set['ax'].plot(xplot, yplot)
    fig.im_set['xdata'] = xplot
    fig.im_set['ydata'] = yplot

    if fig.im_set['ndim'] > 1:
        if fig.im_set['navg'] > 1:
            if fig.im_set['lineplot']:
                txt = "lines"
                i1 = y1
                i2 = y2
                sz = fig.im_set['nlines']
            else:
                txt = "columns"
                i1 = x1
                i2 = x2
                sz = fig.im_set['ncols']
            title = "Average of {0} {1:d} to {2:d} of {3} in\n{4}"
            title = title.format(txt, i1, i2-1, sz,
                                 fig.im_set['image'][fig.im_set['index']])
        else:
            if fig.im_set['lineplot']:
                txt = "Line"
                sz = fig.im_set['nlines']
            else:
                txt = "Column"
                sz = fig.im_set['ncols']
            title = "{0} {1:d} of {2} in\n{3}".format(txt, fig.im_set['line'], sz,
                                                      fig.im_set['image'][
                                                          fig.im_set['index']])
    else:
        title = fig.im_set['image'][fig.im_set['index']]

    fig.im_set['ax'].set_title(title)

    fig.im_set['fig'].canvas.draw_idle()


def implot_plot_line(fig):
    was_set = False
    # we're currently plotting columns
    if not fig.im_set['lineplot']:
        was_set = fig.im_set['overplot']
        # turn overplot off
        if fig.im_set['overplot']:
            fig.im_set['_check'].set_active(0)
            fig.im_set['overplot'] = False
        fig.im_set['lineplot'] = True
        # redo the slider
        implot_remove_slider(fig)
        implot_make_slider(fig)

    implot_plot(fig)
    if was_set:
        fig.im_set['_check'].set_active(0)
        fig.im_set['overplot'] = True


def implot_plot_col(fig):
    was_set = False
    # we're currently plotting lines
    if fig.im_set['lineplot']:
        # turn overplot off
        was_set = fig.im_set['overplot']
        if fig.im_set['overplot']:
            fig.im_set['_check'].set_active(0)
            fig.im_set['overplot'] = False

        fig.im_set['lineplot'] = False
        # redo the slider
        implot_remove_slider(fig)
        implot_make_slider(fig)

    implot_plot(fig)
    if was_set:
        fig.im_set['_check'].set_active(0)
        fig.im_set['overplot'] = True


def implot_keypress(event, fig=None):
    # pressed for the second time
    if fig.im_set['stat'] is not None:
        event.key = 's'

    # user input is complete and was one of the extended options
    if event.key == 'enter' and len(fig.im_set['input']):
        full_inp = fig.im_set['input']
        fig.im_set['input'] = ''
        fig.im_set['iax'].set_text('')
        fig.im_set['fig'].canvas.draw_idle()
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
                fig.im_set['navg'] = int(pieces[1])
                if fig.im_set['navg'] <= 0:
                    fig.im_set['navg'] = 1
            elif pieces[0] == 'i' and len(pieces) > 1:
                infile = file_handler(pieces[1])
                if len(infile) > 0:
                    implot_open_image(fig, infile=infile[0])
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
                        fig.im_set['line'] = int(pieces[1])
                        if pieces[0] == 'l':
                            implot_plot_line(fig)
                        else:
                            implot_plot_col(fig)
                    else:
                        print bell
                else:
                    if int(pieces[1]) >= 0 and int(pieces[2]) >= 0:
                        fig.im_set['line'] = min(int(pieces[1]), int(pieces[2]))
                        oldavg = fig.im_set['navg']
                        fig.im_set['navg'] = np.abs(int(pieces[1]) - int(pieces[2])) - 1
                        if pieces[0] == 'l':
                            implot_plot_line(fig)
                        else:
                            implot_plot_col(fig)
                        fig.im_set['navg'] = oldavg
                    else:
                        print bell
            # toggle log y scale on or off
            elif pieces[0] == 'log+':
                fig.im_set['ax'].set_yscale('log')
            elif pieces[0] == 'log-':
                fig.im_set['ax'].set_yscale('linear')

            else:
                print bell


        except ValueError:
            print bell

        return

    # user is continuing their long input
    if len(fig.im_set['input']):
        # ignore keypress events like 'up', 'down', 'shift'.
        if len(event.key) == 1:
            fig.im_set['input'] += event.key

        fig.im_set['iax'].set_text(fig.im_set['input'])
        fig.im_set['fig'].canvas.draw_idle()
        return
    # user is beginning a long input
    if event.key == ':':
        fig.im_set['input'] += ':'
        fig.im_set['iax'].set_text(fig.im_set['input'])
        # XXX: make a global variable of 'previous inputs' that resets when
        # you hit the enter key. Then interpret that.
        # Also make use of the stdout manipulation to be printing what you're
        # typing on the command line.

    # wants to plot lines
    if event.key == 'l':
        # we're currently plotting columns
        if not fig.im_set['lineplot']:
            # which column to plot
            fracthru = fig.im_set['line'] * 1. / fig.im_set['ncols']
            fig.im_set['line'] = int(fig.im_set['nlines'] * fracthru)

            implot_plot_line(fig)

    # wants to plot columns
    if event.key == 'c':
        # we're currently plotting lines
        if fig.im_set['lineplot']:
            if fig.im_set['nlines'] == 1:
                print bell
                return

            # which line to plot
            fracthru = fig.im_set['line'] * 1. / fig.im_set['nlines']
            fig.im_set['line'] = int(fig.im_set['ncols'] * fracthru)

            implot_plot_col(fig)

    # go to previous image in list
    if event.key == 'm':
        if fig.im_set['index'] > 0:
            fig.im_set['index'] -= 1
            implot_open_image(fig)

    # go to next image in list
    if event.key == 'n':
        if fig.im_set['index'] < len(fig.im_set['image']) - 1:
            fig.im_set['index'] += 1
            implot_open_image(fig)

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
        if (event.inaxes is not fig.im_set['ax'] or event.xdata is None or
                    event.ydata is None):
            fig.im_set['stat'] = None
            return

        if fig.im_set['stat'] is None:
            fig.im_set['stat'] = (event.xdata, event.ydata)
            print 'Again'
            return

        # have 2 successive presses, calculate the statistics
        x1 = min(fig.im_set['stat'][0], event.xdata)
        x2 = max(fig.im_set['stat'][0], event.xdata)

        region = np.where((fig.im_set['xdata'] >= x1) & (fig.im_set['xdata'] <= x2))[0]
        if len(region) == 0:
            ind1 = np.abs(fig.im_set['xdata'] - x1).argmin()
            if ind1 != len(fig.im_set['xdata']):
                region = np.array([ind1, ind1+1])
            else:
                region = np.array([ind1, ind1 - 1])
        mean = fig.im_set['ydata'][region].mean()
        std = fig.im_set['ydata'][region].std()
        sum = fig.im_set['ydata'][region].sum()
        median = np.median(fig.im_set['ydata'][region])
        outstr = "Median={0:g}, mean={1:g}, rms={2:g}, sum={3:g}, npix={4:d}"
        print outstr.format(median, mean, std, sum, len(region))
        fig.im_set['stat'] = None

    if event.key == ' ':
        if event.xdata is not None:
            pixel = np.abs(event.xdata - fig.im_set['xdata']).argmin()
            txt = 'Pixel = [{0}, {1}] Value = {2}'
            if fig.im_set['lineplot']:
                txt = txt.format(fig.im_set['line'], pixel, fig.im_set['ydata'][pixel])
            else:
                txt = txt.format(pixel, fig.im_set['line'], fig.im_set['ydata'][pixel])
            print txt

    if event.key == '?':
        # XXX: implement pagefiles
        """
        # Print command summary.
        call
        gpagefile(gp, KEYSFILE, "implot cursor commands")
        """
        pass

    fig.im_set['fig'].canvas.draw_idle()


def button_click(label=None, fig=None):
    if label == 'Overplot':
        if not fig.im_set['overplot']:
            fig.im_set['overplot'] = True
        else:
            fig.im_set['overplot'] = False


def implot_change_line(value=None, fig=None):
    fig.im_set['line'] = int(value)
    implot_plot(fig)
    pass


def implot_open_image(fig, infile=None):
    if infile is None:
        infile = fig.im_set['image'][fig.im_set['index']]

    # open the image
    try:
        hdulist = fits.open(infile)
    except IOError:
        print "Error reading image {0} ...".format(
            fig.im_set['image'][fig.im_set['index']])
        return
        # XXX: go to next image

    fig.im_set['im'] = hdulist[0].data
    hdulist.close()

    if fig.im_set['im'] is None:
        # XXX: call error (1, "image has no pixels")
        pass

    fig.im_set['ncols'] = fig.im_set['im'].shape[0]
    if len(fig.im_set['im'].shape) > 1:
        fig.im_set['nlines'] = fig.im_set['im'].shape[1]
        fig.im_set['ndim'] = 2
    else:
        fig.im_set['nlines'] = 1
        fig.im_set['ndim'] = 1

    if fig.im_set['line'] is None:
        fig.im_set['line'] = max(0, min(fig.im_set['nlines'],
                                    (fig.im_set['nlines'] + 1) / 2) - 1)
        fig.im_set['lineinit'] = fig.im_set['line']
        fig.im_set['colinit'] = max(0,
                                min(fig.im_set['ncols'], (fig.im_set['ncols'] + 1) / 2) - 1)
    # redo the slider
    implot_remove_slider(fig)
    implot_make_slider(fig)

    implot_plot(fig)


def implot_remove_slider(fig):
    if fig.im_set['_slider'] is None:
        return

    fig.im_set['_slider'].disconnect(fig.im_set['sid'])
    fig.im_set['fig'].delaxes(fig.im_set['sax'])
    fig.im_set['sax'] = None
    fig.im_set['sid'] = None
    fig.im_set['_slider'] = None


def implot_make_slider(fig):
    if fig.im_set['lineplot']:
        txt = 'Line'
    else:
        txt = 'Column'
    # slider can't adjust bounds, so make sure it has room for the largest
    # dimension, then just prevent the smaller dimension from having higher
    # values
    if fig.im_set['lineplot']:
        smax = fig.im_set['nlines'] - 1
    else:
        smax = fig.im_set['ncols'] - 1
    # smax = max(fig.im_set['ncols'], fig.im_set['nlines']) - 1
    fig.im_set['sax'] = plt.axes([0.83, 0.6, 0.15, 0.05], zorder=1)
    slider = Slider(fig.im_set['sax'], txt, 0, smax, valinit=fig.im_set['line'],
                    valfmt='%.0f')
    partial = functools.partial(implot_change_line, fig=fig)
    fig.im_set['sid'] = slider.on_changed(partial)
    slider.vline.set_visible(False)
    fig.im_set['sax'].legend(loc='top')
    # slider.poly.set_fc('r')
    slider.label.set_x(0.5)
    slider.label.set_y(1.4)
    slider.label.set_ha('center')
    slider.valtext.set_x(0.5)
    slider.valtext.set_y(-0.5)
    slider.valtext.set_ha('center')

    fig.im_set['_slider'] = slider


def implot(*args, **kwargs):
    # XXX: where does this go?
    # Disable default Matplotlib shortcut keys:
    keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
    for key in keymaps:
        plt.rcParams[key] = ''

    params = loadparams(*args, **kwargs)
    images = file_handler(params['image'].value)

    # we couldn't find any images to plot
    if len(images) == 0:
        return

    plt.ion()
    fig, ax = plt.subplots()

    fig.im_set = copy.deepcopy(im_set_init)

    fig.im_set['fig'] = fig
    fig.im_set['ax'] = ax

    fig.im_set['image'] = images
    fig.im_set['line'] = params['line'].value
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
    fig.im_set['navg'] = 1

    # number of images in the list to look through
    nim = len(fig.im_set['image'])
    # which image we're examining now
    fig.im_set['index'] = 0

    # open the image and plot it
    implot_open_image(fig)

    partial = functools.partial(implot_keypress, fig=fig)
    fig.im_set['cid'] = fig.canvas.mpl_connect('key_press_event', partial)

    fig.subplots_adjust(right=0.8)

    inptxt = fig.text(0.02, 0.02, '')
    fig.im_set['iax'] = inptxt

    nax = plt.axes([0.85, 0.65, 0.2, 0.25], zorder=-1)
    check = CheckButtons(nax, ('Overplot',), (fig.im_set['overplot'], ))
    partial = functools.partial(button_click, fig=fig)
    check.on_clicked(partial)

    fig.im_set['_check'] = check
    nax.set_axis_off()

    plt.show(block=False)

    return nax
