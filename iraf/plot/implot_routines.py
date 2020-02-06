import copy
import functools

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import CheckButtons, Slider

from iraf.sys import image_open
from iraf.utils import file_handler

__all__ = ['implot']

# all the settings the figure needs to access and their initial values
im_set_init = {'fig': None, 'ax': None, 'sax': None, 'lineplot': True,
               'line': None, 'im': None, 'ndim': None, 'navg': 1,
               'ncols': None, 'helptxt': None, 'wcs': None,
               'nlines': None, 'image': None, 'index': 0, 'cid': None,
               'input': '', 'iax': None, 'overplot': False, 'sid': None,
               '_slider': None, 'stat': None, 'xdata': None, 'ydata': None,
               '_check': None}

bell = '\a'


def implot_plot(fig):
    # reset if we're not overplotting
    if not fig.im_set['overplot']:
        fig.im_set['ax'].cla()
    # plot lines
    if fig.im_set['lineplot']:
        i1 = max(0, min(fig.im_set['nlines'] - 1,
                        fig.im_set['line'] - fig.im_set['navg'] // 2))
        i2 = max(1, min(fig.im_set['nlines'],
                        i1 + fig.im_set['navg']))
        if fig.im_set['ndim'] == 1:
            yplot = fig.im_set['im']
        else:
            yplot = fig.im_set['im'][:, i1:i2].mean(axis=1)
        # XXX: need to add in the WCS transformation here
        # from IRAF:
        # call plt_wcs (im, mw, ct, 1, Memr[axvals], real(x1), real(x2),
        # Memr[x], nx, xlabel, format, SZ_FNAME)
        xplot = np.arange(fig.im_set['ncols'])
    # plot columns
    else:
        i1 = max(0, min(fig.im_set['ncols'] - 1,
                        fig.im_set['line'] - fig.im_set['navg'] // 2))
        i2 = max(1, min(fig.im_set['ncols'],
                        i1 + fig.im_set['navg']))
        # XXX: is 1-D possible here?
        if fig.im_set['ndim'] == 1:
            yplot = fig.im_set['im'][i1:i2]
        else:
            yplot = fig.im_set['im'][i1:i2, :].mean(axis=0)
        # XXX: need to add in the WCS transformation here
        xplot = np.arange(fig.im_set['nlines'])

    fig.im_set['ax'].plot(xplot, yplot)
    fig.im_set['xdata'] = xplot
    fig.im_set['ydata'] = yplot

    if fig.im_set['ndim'] > 1:
        if fig.im_set['lineplot']:
            xlab = 'Column'
        else:
            xlab = 'Line'
        fig.im_set['ax'].set_xlabel(xlab)

        if fig.im_set['navg'] > 1:
            if fig.im_set['lineplot']:
                txt = "lines"
                sz = fig.im_set['nlines'] - 1
            else:
                txt = "columns"
                sz = fig.im_set['ncols'] - 1
            title = "Average of {0} {1:d} to {2:d} of {3} in\n{4}"
            title = title.format(txt, i1, i2 - 1, sz,
                                 fig.im_set['image'][fig.im_set['index']])
        else:
            if fig.im_set['lineplot']:
                txt = "Line"
                sz = fig.im_set['nlines'] - 1
            else:
                txt = "Column"
                sz = fig.im_set['ncols'] - 1
            title = "{0} {1:d} of {2} in\n{3}".format(txt, fig.im_set['line'],
                                                      sz,
                                                      fig.im_set['image'][
                                                          fig.im_set['index']])
    else:
        title = fig.im_set['image'][fig.im_set['index']]

    fig.im_set['ax'].set_title(title)

    fig.im_set['fig'].canvas.draw()


def implot_plot_line(fig, draw=True):
    was_set = False
    # we're currently plotting columns
    if not fig.im_set['lineplot']:
        was_set = fig.im_set['overplot']
        # turn overplot off
        if fig.im_set['overplot']:
            # toggle overplot
            fig.im_set['_check'].set_active(0)
        fig.im_set['lineplot'] = True
        # redo the slider
        implot_remove_slider(fig)
        implot_make_slider(fig)
    if draw:
        implot_plot(fig)
    elif not fig.im_set['overplot']:
        fig.im_set['ax'].cla()
    if was_set:
        # toggle overplot
        fig.im_set['_check'].set_active(0)


def implot_plot_col(fig, draw=True):
    was_set = False
    # we're currently plotting lines
    if fig.im_set['lineplot']:
        # turn overplot off
        was_set = fig.im_set['overplot']
        if fig.im_set['overplot']:
            # toggle overplot
            fig.im_set['_check'].set_active(0)

        fig.im_set['lineplot'] = False
        # redo the slider
        implot_remove_slider(fig)
        implot_make_slider(fig)

    if draw:
        implot_plot(fig)
    elif not fig.im_set['overplot']:
        fig.im_set['ax'].cla()

    if was_set:
        # toggle overplot
        fig.im_set['_check'].set_active(0)


def implot_keypress(event, fig=None):
    """

    Parameters
    ----------
    event
    fig

    Returns
    -------

    """
    # reset stat checking if the user didn't press s
    if event.key != 's':
        fig.im_set['stat'] = None

    helpstr = """
    Implot Command Summary:

    c       plot columns
    l       plot lines
    m       go to the previous image in the input list
    n       go to the next image in the input list
    p       measure profile (mark region and background with 2 pos)
    s       print statistics on a region
    <space> print coordinates and data value
    
    In addition to the above keystrokes, the following ':' escapes
    are recognized by the program:
    
    :a N           set number of lines or columns to be averaged
    :c M [N]       plot column[s] M [to N, inclusive]
    :i image       open a different image
    :l M [N]       plot line[s] M [to N, inclusive]
    :w wcstype     set wcs type (logical|physical|world)
    """

    # user is continuing their long input
    # just add to the input string and return, bypassing single
    # character commands below
    if len(fig.im_set['input']) and event.key != 'enter':
        # ignore keypress events like 'up', 'down', 'shift'.
        if len(event.key) == 1:
            fig.im_set['input'] += event.key

        # let them delete things
        if event.key == 'backspace' or event.key == 'delete':
            fig.im_set['input'] = fig.im_set['input'][:-1]
            if len(fig.im_set['input']) == 0:
                fig.im_set['helptxt'].set_visible(True)

        fig.im_set['iax'].set_text(fig.im_set['input'])
        fig.im_set['fig'].canvas.draw()
        return

    # user input is complete and was one of the extended options
    if event.key == 'enter' and len(fig.im_set['input']):
        full_inp = fig.im_set['input']
        # reset the input
        fig.im_set['input'] = ''
        fig.im_set['iax'].set_text('')
        fig.im_set['helptxt'].set_visible(True)
        fig.im_set['fig'].canvas.draw()
        # nothing input
        if len(full_inp.strip()) == 1:
            return
        # print out for the user's history
        print(full_inp)
        # ignore the leading ':'
        full_inp = full_inp[1:]
        pieces = full_inp.strip().split()
        try:
            # change the number of lines/columns to average together
            if pieces[0] == 'a' and len(pieces) > 1:
                fig.im_set['navg'] = int(pieces[1])
                if fig.im_set['navg'] <= 0:
                    fig.im_set['navg'] = 1
                # make the plot
                implot_plot(fig)

            # open another image
            elif pieces[0] == 'i' and len(pieces) > 1:
                implot_open_image(fig, infile=pieces[1])

            # XXX: Change wcs type.
            elif pieces[0] == 'w' and len(pieces) > 1:
                print('WCS types not implemented yet.')

            # plot a selected line/column or average between 2 lines/columns
            elif (pieces[0] == 'l' or pieces[0] == 'c') and len(pieces) > 1:
                if len(pieces) == 2:
                    if int(pieces[1]) >= 0:
                        fig.im_set['line'] = int(pieces[1])
                        if pieces[0] == 'l':
                            implot_plot_line(fig, draw=False)
                        else:
                            implot_plot_col(fig, draw=False)
                        fig.im_set['_slider'].set_val(fig.im_set['line'])
                    else:
                        print(bell)
                else:
                    if int(pieces[1]) >= 0 and int(pieces[2]) >= 0:
                        i1 = int(pieces[1])
                        i2 = int(pieces[2])
                        # inclusive range
                        diff = np.abs(i1 - i2) + 1
                        # don't use this as the new number of lines to average
                        oldavg = fig.im_set['navg']
                        fig.im_set['navg'] = diff
                        fig.im_set['line'] = min(i1, i2) + diff // 2

                        if pieces[0] == 'l':
                            implot_plot_line(fig, draw=False)
                        else:
                            implot_plot_col(fig, draw=False)
                        fig.im_set['_slider'].set_val(fig.im_set['line'])
                        # reset the number of lines to average
                        fig.im_set['navg'] = oldavg
                    else:
                        print(bell)
            else:
                print(bell)

        except ValueError:
            print(bell)

        return

    # user is beginning a long input
    if event.key == ':':
        fig.im_set['helptxt'].set_visible(False)
        fig.im_set['input'] += ':'
        fig.im_set['iax'].set_text(fig.im_set['input'])

    # wants to plot lines
    # XXX: switching between l/c is not determinative. So it will
    # gradually shift down rather than bouncing between the same 2.
    if event.key == 'l':
        # we're currently plotting columns
        if not fig.im_set['lineplot']:
            # which column to plot
            fracthru = fig.im_set['line'] / fig.im_set['ncols']
            fig.im_set['line'] = int(fig.im_set['nlines'] * fracthru)

            implot_plot_line(fig)

    # wants to plot columns
    if event.key == 'c':
        # we're currently plotting lines
        if fig.im_set['lineplot']:
            if fig.im_set['nlines'] == 1:
                print(bell)
                return

            # which line to plot
            fracthru = fig.im_set['line'] / fig.im_set['nlines']
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
        print("Profile function not currently implemented.")
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

    # print statistics on a region
    if event.key == 's':
        if (event.inaxes is not fig.im_set['ax'] or event.xdata is None or
                event.ydata is None):
            fig.im_set['stat'] = None
            return

        if fig.im_set['stat'] is None:
            fig.im_set['stat'] = (event.xdata, event.ydata)
            print("Press 's' again at other bound of region.")
            return

        # have 2 successive presses, calculate the statistics
        x1 = min(fig.im_set['stat'][0], event.xdata)
        x2 = max(fig.im_set['stat'][0], event.xdata)

        region = np.where((fig.im_set['xdata'] >= x1) &
                          (fig.im_set['xdata'] <= x2))[0]
        if len(region) == 0:
            print('No valid data points found to calculate statistics.')
            """
            ind1 = np.abs(fig.im_set['xdata'] - x1).argmin()
            if ind1 != len(fig.im_set['xdata']):
                region = np.array([ind1, ind1 + 1])
            else:
                region = np.array([ind1, ind1 - 1])
            """
        else:
            mean = fig.im_set['ydata'][region].mean()
            std = fig.im_set['ydata'][region].std()
            isum = fig.im_set['ydata'][region].sum()
            median = np.median(fig.im_set['ydata'][region])
            ostr = "Between {0:g} and {1:g}:\n".format(x1, x2)
            ostr += "Median={0:g}, mean={1:g}, std={2:g}, sum={3:g}, npix={4:d}"
            print(ostr.format(median, mean, std, isum, len(region)))

        fig.im_set['stat'] = None

    if event.key == ' ':
        if event.xdata is not None and event.inaxes is fig.im_set['ax']:
            pixel = np.abs(event.xdata - fig.im_set['xdata']).argmin()
            txt = 'Pixel = [{0}, {1}] Value = {2}'
            if fig.im_set['lineplot']:
                txt = txt.format(fig.im_set['line'], pixel,
                                 fig.im_set['ydata'][pixel])
            else:
                txt = txt.format(pixel, fig.im_set['line'],
                                 fig.im_set['ydata'][pixel])
            print(txt)

    if event.key == '?':
        print(helpstr)

    fig.im_set['fig'].canvas.draw()


def button_click(label=None, fig=None):
    # only one button for now, but still make sure it's the right one.
    if label == 'Overplot':
        # toggle overplot values
        if not fig.im_set['overplot']:
            fig.im_set['overplot'] = True
        else:
            fig.im_set['overplot'] = False


def implot_open_image(fig, infile=None):
    if infile is None:
        infile = fig.im_set['image'][fig.im_set['index']]
    else:
        infile = file_handler(infile)
        if len(infile) > 0:
            infile = infile[0]

    # open the image
    hdulist = image_open(infile)
    if hdulist is None:
        print(f"Error opening image {infile}")
        return
        # XXX: go to next image

    # transpose to match the IRAF definition of line/column
    fig.im_set['im'] = hdulist[0].data.T

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
                                        (fig.im_set['nlines'] + 1) // 2) - 1)

    # redo the slider
    implot_remove_slider(fig)
    implot_make_slider(fig)

    implot_plot(fig)


def implot_change_line(value=None, fig=None):
    fig.im_set['line'] = int(value)
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
                    valfmt='%d')
    partial = functools.partial(implot_change_line, fig=fig)
    fig.im_set['sid'] = slider.on_changed(partial)
    slider.vline.set_visible(False)
    slider.label.set_x(0.5)
    slider.label.set_y(1.4)
    slider.label.set_ha('center')
    slider.valtext.set_x(0.5)
    slider.valtext.set_y(-0.5)
    slider.valtext.set_ha('center')

    fig.im_set['_slider'] = slider


def implot(image, *, line=None, wcs='logical'):
    """

    Parameters
    ----------
    image
    line : starting line to plot. Defaults to the middle.
    wcs

    Returns
    -------

    """
    # XXX: where does this go?
    # Disable default Matplotlib shortcut keys:
    keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
    for key in keymaps:
        plt.rcParams[key] = ''

    images = file_handler(image)

    # we couldn't find any images to plot
    if len(images) == 0:
        return

    plt.ion()
    fig, ax = plt.subplots()

    # initialize parameters and fix them to the image so they persist
    # with the image
    fig.im_set = copy.deepcopy(im_set_init)

    fig.im_set['fig'] = fig
    fig.im_set['ax'] = ax

    fig.im_set['image'] = images
    fig.im_set['line'] = line
    fig.im_set['wcs'] = wcs

    # open the image and plot it
    implot_open_image(fig)

    # set up the keyboard commands
    partial = functools.partial(implot_keypress, fig=fig)
    fig.im_set['cid'] = fig.canvas.mpl_connect('key_press_event', partial)

    fig.subplots_adjust(right=0.8)

    helptxt = fig.text(0.02, 0.02, 'Press ? for help')
    fig.im_set['helptxt'] = helptxt
    # where user input will appear
    inptxt = fig.text(0.02, 0.02, '')
    fig.im_set['iax'] = inptxt

    # set up the slider allowing the user to scan through lines/columns
    nax = plt.axes([0.85, 0.65, 0.2, 0.25], zorder=-1)

    # label and whether to start with True or False
    check = CheckButtons(nax, ('Overplot',), (fig.im_set['overplot'],))
    nax.set_axis_off()

    partial = functools.partial(button_click, fig=fig)
    check.on_clicked(partial)

    fig.im_set['_check'] = check

    plt.show(block=False)

    return
