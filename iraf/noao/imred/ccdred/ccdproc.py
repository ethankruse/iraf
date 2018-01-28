from iraf.utils import file_handler
from .combine import Instrument, ccdtypes, ccdsubset, get_header_value
from .combine import file_new_copy, type_max
import numpy as np
import os
from iraf.sys import image_open, image_close
import tempfile

__all__ = ['ccdproc']


def ccdnscan(hdulist, instrument, ccdtype, scantype, inscan,
             scancor):
    """
    Return the number CCD scan rows.

    If not found in the header return the "nscan" parameter for objects and
    1 for calibration images.

    Parameters
    ----------
    hdulist
    instrument
    ccdtype
    scantype
    inscan
    scancor

    Returns
    -------

    """
    nscan = get_header_value(hdulist, instrument, 'nscanrow')
    if nscan is None:
        if ccdtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
            nscan = 1
        else:
            if scantype == 'shortscan':
                nscan = inscan
            else:
                if scancor:
                    nscan = -np.inf
                else:
                    nscan = 1
    return nscan


def cal_list(list1, listtype, instrument, calimages, nscans, caltypes, subsets,
             scantype, inscan, scancor):
    for image in list1:
        # Open the image.  If an explicit type is given it is an
        # error if the image can't be opened.
        hdulist = image_open(image)
        if hdulist is None:
            if listtype == 'unknown':
                continue
            else:
                raise Exception(f'Error opening {image}')

        # Override image header CCD type if a list type is given.
        if listtype == 'unknown':
            ccdtype = ccdtypes(hdulist, instrument)
        else:
            ccdtype = listtype

        if ccdtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
            if image not in calimages:
                caltypes.append(ccdtype)
                nsc = ccdnscan(hdulist, instrument, ccdtype, scantype, inscan,
                               scancor)
                nscans.append(nsc)
                subsets.append(ccdsubset(hdulist, instrument))
                calimages.append(image)

        image_close(image)


def ccdproc(images, output, *, ccdtype='object', noproc=False, fixpix=True,
            overscan=True, trim=True, zerocor=True, darkcor=True, flatcor=True,
            illumcor=False, fringecor=False, readcor=False, scancor=False,
            readaxis='line', fixfile=None, biassec=None, trimsec=None,
            zero=None, dark=None, flat=None, illum=None, fringe=None,
            minreplace=1., scantype='shortscan', nscan=1, interactive=False,
            function='legendre', order=1, sample='*', naverage=1, niterate=1,
            low_reject=3., high_reject=3., grow=0., instrument=None,
            pixeltype="real real"):
    """

    Parameters
    ----------
    images
    output
    ccdtype
    noproc
    fixpix
    overscan
    trim
    zerocor
    darkcor
    flatcor
    illumcor
    fringecor
    readcor
    scancor
    readaxis
    fixfile
    biassec
    trimsec
    zero
    dark
    flat
    illum
    fringe
    minreplace
    scantype
    nscan
    interactive
    function
    order
    sample
    naverage
    niterate
    low_reject
    high_reject
    grow
    instrument
    pixeltype

    Returns
    -------

    """
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

    # start of cal_open

    # define CCDTYPES "|object|zero|dark|flat|illum|fringe|other|comp|"
    ccdopts = "object|zero|dark|flat|illum|fringe|other|comp".split('|')
    ccdtype = ccdtype.strip().lower()
    if len(ccdtype) == 0:
        ccdtype = 'none'
    elif ccdtype not in ccdopts:
        ccdtype = 'unknown'

    calimages = []
    nscans = []
    caltypes = []
    subsets = []
    scantype = scantype.strip().lower()

    if ccdtype != 'zero' and zerocor:
        list1 = file_handler(zero)
        cal_list(list1, 'zero', instrument, calimages, nscans,
                 caltypes, subsets, scantype, nscan, scancor)
    if ccdtype not in ['zero', 'dark'] and darkcor:
        list1 = file_handler(dark)
        cal_list(list1, 'dark', instrument, calimages, nscans,
                 caltypes, subsets, scantype, nscan, scancor)
    if ccdtype not in ['zero', 'dark', 'flat'] and flatcor:
        list1 = file_handler(flat)
        cal_list(list1, 'flat', instrument, calimages, nscans,
                 caltypes, subsets, scantype, nscan, scancor)
    if ccdtype not in ['zero', 'dark', 'flat', 'illum'] and illumcor:
        list1 = file_handler(illum)
        cal_list(list1, 'illum', instrument, calimages, nscans,
                 caltypes, subsets, scantype, nscan, scancor)
    if ccdtype not in ['zero', 'dark', 'flat', 'fringe'] and fringecor:
        list1 = file_handler(fringe)
        cal_list(list1, 'fringe', instrument, calimages, nscans,
                 caltypes, subsets, scantype, nscan, scancor)
    cal_list(inputs, 'unknown', instrument, calimages, nscans,
             caltypes, subsets, scantype, nscan, scancor)

    # end of cal_open

    origccdtype = ccdtype

    # Process each image.
    for imct, image in enumerate(inputs):
        if noproc:
            print(f'{image}:\n')

        imin = image_open(image)
        if imin is None:
            continue

        ccdtype = ccdtypes(imin, instrument)
        if ccdtype != origccdtype:
            image_close(imin)
            continue

        # Set output image.
        if len(outputs) > 0:
            outtmp = False
            outim = outputs[imct]
        else:
            outtmp = True
            outims = tempfile.NamedTemporaryFile(delete=False)
            outims.close()
            outim = outims.name

        file_new_copy(outim, imin, mode='NEW_COPY', overwrite=True,
                      instrument=instrument)

        if pixeltype is not None and len(pixeltype) > 0:
            otyp = pixeltype.split()[0]

            otypes = "short|ushort|integer|long|real|double".split('|')
            ndtypes = [np.short, np.ushort, np.int_, np.int_,
                       np.single, np.double]
            if otyp in otypes:
                outtype = ndtypes[otypes.index(otyp)]
            else:
                raise Exception(f'Unrecognized pixeltype: {otyp}')
            # make sure the output type will work given the input
            outtype = type_max(imin[0].data.dtype, outtype)
        else:
            outtype = imin[0].data.dtype

        image_close(imin)
        if outtmp:
            os.remove(outim)
