from iraf.utils import file_handler
from .combine import Instrument, ccdtypes, ccdsubset, get_header_value
from .combine import file_new_copy, type_max
import numpy as np
import os
from iraf.sys import image_open, image_close
import tempfile

__all__ = ['ccdproc']


class CCD(object):
    """
    The CCD structure:  This structure is used to communicate processing
    parameters between the package procedures.  It contains pointers to
    data, calibration image IMIO pointers, scaling parameters, and the
    correction flags.  The corrections flags indicate which processing
    operations are to be performed.  The subsection parameters do not
    include a step size.  A step size is assumed.  If arbitrary subsampling
    is desired this would be the next generalization.

    Taken from ccdred.h

    """
    def __init__(self):
        # CCD data coordinates
        self.c1 = 0  # CCD starting column
        self.c2 = 0  # CCD starting column
        self.l1 = 0  # CCD starting line
        self.l2 = 0  # CCD ending line

        # Input data
        self.inim = None  # Input image pointer
        self.inc1 = 0  # Input data starting column
        self.inc2 = 0  # Input data ending column
        self.inl1 = 0  # Input data starting line
        self.inl2 = 0  # Input data ending line

        # Output data
        self.outim = None  # Output image pointer
        self.outc1 = 0  # Output data starting column
        self.outc2 = 0  # Output data ending column
        self.outl1 = 0  # Output data starting line
        self.outl2 = 0  # Output data ending line

        # Mask data
        self.maskim = None  # Mask image pointer
        self.maskc1 = 0  # Mask data starting column
        self.maskc2 = 0  # Mask data ending column
        self.maskl1 = 0  # Mask data starting line
        self.maskl2 = 0  # Mask data ending line
        self.maskpm = None  # Mask pointer
        self.maskfp = None  # Mask fixpix data

        # Zero level data
        self.zeroim = None  # Zero level image pointer
        self.zeroc1 = 0  # Zero level data starting column
        self.zeroc2 = 0  # Zero level data ending column
        self.zerol1 = 0  # Zero level data starting line
        self.zerol2 = 0  # Zero level data ending line

        # Dark count data
        self.darkim = None  # Dark count image pointer
        self.darkc1 = 0
        self.darkc2 = 0
        self.darkl1 = 0
        self.darkl2 = 0

        # Flat field data
        self.flatim = None
        self.flatc1 = 0
        self.flatc2 = 0
        self.flatl1 = 0
        self.flatl2 = 0

        # Illumination data
        self.illumim = None
        self.illumc1 = 0
        self.illumc2 = 0
        self.illuml1 = 0
        self.illuml2 = 0

        # Fringe data
        self.fringeim = None
        self.fringec1 = 0
        self.fringec2 = 0
        self.fringel1 = 0
        self.fringel2 = 0

        # Trim section
        self.trimc1 = 0
        self.trimc2 = 0
        self.triml1 = 0
        self.triml2 = 0

        # Bias section
        self.biasc1 = 0
        self.biasc2 = 0
        self.biasl1 = 0
        self.biasl2 = 0

        self.readaxis = 1  # Read out axis (1=cols, 2=lines)
        self.calctype = np.double  # Calculation data type
        self.overscantype = 0  # Overscan type
        self.overscanvec = None  # Pointer to overscan vector
        self.darkscale = 0  # Dark count scale factor
        self.fringescale = 0  # Fringe scale factor
        self.flatscale = 0  # Flat field scale factor
        self.illumscale = 0  # Illumination scale factor
        self.minreplace = 0  # Minimum replacement value
        self.mean = 0  # Mean of output image
        self.cor = 0  # Overall correction flag
        # Individual correction flags
        self.cors = {'fixpix': False, 'trim': False, 'overscan': False,
                     'zerocor': False, 'darkcor': False, 'flatcor': False,
                     'illumcor': False, 'fringecor': False, 'findmean': False,
                     'minrep': False}


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
            pixeltype="real"):
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

    pixeltype is supposed to be a string of 2 words:
    Output and calculation pixel datatypes, e.g. 'real real'. But for ccdproc
    the only calculation dataypes are real or short, and we are just doing
    everything in real, so ignore the second part.

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
        out = image_open(outim, mode='update')

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

        # Set processing parameters applicable to all images.

        # Create the ccd structure.
        ccd = CCD()
        ccd.inim = imin
        ccd.outim = out
        ccd.cor = False
        ccd.cors['fixpix'] = False
        ccd.cors['overscan'] = False
        ccd.cors['trim'] = False
        readaxis = readaxis.strip().lower()
        # XXX: the readaxis that is set in setproc.x is different from the
        # comment in ccdred.h. Make sure it's doing the right thing.
        if readaxis in ['line', 'column']:
            ccd.readaxis = readaxis
        else:
            raise Exception(f'Invalid readaxis parameter {readaxis}')
        ccd.minreplace = minreplace

        # begin set_sections

        # python and IRAF axis order is backwards. first axis in IRAF is last
        # in python
        nc = ccd.inim[0].data.shape[-1]
        nl = ccd.inim[0].data.shape[-2]
        # The default data section is the entire image.
        c1 = 0
        c2 = nc-1
        l1 = 0
        l2 = nl-1
        cs = 1
        ls = 1

        datasec = get_header_value(ccd.inim, instrument, 'datasec')

        # end set_sections

        image_close(imin)
        image_close(out)
        if outtmp:
            # XXX: will eventually need to copy to input parameters
            os.remove(outim)
