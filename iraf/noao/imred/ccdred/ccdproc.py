from iraf.utils import file_handler
from .instruments import Instrument, ccdtypes, ccdsubset, get_header_value
from .instruments import file_new_copy, set_header_value, delete_header_value
from .combine import type_max
from ..ccdred import imagetypes
import numpy as np
import os
from iraf.sys import image_open
import tempfile
import datetime
import time

__all__ = ['ccdproc', 'ccd_section']


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

    Notes about what these values mean and where they come from.
    see IRAF documention for ccdgeometry to explain some of these terms.

    in* : input image 'datasec', defaults to size of the image
    ccd* : input image 'ccdsec', defaults to size of in*, but with starting
            point at 1. so if in* = [5:50], ccd* = [1:46]
    trim* : 'trimsec' parameter which might be 'image' so grabs from the image
            header. defaults to the same as in*.
    bias* : input image 'biassec', defaults to size of the image
    out * : defaults to in *

    in* and ccd* must have the same size in each dimension
    in* and ccd* are both adjusted to be within the trim limits if necessary.
    so if trim = [5:30] and in = [3:27], the new in = [5:27].
    ccd is updated accordingly so if ccd = [10:34], the new ccd = [12:34]

    out* is adjusted by the trim starting points, so that out* becomes
    appropriate way to index into the outarray with size of trim

    so if out*/in* = [10:30] and trim = [5:50], the new out* = [6:26]

    zero*/dark*/flat*/illum*/fringe* :
    grab the image's datasec (defaults to size of the image) starting values
    (c1/l1)
    grab the image's ccdsec starting values (defaults to whatever the
     datasec was)
    set the zero* value to ccd* - ccdsec1 + datasec1

    """
    def __init__(self):
        # NOTE: IRAF treats 2D arrays as [column, line], but in Python/numpy
        # this translates to array[line, column]

        # CCD data coordinates
        self.ccdc1 = 0  # CCD starting column
        self.ccdc2 = 0  # CCD starting column
        self.ccdl1 = 0  # CCD starting line
        self.ccdl2 = 0  # CCD ending line

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

        self.readaxis = 'line'  # Read out axis
        self.calctype = np.double  # Calculation data type
        self.overscantype = 0  # Overscan type
        self.overscanvec = None  # Pointer to overscan vector
        self.darkscale = 0  # Dark count scale factor
        self.fringescale = 0  # Fringe scale factor
        self.flatscale = 0  # Flat field scale factor
        self.illumscale = 0  # Illumination scale factor
        self.minreplace = 0  # Minimum replacement value
        self.mean = 0  # Mean of output image
        self.cor = False  # Overall correction flag
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

        hdulist.close()


def cal_scan(image, nscan, scancor):
    # CAL_SCAN -- Generate name for scan corrected calibration image.
    if not scancor or nscan == 1:
        return image
    root, ext = os.path.splitext(image)
    if not np.isfinite(nscan):
        return f"{root}.1d{ext}"
    else:
        return f"{root}.{nscan:d}{ext}"


def cal_image(hdulist, instrument, ccdtype, nscan, calibs, scancor):
    calimages, nscans, caltypes, subsets = calibs

    useind = None
    if ccdtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
        ct = 0
        for ii, iimage in enumerate(calimages):
            if caltypes[ii] != ccdtype:
                continue
            if ccdtype in ['flat', 'illum', 'fringe']:
                usub = ccdsubset(hdulist, instrument)
                if subsets[ii] != usub:
                    continue
            ct += 1
            if ct == 1:
                useind = ii
            else:
                if nscans[ii] == nscans[useind]:
                    estr = f"Warning: Extra calibration image " \
                           f"{calimages[ii]} ignored"
                    print(estr)
                    # Reset the image type to eliminate further warnings.
                    caltypes[ii] = 'unknown'
                elif nscans[useind] != nscan and (nscans[ii] == nscan or
                                                  nscans[ii] == 1):
                    useind = ii

    # If no calibration image is found then it is an error.
    if useind is None:
        if ccdtype == 'zero':
            estr = "No zero level calibration image found"
        elif ccdtype == 'dark':
            estr = "No dark count calibration image found"
        elif ccdtype == 'flat':
            estr = "No flat field calibration image of subset %s found"
        elif ccdtype in ['flat', 'illum', 'fringe']:
            usub = ccdsubset(hdulist, instrument)
            estr = f"No {ccdtype} calibration image of subset {usub} found"
        else:
            estr = f"Unrecognized ccdtype: {ccdtype}"
        raise Exception(estr)

    useim = calimages[useind]
    if nscan != nscans[useind]:
        if nscan != 1 and nscans[useind] == 1:
            useim = cal_scan(useim, nscan, scancor)
        else:
            estr = f"Cannot find or create calibration with nscan of {nscan}"
            raise Exception(estr)

    # Check that the input image is not the same as the calibration image.
    if os.path.samefile(useim, hdulist.filename()):
        estr = f"Calibration image {useim} is the same as the input image"
        raise Exception(estr)
    return useim


def ccd_section(section, defaults=(None, None, 1, None, None, 1)):
    """
    CCD_SECTION -- Parse a 2D image section into its elements.
    1. The default values must be set by the caller.
    2. A null image section is OK.
    3. The first nonwhitespace character must be '['.
    4. The last interpreted character must be ']'.

    The default values should be given as a list/tuple with order:
    (x_start, x_end, x_step, y_start, y_end, y_step).

    If they are not explicitly given, default start and end values are None,
    default step size is 1.

    E.g. [:, 1:20] will return (None, None, None, 1, 20, None)

    Parameters
    ----------
    section
    defaults

    Returns
    -------
    Tuple of 6 parameters: first dimension start, stop, step and second
    dimension start, stop, step

    """
    if len(defaults) != 6:
        raise Exception('defaults given to ccd_section must have length 6')

    if section is None:
        return defaults

    section = section.strip()
    if len(section) == 0:
        return defaults

    if section[0] != '[' and section[-1] != ']':
        raise Exception(f"Error in 2D image section specification {section}")

    osection = section
    # remove the brackets
    section = section[1:-1]
    dims = section.split(',')
    if len(dims) != 2:
        raise Exception(f"Error in 2D image section specification {osection}")

    ret = []
    defs = [defaults[:3], defaults[3:]]

    for ii, dim in enumerate(dims):
        # get the default values in this dimension
        d1, d2, ds = defs[ii]

        dim = dim.strip()
        split = dim.split(':')
        if len(split) == 1:
            try:
                start = int(split[0])
                end = int(split[0])
                step = ds
            except ValueError:
                start = d1
                end = d2
                step = ds
        elif len(split) == 2 or len(split) == 3:
            step = ds
            try:
                start = int(split[0])
            except ValueError:
                start = d1
            try:
                end = int(split[1])
            except ValueError:
                end = d2
            if len(split) == 3:
                try:
                    step = int(split[2])
                except ValueError:
                    step = ds
        else:
            raise Exception(
                f"Error in 2D image section specification {osection}")
        ret.append(start)
        ret.append(end)
        ret.append(step)
    return ret


def logstring(instring, inim, verbose, logfd):
    """
    Meant to be the same as calling timelog() then ccdlog() in ccdproc.
    These prepend the time and image name to the input string

    Parameters
    ----------
    instring
    inim
    verbose
    logfd

    Returns
    -------

    """
    # the timelog() part
    now = datetime.datetime.now()
    now = now.strftime('%Y-%m-%d %H:%M:%S')
    ostr = f'{now} {instring}'

    # the ccdlog() part
    printstr = f'{inim.filename()}: {ostr}'
    if verbose:
        print(printstr)
    if logfd is not None:
        logfd.write(printstr + '\n')
    return ostr


def process(ccd):
    # the equivalent of proc1r/proc2r.
    # the differences occur where the if readaxis == 'line' bits are.
    # mostly the same function

    # translate into numpy space. that means instead of [c, l] you do [l, c]
    # and also that you're 0-indexed and not inclusive of the end anymore.
    fulloutarr = ccd.inim[0].data * 1
    if ccd.cors['trim']:
        fulloutarr = fulloutarr[ccd.triml1-1:ccd.triml2,
                                ccd.trimc1-1:ccd.trimc2]
    # XXX: need to deal with xt_fpsr
    if ccd.maskfp is not None:
        raise Exception('maskfp not yet implemented')
    # grab the bit we want to correct
    outarr = fulloutarr[ccd.outl1-1:ccd.outl2, ccd.outc1-1:ccd.outc2]

    # make the overscan correction
    if ccd.cors['overscan']:
        overscanc1 = ccd.biasc1 - 1
        noverscan = ccd.biasc2 - overscanc1
        # XXX: get the overscan value/vector. different for line/column
        if ccd.overscantype in ['mean', 'median', 'minmax']:
            overscan = 1
        else:
            overscan = 1
        outarr -= overscan

    # make the zero correction
    if ccd.cors['zerocor']:
        # XXX: for zero and dark and flat, need to check which dimension is 1D
        # and make sure to broadcast in the appropriate direction
        # that's the main difference between line and column readaxis
        if len(ccd.zeroim[0].data.shape) == 1:
            if ccd.readaxis == 'line':
                zerobuf = ccd.zeroim[0].data[ccd.zeroc1:ccd.zeroc2]
            else:
                zerobuf = ccd.zeroim[0].data[ccd.zerol1:ccd.zerol2]
        else:
            zerobuf = ccd.zeroim[0].data[ccd.zerol1-1:ccd.zerol2,
                                         ccd.zeroc1-1:ccd.zeroc2]
        outarr -= zerobuf

    # make the dark correction
    if ccd.cors['darkcor']:
        # XXX: for zero and dark and flat, need to check which dimension is 1D
        # and make sure to broadcast in the appropriate direction
        # that's the main difference between line and column readaxis
        if len(ccd.darkim[0].data.shape) == 1:
            if ccd.readaxis == 'line':
                darkbuf = ccd.darkim[0].data[ccd.darkc1:ccd.darkc2]
            else:
                darkbuf = ccd.darkim[0].data[ccd.darkl1:ccd.darkl2]
            # XXX: Bug in IRAF??
            if ccd.readaxis == 'line':
                darkscale = ccd.flatscale
            else:
                darkscale = ccd.darkscale
        else:
            darkbuf = ccd.darkim[0].data[ccd.darkl1-1:ccd.darkl2,
                                         ccd.darkc1-1:ccd.darkc2]
            darkscale = ccd.darkscale
        outarr -= darkbuf * darkscale

    # make the flat field correction
    if ccd.cors['flatcor']:
        # XXX: for zero and dark and flat, need to check which dimension is 1D
        # and make sure to broadcast in the appropriate direction
        # that's the main difference between line and column readaxis
        if len(ccd.flatim[0].data.shape) == 1:
            if ccd.readaxis == 'line':
                flatbuf = ccd.flatim[0].data[ccd.flatc1:ccd.flatc2]
            else:
                flatbuf = ccd.flatim[0].data[ccd.flatl1:ccd.flatl2]
        else:
            flatbuf = ccd.flatim[0].data[ccd.flatl1-1:ccd.flatl2,
                                         ccd.flatc1-1:ccd.flatc2]
        outarr *= ccd.flatscale / flatbuf

    # do the illumination correction
    if ccd.cors['illumcor']:
        illumbuf = ccd.illumim[0].data[ccd.illuml1-1:ccd.illuml2,
                                       ccd.illumc1-1:ccd.illumc2]
        outarr *= ccd.illumscale / illumbuf

    # do the fringe correction
    if ccd.cors['fringecor']:
        fringebuf = ccd.fringeim[0].data[ccd.fringel1-1:ccd.fringel2,
                                         ccd.fringec1-1:ccd.fringec2]
        outarr -= ccd.fringescale * fringebuf

    if ccd.cors['minrep']:
        outarr[np.where(outarr < ccd.minreplace)] = ccd.minreplace

    if ccd.cors['findmean']:
        ccd.mean = outarr.mean()

    # XXX: is this needed or is fulloutarr updated as a view of outarr?
    fulloutarr[ccd.outl1-1:ccd.outl2, ccd.outc1-1:ccd.outc2] = outarr * 1
    ccd.outim[0].data = fulloutarr * 1


def already_processed(image, instrument, key):
    """
    Meant to be a less confusing version of the IRAF function ccdflag.

    Returns True if the input image header has a value other than the instrument
    default for the key, often indicative that IRAF has already
    processed the image in this way before.

    Returns False if the input image header
    has the default value (or key is not in the header at all), usually
    indicating IRAF has not processed the image for 'key' before.

    Parameters
    ----------
    image
    instrument
    key

    Returns
    -------

    """
    inp = get_header_value(image, instrument, key)
    default = get_header_value(image, instrument, key, default=True)
    return inp != default


def ccdproc(images, *, output=None, ccdtype='object', noproc=False, fixpix=True,
            overscan=True, trim=True, zerocor=True, darkcor=True, flatcor=True,
            illumcor=False, fringecor=False, readcor=False, scancor=False,
            readaxis='line', fixfile=None, biassec=None, trimsec=None,
            zero=None, dark=None, flat=None, illum=None, fringe=None,
            minreplace=1., scantype='shortscan', nscan=1, interactive=False,
            overscan_function='legendre', order=1, sample='*', naverage=1,
            niterate=1, low_reject=3., high_reject=3., grow=0., instrument=None,
            pixeltype="real", logfile=None, verbose=False):
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
    overscan_function
    order
    sample
    naverage
    niterate
    low_reject
    high_reject
    grow
    instrument
    pixeltype
    logfile
    verbose

    pixeltype is supposed to be a string of 2 words:
    Output and calculation pixel datatypes, e.g. 'real real'. But for ccdproc
    the only calculation dataypes are real or short, and we are just doing
    everything in real, so ignore the second part.

    Returns
    -------

    Other Parameters
    ----------------


    """
    inputs = file_handler(images)
    # can't assume the output files exist
    outputs = file_handler(output, exists=False)

    # if the output isn't empty but doesn't match the input length
    if 0 < len(outputs) != len(inputs):
        raise Exception("Input and output lists do not match")

    # was given a string or something else, so set up the instrument object
    if not isinstance(instrument, Instrument):
        instrument = Instrument(instrument)

    # this allows interactive to be 4 valued (yes, no, always yes, always no)
    # to allow for not prompting for every image
    # XXX: set_interactive("", interactive)

    # start of cal_open
    # XXX: allow for None or empty string to process all ccdtypes?
    ccdtype = ccdtype.strip().lower()
    if len(ccdtype) == 0:
        ccdtype = 'none'
    elif ccdtype not in imagetypes:
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

    calibs = (calimages, nscans, caltypes, subsets)

    # end of cal_open

    origccdtype = ccdtype

    # Open the log file.
    logfd = None
    if logfile is not None:
        logfd = open(logfile, 'a')

    # Process each image.
    for imct, image in enumerate(inputs):
        if noproc:
            print(f'{image}:\n')

        imin = image_open(image)
        if imin is None:
            continue

        ccdtype = ccdtypes(imin, instrument)
        if ccdtype != origccdtype:
            imin.close()
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

        # XXX: if noproc = True, IRAF still creates & overwrites an output file?
        # Then deletes it at the end? Is this true? Should I make the temp file
        # if noproc somehow?
        file_new_copy(outim, imin, mode='NEW_COPY', overwrite=True,
                      instrument=instrument)
        out = image_open(outim, mode='update')

        if pixeltype is not None and len(pixeltype) > 0:
            otyp = pixeltype.strip().split()[0]

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
        # in IRAF readaxis line == 1, readaxis column == 2
        if readaxis in ['line', 'column']:
            ccd.readaxis = readaxis
        else:
            raise Exception(f'Invalid readaxis parameter {readaxis}')
        ccd.minreplace = minreplace

        # begin set_sections

        # XXX: need to think about this
        """
        How do I want to handle IRAF/FITS vs numpy axis order? What if we
        want to extend to hdf5 or other formats? C-order vs F-order arrays
        is going to be tough to handle. Same with FITS 1-indexing vs Python
        0-indexing. This ccd_section is designed for 1-indexing and I later
        adjust to 0-indexing, but what about in the future when I want to save
        as hdf5 and its datasec parameter is already 0-indexed?
        
        For now we're just assuming only working with FITS files, but will
        need to be very careful in the future and revisit this and probably 
        rewrite chunks to make less assumptions.
        """
        # python and IRAF axis order is backwards. first axis in IRAF is last
        # in python
        nc = ccd.inim[0].data.shape[-1]
        nl = ccd.inim[0].data.shape[-2]

        # The default data section is the entire image.
        datasec = get_header_value(ccd.inim, instrument, 'datasec')
        c1, c2, cs, l1, l2, ls = ccd_section(datasec, defaults=(1, nc, 1,
                                                                1, nl, 1))
        if c1 < 1 or c2 > nc or cs != 1 or l1 < 1 or l2 > nl or ls != 1:
            raise Exception(f"Error in datasec parameter: {datasec}")

        ccd.inc1 = c1
        ccd.inc2 = c2
        ccd.inl1 = l1
        ccd.inl2 = l2

        # The default trim section is the data section.
        # Defer limit checking until actually used.
        # XXX: set default value to 'image'?
        ts = trimsec
        if trimsec == 'image':
            ts = get_header_value(ccd.inim, instrument, 'trimsec')
        c1, c2, cs, l1, l2, ls = ccd_section(ts, defaults=(c1, c2, 1,
                                                           l1, l2, 1))
        if cs != 1 or ls != 1:
            raise Exception(f"Error in trimsec parameter: {ts}")

        ccd.trimc1 = c1
        ccd.trimc2 = c2
        ccd.triml1 = l1
        ccd.triml2 = l2

        # The default bias section is the whole image.
        # XXX: set default value to 'image'?
        bs = biassec
        if biassec == 'image':
            bs = get_header_value(ccd.inim, instrument, 'biassec')
        c1, c2, cs, l1, l2, ls = ccd_section(bs, defaults=(1, nc, 1, 1, nl, 1))
        if cs != 1 or ls != 1:
            raise Exception(f"Error in biassec parameter: {bs}")

        ccd.biasc1 = c1
        ccd.biasc2 = c2
        ccd.biasl1 = l1
        ccd.biasl2 = l2

        # The default ccd section is the size of the data section.
        c2 = ccd.inc2 - ccd.inc1 + 1
        l2 = ccd.inl2 - ccd.inl1 + 1

        ccs = get_header_value(ccd.inim, instrument, 'ccdsec')
        c1, c2, cs, l1, l2, ls = ccd_section(ccs, defaults=(1, c2, 1, 1, l2, 1))

        if cs != 1 or ls != 1:
            raise Exception(f"Error in ccdsec parameter: {ccs}")

        ccd.ccdc1 = c1
        ccd.ccdc2 = c2
        ccd.ccdl1 = l1
        ccd.ccdl2 = l2

        if (ccd.inc2 - ccd.inc1 != ccd.ccdc2 - ccd.ccdc1 or
                ccd.inl2 - ccd.inl1 != ccd.ccdl2 - ccd.ccdl1):
            raise Exception("Size of DATASEC and CCDSEC do not agree")

        # The default output data section is the input data section.
        ccd.outc1 = ccd.inc1
        ccd.outc2 = ccd.inc2
        ccd.outl1 = ccd.inl1
        ccd.outl2 = ccd.inl2

        # Set the physical WCS to be CCD coordinates.
        # XXX: need to implement a bunch of WCS stuff here.

        # end set_sections

        # begin set_trim

        # also, this all just deals with the ccd object, could pull this out
        # into its own function.
        if trim and not already_processed(ccd.inim, instrument, 'trim'):
            # Check trim section.
            if (ccd.trimc1 < 1 or ccd.trimc2 > nc or ccd.triml1 < 1 or
                    ccd.triml2 > nl):
                estr = f"Error in trim section: image={ccd.inim.filename()}" \
                        f"[{nc},{nl}], trimsec=[{ccd.trimc1}:{ccd.trimc2}," \
                        f"{ccd.triml1}:{ccd.triml2}]"
                raise Exception(estr)
            # If no processing is desired print trim section and return.
            if noproc:
                ostr = f"  [TO BE DONE] Trim section is [{ccd.trimc1}:" \
                       f"{ccd.trimc2},{ccd.triml1}:{ccd.triml2}]."
                print(ostr)
            else:
                # if trim limits are inside inc1, clip inc1 to the trim limits,
                # otherwise leave inc1 alone if trim is wider on either side.
                # then make sure ccdc1 is adjusted the same so they remain the
                # same size
                xt1 = max(0, ccd.trimc1 - ccd.inc1)
                xt2 = min(0, ccd.trimc2 - ccd.inc2)
                yt1 = max(0, ccd.triml1 - ccd.inl1)
                yt2 = min(0, ccd.triml2 - ccd.inl2)

                ccd.ccdc1 += xt1
                ccd.ccdc2 += xt2
                ccd.ccdl1 += yt1
                ccd.ccdl2 += yt2

                ccd.inc1 += xt1
                ccd.inc2 += xt2
                ccd.inl1 += yt1
                ccd.inl2 += yt2

                # output image has the size of the trim limits.
                # even if the input data section does not cover it
                # also this gives outc2 - outc1 as the same size as the
                # inc2 - inc1, but then sets IM_LEN(out) to be the trim lengths
                ccd.outc1 = ccd.inc1 - ccd.trimc1 + 1
                ccd.outc2 = ccd.inc2 - ccd.trimc1 + 1
                ccd.outl1 = ccd.inl1 - ccd.triml1 + 1
                ccd.outl2 = ccd.inl2 - ccd.triml1 + 1

                ccd.cors['trim'] = True
                ccd.cor = True

                ostr = f"Trim data section is [{ccd.trimc1:d}:{ccd.trimc2:d}," \
                       f"{ccd.triml1:d}:{ccd.triml2:d}]"

                logstr = logstring(ostr, ccd.inim, verbose, logfd)
                # XXX: this can't be what is actually put in the header right?
                set_header_value(ccd.outim, instrument, 'trim', logstr)

        # end set_trim

        # begin set_fixpix
        if fixpix and not already_processed(ccd.inim, instrument, 'fixpix'):
            # Get the bad pixel file.  If the name is "image" then get the file
            # name from the image header or symbol table.
            fx = fixfile
            if fixfile == 'image':
                fx = get_header_value(ccd.inim, instrument, 'fixfile')

            # If no processing is desired print message and return.
            if noproc:
                print(f"  [TO BE DONE] Bad pixel file is {fx}")
            else:
                # Map the bad pixel image and return on an error.
                # XXX: this function does some complicated stuff. need to work
                #     it all out
                # im = xt_pmmap (Memc[image], IN_IM(ccd), Memc[image], SZ_FNAME)
                bpm = image_open(fx)

                ccd.maskim = bpm
                # XXX: not entirely sure what imstati (im, IM_PMDES) is doing
                ccd.maskpm = bpm
                # XXX: need to interpret xt_fpinit (MASK_PM(ccd), 2, 3)
                ccd.maskfp = None

                ccd.cor = True
                ccd.cors['fixpix'] = True

                # Log the operation.
                ostr = f"Bad pixel file is {fx}"
                logstr = logstring(ostr, ccd.inim, verbose, logfd)
                # XXX: this can't be what is actually put in the header right?
                set_header_value(ccd.outim, instrument, 'fixpix', logstr)

        # end set_fixpix

        # begin set_overscan
        if overscan and not already_processed(ccd.inim, instrument, 'overscan'):
            # Check bias section.
            if (ccd.biasc1 < 1 or ccd.biasc2 > nc or ccd.biasl1 < 1 or
                    ccd.biasl2 > nl):
                estr = f"Error in bias section: image={ccd.inim.filename()}" \
                       f"[{nc},{nl}], biassec=[{ccd.biasc1}:{ccd.biasc2}," \
                       f"{ccd.biasl1}:{ccd.biasl2}]"
                raise Exception(estr)
            if (ccd.biasc1 == 1 and ccd.biasc2 == nc and ccd.biasl1 == 1 and
                    ccd.biasl2 == nl):
                estr = "Bias section not specified or given as full image"
                raise Exception(estr)
            # If no processing is desired then print overscan strip and return.
            if noproc:
                ostr = f"  [TO BE DONE] Overscan section is [{ccd.biasc1}:" \
                       f"{ccd.biasc2},{ccd.biasl1}:{ccd.biasl2}]."
                print(ostr)
            else:
                ostypes = ["mean", "median", "minmax", "chebyshev",
                           "legendre", "spline3", "spline1"]
                overscan_function = overscan_function.strip().lower()
                if overscan_function not in ostypes:
                    raise Exception(f'Could not recognize overscan function '
                                    f'{overscan_function}')
                # Determine the overscan section parameters. The readout axis
                # determines the type of overscan.  The step sizes are ignored.
                # The limits in the long dimension are replaced by the trim
                # limits.
                if overscan_function in ['mean', 'median', 'minmax']:
                    fitoscan = None
                    if ccd.readaxis == 'column':
                        estr = "Overscan function type not allowed with" \
                               " readaxis of column"
                        raise Exception(estr)
                else:
                    if ccd.readaxis == 'line':
                        first = ccd.biasc1
                        last = ccd.biasc2
                        # it's supposed to be the mean in every line between
                        # c1 and c2
                        amean = ccd.inim[0].data[:, first-1:last].mean(axis=1)
                        # Trim the overscan vector and set the pixel coordinate.
                        trimoscan = amean[ccd.inl1-1:ccd.inl2]
                    else:
                        first = ccd.biasl1
                        last = ccd.biasl2
                        # it's supposed to be the mean in every column between
                        # l1 and l2
                        amean = ccd.inim[0].data[first-1:last, :].mean(axis=0)
                        # Trim the overscan vector and set the pixel coordinate.
                        trimoscan = amean[ccd.inc1-1:ccd.inc2]
                    # XXX: fit_overscan() goes here. needs to be implemented
                    fitoscan = trimoscan * 1
                # Set the CCD structure overscan parameters.
                ccd.cors['overscan'] = True
                ccd.cor = True
                ccd.overscantype = overscan_function
                ccd.overscanvec = fitoscan

                # Log the operation.
                if overscan_function in ['mean', 'median', 'minmax']:
                    ostr = f"Overscan section is [{ccd.biasc1}:{ccd.biasc2}," \
                           f"{ccd.biasl1}:{ccd.biasl2}] with " \
                           f"function={overscan_function}"
                else:
                    ostr = f"Overscan section is [{ccd.biasc1}:{ccd.biasc2}," \
                           f"{ccd.biasl1}:{ccd.biasl2}] with " \
                           f"mean={fitoscan.mean():g}"
                logstr = logstring(ostr, ccd.inim, verbose, logfd)
                # XXX: this can't be what is actually put in the header right?
                set_header_value(ccd.outim, instrument, 'overscan', logstr)

        # end set_overscan

        # Set processing parameters for the standard CCD image types.
        if ccdtype not in ['zero']:
            # begin set_zero
            if zerocor and not already_processed(ccd.inim, instrument,
                                                 'zerocor'):
                # Get the zero level correction image.
                if scancor:
                    znscan = ccdnscan(ccd.inim, instrument, ccdtype, scantype,
                                      nscan, scancor)
                else:
                    znscan = 1

                cal = cal_image(ccd.inim, instrument, 'zero', nscan, calibs,
                                scancor)
                # If no processing is desired print zero correction image
                #  and return.
                if noproc:
                    ot = "  [TO BE DONE] Zero level correction image is {cal}."
                    print(ot)
                else:
                    # Map the image and return on an error.
                    # Process the zero image if necessary.
                    # If nscan > 1 then the zero may not yet exist so create it
                    # from the unscanned zero.
                    # XXX: this bit is more complicated. can call ccdproc
                    # recursively to create this image.
                    zeroim = image_open(cal)

                    # Set the processing parameters in the CCD structure.
                    znc = zeroim[0].data.shape[-1]
                    znl = zeroim[0].data.shape[-2]

                    zdatasec = get_header_value(zeroim, instrument, 'datasec')
                    zc1, zc2, zcs, zl1, zl2, zls = ccd_section(zdatasec, defaults=(1, znc, 1, 1, znl, 1))

                    if (zc1 < 1 or zc2 > znc or zcs != 1 or zl1 < 1 or
                            zl2 > znl or zls != 1):
                        estr = f"Data section error: image={cal}[{znc},{znl}]" \
                               f", datasec=[{zc1}:{zc2},{zl1}:{zl2}]"
                        raise Exception(estr)
                    # save the datasec starting points
                    datac1 = zc1
                    datal1 = zl1

                    # XXX: for zero/dark/flat/illum/fringe, need to make sure
                    # datasec and ccdsec have the same size/shape.
                    zcsec = get_header_value(zeroim, instrument, 'ccdsec')
                    zc1, zc2, zcs, zl1, zl2, zls = ccd_section(zcsec, defaults=(zc1, zc2, zcs, zl1, zl2, zls))

                    if znc == 1:
                        zc1 = ccd.ccdc1
                        zc2 = ccd.ccdc2
                    if znl == 1:
                        zl1 = ccd.ccdl1
                        zl2 = ccd.ccdl2
                    # save the ccdsec starting points
                    ccdc1 = zc1
                    ccdl1 = zl1

                    if (zc1 > ccd.ccdc1 or zc2 < ccd.ccdc2 or
                            zl1 > ccd.ccdl1 or zl2 < ccd.ccdl2):
                        estr = f'CCD section error: input=[{ccd.ccdc1}:' \
                               f'{ccd.ccdc2},{ccd.ccdl1}:{ccd.ccdl2}], ' \
                               f'{cal}=[{zc1}:{zc2},{zl1}:{zl2}]'
                        raise Exception(estr)

                    # make sure this and ccd* are starting at the same physical
                    # pixel on the detector. ensures we're lining up the same
                    # regions on the detector for calibration
                    ccd.zeroim = zeroim
                    ccd.zeroc1 = ccd.ccdc1 - ccdc1 + datac1
                    ccd.zeroc2 = ccd.ccdc2 - ccdc1 + datac1
                    ccd.zerol1 = ccd.ccdl1 - ccdl1 + datal1
                    ccd.zerol2 = ccd.ccdl2 - ccdl1 + datal1

                    ccd.cors['zerocor'] = True
                    ccd.cor = True

                    # Log the operation.
                    ostr = f"Zero level correction image is {cal}"
                    logstr = logstring(ostr, ccd.inim, verbose, logfd)
                    # XXX: can't be what is actually put in the header right?
                    set_header_value(ccd.outim, instrument, 'zerocor', logstr)

            # end set_zero

        if ccdtype not in ['zero', 'dark']:
            # begin set_dark
            if darkcor and not already_processed(ccd.inim, instrument,
                                                 'darkcor'):
                # Get the dark count correction image name.
                if scancor:
                    znscan = ccdnscan(ccd.inim, instrument, ccdtype, scantype,
                                      nscan, scancor)
                else:
                    znscan = 1

                cal = cal_image(ccd.inim, instrument, 'dark', nscan, calibs,
                                scancor)

                # If no processing is desired print dark count image and return.
                if noproc:
                    ot = f"  [TO BE DONE] Dark count correction image is {cal}."
                    print(ot)
                else:
                    # Map the image and return on an error.
                    # Process the dark count image if necessary.
                    # If nscan > 1 then the dark may not yet exist so create it
                    # from the unscanned dark.

                    # XXX: this bit is more complicated. can call ccdproc
                    # recursively to create this image.
                    darkim = image_open(cal)

                    # Set the processing parameters in the CCD structure.
                    dnc = darkim[0].data.shape[-1]
                    dnl = darkim[0].data.shape[-2]

                    ddatasec = get_header_value(darkim, instrument, 'datasec')
                    dc1, dc2, dcs, dl1, dl2, dls = ccd_section(ddatasec, defaults=(1, dnc, 1, 1, dnl, 1))

                    if (dc1 < 1 or dc2 > dnc or dcs != 1 or dl1 < 1 or
                            dl2 > dnl or dls != 1):
                        estr = f"Data section error: image={cal}[{dnc},{dnl}]" \
                               f", datasec=[{dc1}:{dc2},{dl1}:{dl2}]"
                        raise Exception(estr)

                    datac1 = dc1
                    datal1 = dl1

                    dcsec = get_header_value(darkim, instrument, 'ccdsec')
                    dc1, dc2, dcs, dl1, dl2, dls = ccd_section(dcsec, defaults=(dc1, dc2, dcs, dl1, dl2, dls))

                    if dnc == 1:
                        dc1 = ccd.ccdc1
                        dc2 = ccd.ccdc2
                    if dnl == 1:
                        dl1 = ccd.ccdl1
                        dl2 = ccd.ccdl2
                    ccdc1 = dc1
                    ccdl1 = dl1

                    if (dc1 > ccd.ccdc1 or dc2 < ccd.ccdc2 or
                            dl1 > ccd.ccdl1 or dl2 < ccd.ccdl2):
                        estr = f'CCD section error: input=[{ccd.ccdc1}:' \
                               f'{ccd.ccdc2},{ccd.ccdl1}:{ccd.ccdl2}], ' \
                               f'{cal}=[{dc1}:{dc2},{dl1}:{dl2}]'
                        raise Exception(estr)

                    # make sure this and ccd* are starting at the same physical
                    # pixel on the detector. ensures we're lining up the same
                    # regions on the detector for calibration
                    ccd.darkim = darkim
                    ccd.darkc1 = ccd.ccdc1 - ccdc1 + datac1
                    ccd.darkc2 = ccd.ccdc2 - ccdc1 + datac1
                    ccd.darkl1 = ccd.ccdl1 - ccdl1 + datal1
                    ccd.darkl2 = ccd.ccdl2 - ccdl1 + datal1

                    # Get the dark count integration times.
                    # Return an error if not found.
                    dt1 = get_header_value(ccd.inim, instrument, 'darktime')
                    if dt1 is None:
                        dt1 = get_header_value(ccd.inim, instrument, 'exptime')
                    dt2 = get_header_value(darkim, instrument, 'darktime')
                    if dt2 is None:
                        dt2 = get_header_value(darkim, instrument, 'exptime')
                    if dt2 is None or dt2 <= 0.:
                        raise Exception(f"Dark time is zero for {cal}")

                    ccd.darkscale = dt1 / dt2
                    ccd.cors['darkcor'] = True
                    ccd.cor = True

                    # Log the operation.
                    ostr = f"Dark count correction image is {cal} with " \
                           f"scale={ccd.darkscale:g}"
                    logstr = logstring(ostr, ccd.inim, verbose, logfd)
                    # XXX: can't be what is actually put in the header right?
                    set_header_value(ccd.outim, instrument, 'darkcor', logstr)

            # end set_dark

        if ccdtype == 'flat':
            ccd.cors['findmean'] = True
            ccd.cors['minrep'] = True

        if ccdtype not in ['zero', 'dark', 'flat']:
            # begin set_flat
            if flatcor and not already_processed(ccd.inim, instrument,
                                                 'flatcor'):
                # Get the flat field correction image name.
                if scancor:
                    znscan = ccdnscan(ccd.inim, instrument, ccdtype, scantype,
                                      nscan, scancor)
                else:
                    znscan = 1

                cal = cal_image(ccd.inim, instrument, 'flat', nscan, calibs,
                                scancor)

                # If no processing is desired print flat field image and return.
                if noproc:
                    ot = f"  [TO BE DONE] Flat correction image is {cal}."
                    print(ot)
                else:
                    # Map the image and return on an error.
                    # Process the flat field image if necessary.
                    # If nscan > 1 then the flat field may not yet exist
                    # so create it from the unscanned flat field.

                    # XXX: this bit is more complicated. can call ccdproc
                    # recursively to create this image.
                    flatim = image_open(cal)

                    # Set the processing parameters in the CCD structure.
                    fnc = flatim[0].data.shape[-1]
                    fnl = flatim[0].data.shape[-2]

                    fdatasec = get_header_value(flatim, instrument, 'datasec')
                    fc1, fc2, fcs, fl1, fl2, fls = ccd_section(fdatasec, defaults=(1, fnc, 1, 1, fnl, 1))

                    if (fc1 < 1 or fc2 > fnc or fcs != 1 or fl1 < 1 or
                            fl2 > fnl or fls != 1):
                        estr = f"Data section error: image={cal}[{fnc},{fnl}]" \
                               f", datasec=[{fc1}:{fc2},{fl1}:{fl2}]"
                        raise Exception(estr)

                    datac1 = fc1
                    datal1 = fl1

                    fcsec = get_header_value(flatim, instrument, 'ccdsec')
                    fc1, fc2, fcs, fl1, fl2, fls = ccd_section(fcsec, defaults=(fc1, fc2, fcs, fl1, fl2, fls))

                    if fnc == 1:
                        fc1 = ccd.ccdc1
                        fc2 = ccd.ccdc2
                    if fnl == 1:
                        fl1 = ccd.ccdl1
                        fl2 = ccd.ccdl2
                    ccdc1 = fc1
                    ccdl1 = fl1

                    if (fc1 > ccd.ccdc1 or fc2 < ccd.ccdc2 or
                            fl1 > ccd.ccdl1 or fl2 < ccd.ccdl2):
                        estr = f'CCD section error: input=[{ccd.ccdc1}:' \
                               f'{ccd.ccdc2},{ccd.ccdl1}:{ccd.ccdl2}], ' \
                               f'{cal}=[{fc1}:{fc2},{fl1}:{fl2}]'
                        raise Exception(estr)

                    # make sure this and ccd* are starting at the same physical
                    # pixel on the detector. ensures we're lining up the same
                    # regions on the detector for calibration
                    ccd.flatim = flatim
                    ccd.flatc1 = ccd.ccdc1 - ccdc1 + datac1
                    ccd.flatc2 = ccd.ccdc2 - ccdc1 + datac1
                    ccd.flatl1 = ccd.ccdl1 - ccdl1 + datal1
                    ccd.flatl2 = ccd.ccdl2 - ccdl1 + datal1

                    # If no mean value use 1 as the scale factor.
                    fscale = get_header_value(flatim, instrument, 'ccdmean')
                    if fscale is None:
                        fscale = 1.
                    ccd.flatscale = fscale
                    ccd.cors['flatcor'] = True
                    ccd.cor = True

                    # Log the operation.
                    ostr = f"Flat field image is {cal} with " \
                           f"scale={ccd.flatscale:g}"
                    logstr = logstring(ostr, ccd.inim, verbose, logfd)
                    # XXX: can't be what is actually put in the header right?
                    set_header_value(ccd.outim, instrument, 'flatcor', logstr)
            # end set_flat

        if ccdtype not in ['zero', 'dark', 'flat', 'illum']:
            # begin set_illum
            if illumcor and not already_processed(ccd.inim, instrument,
                                                  'illumcor'):
                # Get the illumcor correction image.
                cal = cal_image(ccd.inim, instrument, 'illum', 1, calibs,
                                scancor)

                # If no processing is desired print illumination image
                # name and return.
                if noproc:
                    ot = f"  [TO BE DONE] Illumination correction " \
                         f"image is {cal}."
                    print(ot)
                else:
                    # Return a warning if the illumination flag is missing.
                    illumim = image_open(cal)
                    if not already_processed(illumim, instrument, 'mkillum'):
                        estr = "MKILLUM flag missing from illumination image"
                        raise Exception(estr)

                    # If no mean value for the scale factor compute it.
                    iscale = get_header_value(illumim, instrument, 'ccdmean')
                    ccd.illumscale = iscale

                    itime = get_header_value(illumim, instrument, 'ccdmeant')
                    if itime is None:
                        # XXX: itime = IM_MTIME(im)
                        # supposed to be time of last modification?
                        pass

                    # XXX: itime < IM_MTIME(im)
                    if iscale is None or itime < 0:
                        pass

                    iscale = get_header_value(illumim, instrument, 'ccdmean')
                    if iscale is None:
                        iscale = 1.
                    ccd.illumscale = iscale

                    # Set the processing parameters in the CCD structure.
                    inc = illumim[0].data.shape[-1]
                    inl = illumim[0].data.shape[-2]

                    idatasec = get_header_value(illumim, instrument, 'datasec')
                    ic1, ic2, ics, il1, il2, ils = ccd_section(idatasec, defaults=(1, inc, 1, 1, inl, 1))

                    if (ic1 < 1 or ic2 > inc or ics != 1 or il1 < 1 or
                            il2 > inl or ils != 1):
                        estr = f"Data section error: image={cal}[{inc},{inl}]" \
                               f", datasec=[{ic1}:{ic2},{il1}:{il2}]"
                        raise Exception(estr)

                    datac1 = ic1
                    datal1 = il1

                    icsec = get_header_value(illumim, instrument, 'ccdsec')
                    ic1, ic2, ics, il1, il2, ils = ccd_section(icsec, defaults=(ic1, ic2, ics, il1, il2, ils))

                    ccdc1 = ic1
                    ccdl1 = il1

                    if (ic1 > ccd.ccdc1 or ic2 < ccd.ccdc2 or
                            il1 > ccd.ccdl1 or il2 < ccd.ccdl2):
                        estr = f'CCD section error: input=[{ccd.ccdc1}:' \
                               f'{ccd.ccdc2},{ccd.ccdl1}:{ccd.ccdl2}], ' \
                               f'{cal}=[{ic1}:{ic2},{il1}:{il2}]'
                        raise Exception(estr)

                    # make sure this and ccd* are starting at the same physical
                    # pixel on the detector. ensures we're lining up the same
                    # regions on the detector for calibration
                    ccd.illumim = illumim
                    ccd.illumc1 = ccd.ccdc1 - ccdc1 + datac1
                    ccd.illumc2 = ccd.ccdc2 - ccdc1 + datac1
                    ccd.illuml1 = ccd.ccdl1 - ccdl1 + datal1
                    ccd.illuml2 = ccd.ccdl2 - ccdl1 + datal1

                    ccd.cors['illumcor'] = True
                    ccd.cor = True

                    # Log the operation.
                    ostr = f"Illumination image is {cal} with " \
                           f"scale={ccd.illumscale:g}"
                    logstr = logstring(ostr, ccd.inim, verbose, logfd)
                    # XXX: can't be what is actually put in the header right?
                    set_header_value(ccd.outim, instrument, 'illumcor', logstr)
            # end set_illum

            # begin set_fringe
            if fringecor and not already_processed(ccd.inim, instrument,
                                                   'fringcor'):
                # Get the fringe correction image.
                cal = cal_image(ccd.inim, instrument, 'fringe', 1, calibs,
                                scancor)

                # If no processing is desired print fringe image name
                # and return.
                if noproc:
                    ot = f"  [TO BE DONE] Fringe correction image is {cal}."
                    print(ot)
                else:
                    # Return an error if the fringe flag is missing.
                    fringeim = image_open(cal)
                    if not already_processed(fringeim, instrument, 'mkfringe'):
                        estr = "MKFRINGE flag missing from fringe image"
                        raise Exception(estr)

                    # Set the processing parameters in the CCD structure.
                    fnc = fringeim[0].data.shape[-1]
                    fnl = fringeim[0].data.shape[-2]

                    fdatasec = get_header_value(fringeim, instrument, 'datasec')
                    fc1, fc2, fcs, fl1, fl2, fls = ccd_section(fdatasec, defaults=(1, fnc, 1, 1, fnl, 1))

                    if (fc1 < 1 or fc2 > fnc or fcs != 1 or fl1 < 1 or
                            fl2 > fnl or fls != 1):
                        estr = f"Data section error: image={cal}[{fnc},{fnl}]" \
                               f", datasec=[{fc1}:{fc2},{fl1}:{fl2}]"
                        raise Exception(estr)

                    datac1 = fc1
                    datal1 = fl1

                    fcsec = get_header_value(fringeim, instrument, 'ccdsec')
                    fc1, fc2, fcs, fl1, fl2, fls = ccd_section(fcsec, defaults=(fc1, fc2, fcs, fl1, fl2, fls))

                    ccdc1 = fc1
                    ccdl1 = fl1

                    if (fc1 > ccd.ccdc1 or fc2 < ccd.ccdc2 or
                            fl1 > ccd.ccdl1 or fl2 < ccd.ccdl2):
                        estr = f'CCD section error: input=[{ccd.ccdc1}:' \
                               f'{ccd.ccdc2},{ccd.ccdl1}:{ccd.ccdl2}], ' \
                               f'{cal}=[{fc1}:{fc2},{fl1}:{fl2}]'
                        raise Exception(estr)

                    # Get the scaling factors.
                    # If no fringe scale factor assume 1.
                    fexp1 = get_header_value(ccd.inim, instrument, 'exptime')
                    fexp2 = get_header_value(fringeim, instrument, 'exptime')
                    fscl = get_header_value(fringeim, instrument, 'fringscl')
                    if fscl is None:
                        fscl = 1.

                    ccd.fringescale = fexp1 / fexp2 * fscl

                    # make sure this and ccd* are starting at the same physical
                    # pixel on the detector. ensures we're lining up the same
                    # regions on the detector for calibration
                    ccd.fringeim = fringeim
                    ccd.fringec1 = ccd.ccdc1 - ccdc1 + datac1
                    ccd.fringec2 = ccd.ccdc2 - ccdc1 + datac1
                    ccd.fringel1 = ccd.ccdl1 - ccdl1 + datal1
                    ccd.fringel2 = ccd.ccdl2 - ccdl1 + datal1

                    ccd.cors['fringecor'] = True
                    ccd.cor = True

                    # Log the operation.
                    ostr = f"Fringe image is {cal} with " \
                           f"scale={ccd.fringescale:g}"
                    logstr = logstring(ostr, ccd.inim, verbose, logfd)
                    # XXX: can't be what is actually put in the header right?
                    set_header_value(ccd.outim, instrument, 'fringcor', logstr)
            # end set_fringe

        if ccdtype not in ['zero', 'dark', 'flat', 'illum', 'object', 'comp']:
            ccd.cors['findmean'] = True

        # Do the processing if the COR flag is set.
        if ccd.cor:
            process(ccd)
            # begin set_header

            # Set the data section if it is not the whole image.
            nc = ccd.outim[0].data.shape[-1]
            nl = ccd.outim[0].data.shape[-2]
            if (ccd.outc1 != 1 or ccd.outc2 != nc or ccd.outl1 != 1 or
                    ccd.outl2 != nl):
                hstr = f'[{ccd.outc1}:{ccd.outc2},{ccd.outl1}:{ccd.outl2}]'
                set_header_value(ccd.outim, instrument, 'datasec', hstr)
            else:
                delete_header_value(ccd.outim, instrument, 'datasec')

            # Set the CCD section.
            hstr = f'[{ccd.ccdc1}:{ccd.ccdc2},{ccd.ccdl1}:{ccd.ccdl2}]'
            set_header_value(ccd.outim, instrument, 'ccdsec', hstr)

            # If trimming update the trim and bias section parameters.
            if ccd.cors['trim']:
                delete_header_value(ccd.outim, instrument, 'trimsec')
                ccd.biasc1 = max(1, ccd.biasc1 - ccd.trimc1 + 1)
                ccd.biasc2 = min(nc, ccd.biasc2 - ccd.trimc1 + 1)
                ccd.biasl1 = max(1, ccd.biasl1 - ccd.triml1 + 1)
                ccd.biasl2 = min(nl, ccd.biasl2 - ccd.triml1 + 1)
                if ccd.biasc1 <= ccd.biasc2 and ccd.biasl1 <= ccd.biasl2:
                    hstr = f'[{ccd.biasc1}:{ccd.biasc2},{ccd.biasl1}:{ccd.biasl2}]'
                    set_header_value(ccd.outim, instrument, 'biassec', hstr)
                else:
                    delete_header_value(ccd.outim, instrument, 'biassec')
                # XXX: update some WCS stuff here.

            # Set mean value if desired.
            if ccd.cors['findmean']:
                set_header_value(ccd.outim, instrument, 'ccdmean', ccd.mean)
                set_header_value(ccd.outim, instrument, 'ccdmeant',
                                 int(time.time()))
            # Mark image as processed.
            # timelog()
            now = datetime.datetime.now()
            now = now.strftime('%Y-%m-%d %H:%M:%S')
            ostr = f'{now} CCD processing done'
            set_header_value(ccd.outim, instrument, 'ccdproc', ostr)
            # end set_header

        imin.close()
        out.close()
        if outtmp:
            if ccd.cor:
                # Replace the input image by the corrected image.
                os.replace(outim, image)
            else:
                os.remove(outim)

        if ccd.maskim is not None:
            ccd.maskim.close()
        if ccd.zeroim is not None:
            ccd.zeroim.close()
        if ccd.darkim is not None:
            ccd.darkim.close()
        if ccd.flatim is not None:
            ccd.flatim.close()
        if ccd.illumim is not None:
            ccd.illumim.close()
        if ccd.fringeim is not None:
            ccd.fringeim.close()

        # Do special processing on certain image types.
        if ccdtype == 'zero' and readcor:
            if outtmp:
                readim = image
            else:
                readim = outim
            # begin readcor
            rim = image_open(readim)
            if not already_processed(rim, instrument, 'readcor'):
                if noproc:
                    ostr = f"[TO BE DONE] Convert {readim} to readout " \
                           f"correction "
                    print(ostr)
                    rim.close()
                else:
                    # The default data section is the entire image.
                    nc = rim[0].data.shape[-1]
                    nl = rim[0].data.shape[-2]

                    rdatasec = get_header_value(rim, instrument, 'datasec')
                    inc1, inc2, incs, inl1, inl2, inls = ccd_section(rdatasec, defaults=(1, nc, 1, 1, nl, 1))

                    if (inc1 < 1 or inc2 > nc or inl1 < 1 or inl2 > nl or
                            incs != 1 or inl2 != 1):
                        raise Exception('Error in DATASEC parameter')

                    # The default ccd section is the data section.
                    rccdsec = get_header_value(rim, instrument, 'ccdsec')
                    (ccdc1, ccdc2, ccdcs,
                     ccdl1, ccdl2, ccdls) = ccd_section(rccdsec, defaults=(inc1, inc2, incs, inl1, inl2, inls))
                    if ccdcs != 1 or ccdls != 1:
                        raise Exception('Error in CCDSEC parameter')
                    if (inc2 - inc1 != ccdc2 - ccdc1 or
                            inl2 - inl1 != ccdl2 - ccdl1):
                        raise Exception('Size of DATASEC and CCDSEC do not '
                                        'agree')

                    # XXX: make sure the output file deals with the
                    # 'pixeltype' correctly as in set_output though
                    outdata = rim[0].data * 1

                    rtmp = tempfile.NamedTemporaryFile(delete=False)
                    rtmp.close()
                    newout = rtmp.name
                    file_new_copy(newout, rim, instrument=instrument)
                    newhdr = image_open(newout, mode='update')
                    rim.close()
                    # Average across the readout axis.

                    # zero out the parts not in the data section
                    outdata[:, inc1-1] = 0.
                    outdata[:inl1 - 1:, ] = 0.
                    outdata[:, inc2:] = 0.
                    outdata[inl2:, :] = 0.

                    if readaxis == 'line':
                        # can't just use mean because there might be fewer
                        # rows we actually want the mean of and the rest 0s
                        outdata = outdata.sum(axis=0) / (inl2 - inl1 + 1)
                        ostr = f'[{inc1}:{inc2}, 1:1]'
                        set_header_value(newhdr, instrument, 'datasec', ostr)
                        ostr = f'[{ccdc1}:{ccdc2}, :]'
                        set_header_value(newhdr, instrument, 'ccdsec', ostr)
                    else:
                        outdata = outdata.sum(axis=1) / (inc2 - inc1 + 1)
                        ostr = f'[1:1, {inl1}:{inl2}]'
                        set_header_value(newhdr, instrument, 'datasec', ostr)
                        ostr = f'[:, {ccdl1}:{ccdl2}]'
                        set_header_value(newhdr, instrument, 'ccdsec', ostr)

                    newhdr.close()
                    os.replace(newout, readim)

                    newout = image_open(readim, mode='update')
                    # Log the operation.
                    ostr = "Converted to readout format"
                    logstr = logstring(ostr, newout, verbose, logfd)
                    # XXX: can't be what is actually put in the header right?
                    set_header_value(newout, instrument, 'readcor', logstr)

                    newout.close()
            else:
                rim.close()
            # end readcor

        if ccdtype == 'flat':
            if outtmp:
                readim = image
            else:
                readim = outim
            # begin ccdmean

            # Check if this operation has been done.
            meanim = image_open(readim, mode='update')
            cm = get_header_value(meanim, instrument, 'ccdmean')
            domean = True

            if cm is not None:
                cmt = get_header_value(meanim, instrument, 'ccdmeant')
                lastmod = get_header_value(meanim, instrument, 'iraf-tlm')
                # XXX: somehow compare these times? Or just don't bother and
                # do it anyway since it won't take any time at all
                if cmt is None and lastmod is not None and cmt > lastmod:
                    domean = False

            if domean:
                if noproc:
                    ostr = f"  [TO BE DONE] Compute mean of image {readim}"
                    print(ostr)
                else:
                    # Compute and record the mean.
                    mean = meanim[0].data.mean()

                    set_header_value(meanim, instrument, 'ccdmean', mean)
                    set_header_value(meanim, instrument, 'ccdmeant',
                                     int(time.time()))
            meanim.close()
            # end ccdmean

    if logfd is not None:
        logfd.close()
