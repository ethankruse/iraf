import copy
import os
import tempfile
import time
import warnings
from datetime import datetime

import numpy as np

from iraf.sys import image_open
from iraf.utils import file_handler
from . import Instrument
from .utils import CCDProcError, CCDProcWarning, ccdsubset, ccdtypes, \
    delete_header_value, file_new_copy, get_header_value, set_header_value, \
    type_max
from ..ccdred import _imagetypes

__all__ = ['ccdproc']


# XXX: redo this documentation
class CCD(object):
    """
    Taken from ccdred.h:
    The CCD structure:  This structure is used to communicate processing
    parameters between the package procedures.  It contains pointers to
    data, calibration image IMIO pointers, scaling parameters, and the
    correction flags.  The corrections flags indicate which processing
    operations are to be performed.  The subsection parameters do not
    include a step size.  A step size is assumed.  If arbitrary subsampling
    is desired this would be the next generalization.


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
        self.outtype = None  # Data type of the output image
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


def ccd_section(section, defaults=(None, None, 1, None, None, 1)):
    """
    Parse a 2D image section string into its elements.

    Convert section information in data headers into useful bounds. For
    example, convert the header 'biassec' = '[3:10,5:50]' into the bounds
    (3, 10, 1, 5, 50, 1). The 1s indicate a step size of 1 in each dimension.

    The first nonwhitespace character must be '['. The last
    character must be ']'. A None or empty-string section is ok and will
    return the defaults.

    The required section format is [x1:x2:xstep, y1:y2:ystep].
    If they are not explicitly given, default start and end values are None,
    default step size is 1. The second : in either dimension can be ignored
    if not giving a step size. Custom defaults can be passed as an argument.

    E.g. [:, 2:] will return (None, None, 1, 2, None, 1)

    NOTE: this is not an exact copy of IRAF's ccdsection. In the original
    version, you could not use pythonic things like '5:'. Either both start
    and stop had to be specified or else the wildcard * was used. If * is still
    widely in use, this functionality can be added back in. The original
    also flipped the step to be negative if start > stop. We don't do that
    because of Python's use of negative indices to indicate the end of a list.

    Parameters
    ----------
    section : str or None
        The input data section string to be interpreted
    defaults : tuple
        Values filled in if not specified (e.g. with [:, 5:10]).
        The default values should be given as a list/tuple with order:
        (x_start, x_end, x_step, y_start, y_end, y_step).

    Returns
    -------
    tuple
        Tuple of 6 parameters: first dimension start, stop, step and second
        dimension start, stop, step.

    Raises
    ------
    ValueError
        If the input section cannot be parsed.

    """
    if len(defaults) != 6:
        raise IndexError('defaults given to ccd_section must have length 6')

    if section is None:
        return defaults

    section = section.strip()
    if len(section) == 0:
        return defaults

    if section[0] != '[' or section[-1] != ']':
        raise ValueError(f"Error in 2D image section specification {section}")

    osection = section
    # remove the brackets
    section = section[1:-1]
    dims = section.split(',')
    if len(dims) != 2:
        raise ValueError(f"Error in 2D image section specification {osection}")

    ret = []
    defs = [defaults[:3], defaults[3:]]

    for ii, dim in enumerate(dims):
        # get the default values in this dimension
        d1, d2, ds = defs[ii]

        dim = dim.strip()
        split = dim.split(':')
        if len(split) == 1:
            step = ds
            if len(split[0].strip()) == 0:
                start = d1
                end = d2
            else:
                start = int(split[0])
                end = int(split[0])
        elif len(split) == 2 or len(split) == 3:
            step = ds
            if len(split[0].strip()) == 0:
                start = d1
            else:
                start = int(split[0])

            if len(split[1].strip()) == 0:
                end = d2
            else:
                end = int(split[1])

            if len(split) == 3:
                if len(split[2].strip()) > 0:
                    step = int(split[2])
        else:
            raise ValueError(
                f"Error in 2D image section specification {osection}")
        ret.append(start)
        ret.append(end)
        ret.append(step)
    return ret


def ccdnscan(hdulist, instrument, ccdtype, scantype, nscan,
             scancor):
    """
    Return the number of CCD scan rows.

    If not found in the header ('nscanrow' or instrument equivalent),
    return the "nscan" parameter for objects or 1 for calibration images
    (calibration defined as zero, dark, flat, illum, or fringe image types).

    If scantype is longscan, do not return "nscan"; instead return None if
    scancor is True, 1 otherwise.

    Parameters
    ----------
    hdulist : IRAF image
    instrument : Instrument
    ccdtype : str
        "|object|zero|dark|flat|illum|fringe|other|comp|none|unknown|"
    scantype : {'shortscan', 'longscan'}
    nscan : int
    scancor : bool

    Returns
    -------
    float or None

    """
    retscan = get_header_value(hdulist, instrument, 'nscanrow')

    if retscan is None:
        if ccdtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
            retscan = 1
        else:
            if scantype == 'shortscan':
                retscan = nscan * 1
            elif scantype == 'longscan':
                if scancor:
                    retscan = None
                else:
                    retscan = 1
            else:
                raise ValueError(f"Unrecognized scantype: '{scantype}'")
    return retscan


def cal_list(inlist, listtype, instrument, calimages, nscans, caltypes, subsets,
             scantype, nscan, scancor):
    """
    Add calibration images to an internal list.

    Get each image's subset, nscan, and  image type. If listtype is given,
    overwrite whatever image type is found in the header, as this image was
    fed to ccdproc by the user as that image type regardless of its header info.

    Parameters
    ----------
    inlist : list[str]
        Images to be added to the lists.
    listtype : str
        One of the 10 valid image types.
    instrument : Instrument
    calimages : list[str]
        Output list of calibration image file names.
    nscans : list[int]
        Output list of the nscan related to the output images.
    caltypes : list[str]
        Output list of the image types corrsponding to the output images.
    subsets : list[str]
        Subset each output image belongs to.
    scantype : {'shortscan', 'longscan'}
        Scan type for the processing
    nscan : int
    scancor : bool

    Returns
    -------
    None
        The input images are added to the four appropriate output lists given
        as arguments to the function.

    """
    for image in inlist:
        # Open the image.  Unlike IRAF, always an error if it can't be found
        # or opened. In IRAF it's an error unless listtype='unknown'. See note
        # in main ccdproc about this difference.
        # XXX: if we remove the cal_list(inputs, 'unknown'...) from the main
        # ccd_proc, then we can have this always raise the error and use
        # the comment above (because then cal_list should never be called with
        # "unknown" listtype in this package)

        # because we want to gracefully skip over bad input images rather than
        # halting everything completely, allow us to pass by images we can't
        # open if the image type is "unknown", i.e. coming from the input image
        # list.
        try:
            with image_open(image) as hdulist:
                # Override image header CCD type if a list type is given.
                if listtype == 'unknown':
                    ccdtype = ccdtypes(hdulist, instrument)
                else:
                    ccdtype = listtype

                # only add calibration images to our calibration lists
                if ccdtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
                    if image not in calimages:
                        caltypes.append(ccdtype)
                        nsc = ccdnscan(hdulist, instrument, ccdtype, scantype,
                                       nscan, scancor)
                        nscans.append(nsc)
                        subsets.append(ccdsubset(hdulist, instrument))
                        calimages.append(image)
        except OSError:
            if listtype == 'unknown':
                continue
            else:
                raise


def cal_scan(image, nscan, scancor):
    """
    Generate a name for a scan corrected calibration image.

    Parameters
    ----------
    image : str
    nscan : int or None
    scancor : bool

    Returns
    -------
    str

    """
    # make sure we even really want to be here
    if not scancor or nscan == 1:
        return image
    root, ext = os.path.splitext(image)
    # NOTE: Bug in IRAF? This is almost certainly supposed to return *.1.ext,
    # but the actual IRAF does indeed return *.1d.ext
    if nscan is None:
        return f"{root}.1d{ext}"
    else:
        return f"{root}.{nscan:d}{ext}"


def cal_image(hdulist, instrument, ccdtype, nscan, calibs, scancor):
    """
    Find the appropriate calibration image in our list of calibration files.
    Compare the input image's subset and requested calibration image type to
    the images in the calibration list and pick out the one that matches.
    Raises an error if it can't find exactly one appropriate calibration
    image to use.

    Parameters
    ----------
    hdulist : IRAF image
    instrument : Instrument
    ccdtype : str
        image type of calibration we're looking to do
    nscan : int or None
        nscan value of the image to be calibrated
    calibs : (List[str], List[int], List[str], List[str])
        The calibration file paths, nscan values, image types, and subsets
    scancor : bool

    Returns
    -------
    str
    """
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
                    # NOTE: original IRAF doesn't error here and just ignores
                    # the additional calibraiton images and uses the first.
                    estr = f"Multiple calibration images of type '{ccdtype}'" \
                           f" found:\n{calimages[useind]}\n{calimages[ii]}"
                    raise CCDProcError(estr)
                elif nscans[useind] != nscan and (nscans[ii] == nscan or
                                                  nscans[ii] == 1):
                    useind = ii
    else:
        raise ValueError(f"Image type '{ccdtype}' is not a calibration type.")

    # If no calibration image is found then it is an error.
    if useind is None:
        if ccdtype == 'zero':
            estr = "No zero level calibration image found"
        elif ccdtype == 'dark':
            estr = "No dark count calibration image found"
        else:
            usub = ccdsubset(hdulist, instrument)
            estr = f"No {ccdtype} calibration image of subset '{usub}' found"
        raise CCDProcError(estr)

    useim = calimages[useind]
    if nscan != nscans[useind]:
        if nscan != 1 and nscans[useind] == 1:
            useim = cal_scan(useim, nscan, scancor)
        else:
            estr = f"Cannot find or create calibration with nscan of {nscan}"
            raise CCDProcError(estr)

    # Check that the input image is not the same as the calibration image.
    try:
        if os.path.samefile(useim, hdulist.filename()):
            estr = f"Calibration image {useim} is the same as the input image"
            raise CCDProcError(estr)
    # if the second file doesn't exist, that can be ok (e.g. if cal_scan says
    # that we have to create the scan corrected file later.)
    except FileNotFoundError:
        pass
    return useim


def logstring(instring, inim, verbose, logfile):
    """
    Meant to be the same as calling timelog() then ccdlog() in IRAF's ccdproc.
    These prepend the image name and time to the input string.

    Parameters
    ----------
    instring : str
    inim : IRAF image
    verbose : bool
    logfile : str or None

    Returns
    -------
    str
    """
    # the timelog() part
    now = datetime.now()
    now = now.strftime('%Y-%m-%d %H:%M:%S')
    ostr = f'{now} {instring}'

    # the ccdlog() part
    printstr = f'{inim.filename()}: {ostr}'
    if verbose:
        print(printstr)
    if logfile is not None:
        with open(logfile, 'a') as logfd:
            logfd.write(printstr + '\n')
    return ostr


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
    image : IRAF image
    instrument : Instrument
    key : str

    Returns
    -------
    bool
    """
    inp = get_header_value(image, instrument, key)
    default = get_header_value(image, instrument, key, default=True)
    return inp != default


def ccdcheck(image, instrument, ccdtype, flags):
    """
    Check if a calibration image needs further processing before being used.

    Parameters
    ----------
    image : IRAF image
    instrument : Instrument
    ccdtype : str
    flags : dict

    Returns
    -------
    bool

    """
    if flags['trim'] and not already_processed(image, instrument, 'trim'):
        return True
    if flags['fixpix'] and not already_processed(image, instrument, 'fixpix'):
        return True
    if (flags['overscan'] and
            not already_processed(image, instrument, 'overscan')):
        return True

    if ccdtype == 'zero':
        if flags['readcor'] and not already_processed(image, instrument,
                                                      'readcor'):
            return True
    elif flags['zerocor'] and not already_processed(image, instrument,
                                                    'zerocor'):
        return True

    if ccdtype not in ['zero', 'dark']:
        if flags['darkcor'] and not already_processed(image, instrument,
                                                      'darkcor'):
            return True

    if ccdtype == 'flat':
        if flags['scancor'] and not already_processed(image, instrument,
                                                      'scancor'):
            return True
        # XXX: if the ccdmean hasn't been done or is out of date.

    if ccdtype == 'illum':
        if flags['flatcor'] and not already_processed(image, instrument,
                                                      'flatcor'):
            return True
        # XXX: needs the ccdmean to be there

    if ccdtype not in ['zero', 'dark', 'flat', 'illum']:
        if flags['flatcor'] and not already_processed(image, instrument,
                                                      'flatcor'):
            return True
        if flags['illumcor'] and not already_processed(image, instrument,
                                                       'illumcor'):
            return True
        if flags['fringecor'] and not already_processed(image, instrument,
                                                        'fringcor'):
            return True

    return False


def process(ccd):
    """
    Do the actual processing that has been set up by ccdproc.

    Parameters
    ----------
    ccd : CCD

    Returns
    -------

    """
    # the equivalent of proc1r/proc2r plus cor1r/cor2r.
    # Read out axis (1=lines, 2=cols)
    # READAXIS(ccd) = clgwrd ("readaxis",Memc[str],SZ_LINE,"|line|columns|")

    # the differences in those occur where the if readaxis == 'line' bits are.
    # mostly the same function

    # translate into numpy coords. that means instead of [c, l] you do [l, c]
    # and also that you're 0-indexed and not inclusive of the end anymore.
    fulloutarr = ccd.inim[0].data * 1
    fulloutarr = fulloutarr.astype(ccd.calctype)

    if ccd.cors['trim']:
        fulloutarr = fulloutarr[ccd.triml1 - 1:ccd.triml2,
                                ccd.trimc1 - 1:ccd.trimc2]

    # XXX: need to deal with xt_fpsr.
    # it uses maskfp to identify bad pixels and then linearly interpolates
    # over them to "fix" them. Comes from fixpix
    if ccd.maskfp is not None:
        raise NotImplementedError('maskfp and fixpix not yet implemented')

    # grab the bit we want to correct
    outarr = fulloutarr[ccd.outl1 - 1:ccd.outl2, ccd.outc1 - 1:ccd.outc2]

    # make the overscan correction
    if ccd.cors['overscan']:
        if ccd.readaxis == 'line':
            if ccd.overscantype in ['mean', 'median', 'minmax']:
                overscanc1 = ccd.biasc1 - 1
                noverscan = ccd.biasc2 - overscanc1
                oscanreg = outarr[:, overscanc1:overscanc1 + noverscan]
                # begin the find_overscanr function
                if ccd.overscantype == 'minmax':
                    # average with min and max removed.
                    # XXX: figure out the axis bit
                    overscan = np.sum(oscanreg, axis=0)
                    overscan -= oscanreg.min(axis=0) - oscanreg.max(axis=0)
                    # XXX: what happens in IRAF when the overscan region has
                    # a size of 0-2?
                    overscan /= oscanreg.shape[1] - 2
                elif ccd.overscantype == 'median':
                    # NOTE: technically IRAF returns the lower of the middle 2
                    # in the even median case. numpy returns mean of middle 2.
                    # XXX: figure out the axis bit
                    overscan = np.median(oscanreg, axis=0)
                else:
                    # XXX: figure out the axis bit
                    overscan = np.mean(oscanreg, axis=0)
            else:
                overscan = ccd.overscanvec

            # XXX: outarr is 2D, overscan is 1D and need to figure out which way
            # to broadcast
            outarr -= overscan
        else:
            # XXX: outarr is 2D, overscan is 1D and need to figure out which way
            # to broadcast
            outarr -= ccd.overscanvec

    # make the zero correction
    if ccd.cors['zerocor']:
        # XXX: for zero and dark and flat, need to check which dimension is 1D
        # and make sure to broadcast in the appropriate direction
        # that's the main difference between line and column readaxis
        if len(ccd.zeroim[0].data.shape) == 1:
            if ccd.readaxis == 'line':
                zerobuf = ccd.zeroim[0].data[ccd.zeroc1 - 1:ccd.zeroc2]
            else:
                zerobuf = ccd.zeroim[0].data[ccd.zerol1 - 1:ccd.zerol2]
        else:
            zerobuf = ccd.zeroim[0].data[ccd.zerol1 - 1:ccd.zerol2,
                                         ccd.zeroc1 - 1:ccd.zeroc2]
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
            darkbuf = ccd.darkim[0].data[ccd.darkl1 - 1:ccd.darkl2,
                                         ccd.darkc1 - 1:ccd.darkc2]
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
            flatbuf = ccd.flatim[0].data[ccd.flatl1 - 1:ccd.flatl2,
                                         ccd.flatc1 - 1:ccd.flatc2]
        outarr *= ccd.flatscale / flatbuf

    # do the illumination correction
    if ccd.cors['illumcor']:
        illumbuf = ccd.illumim[0].data[ccd.illuml1 - 1:ccd.illuml2,
                                       ccd.illumc1 - 1:ccd.illumc2]
        outarr *= ccd.illumscale / illumbuf

    # do the fringe correction
    if ccd.cors['fringecor']:
        fringebuf = ccd.fringeim[0].data[ccd.fringel1 - 1:ccd.fringel2,
                                         ccd.fringec1 - 1:ccd.fringec2]
        outarr -= ccd.fringescale * fringebuf

    if ccd.cors['minrep']:
        # XXX: use np.clip?
        outarr[np.where(outarr < ccd.minreplace)] = ccd.minreplace

    if ccd.cors['findmean']:
        ccd.mean = outarr.mean()
    ccd.outim[0].data = fulloutarr.astype(ccd.outtype)


def ccdproc(images, *, output=None, ccdtype='object', noproc=False, fixpix=True,
            overscan=True, trim=True, zerocor=True, darkcor=True, flatcor=True,
            illumcor=False, fringecor=False, readcor=False, scancor=False,
            readaxis='line', fixfile=None, biassec='image', trimsec='image',
            zero=None, dark=None, flat=None, illum=None, fringe=None,
            minreplace=1., scantype='shortscan', nscan=1, interactive=False,
            overscan_function='legendre', order=1, sample='*', naverage=1,
            niterate=1, low_reject=3., high_reject=3., grow=0., instrument=None,
            pixeltype='real', logfile=None, verbose=False, overwrite=True):
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
    overwrite

    pixeltype is supposed to be a string of 2 words:
    Output and calculation pixel datatypes, e.g. 'real real'. But for ccdproc
    the only calculation dataypes are real or short, and we are just doing
    everything in real, so ignore the second part. Note that we removed the type
    "long" that is present in iraf since on modern systems that means int64
    which isn't allowed by IRAF since it hasn't updated to the newest FITS
    standard.

    Returns
    -------

    Other Parameters
    ----------------


    Differences from IRAF
    ---------------------
    Input/output lists of different sizes. Or if an input file is corrupted or
    has a different data type, the outputs don't match up

    bug in setoutput.x: if pixeltype is short, it will allow input ushort to be
    cast to short, and if pixeltype is ushort it will allow a short image to be
    cast to ushort. This can be dangerous.
    In IRAF, an int that is supposed to be negative gets set to 0 as a ushort
    rather than wrapping around positive, and same for a ushort int that gets
    cast to an int is clipped at 32767 rather than wrapping negative.
    Similarly, IRAF allows integer to be converted unsafely to real, when
    it should only allow upscaling to double safely.

    """
    # XXX: at the end, comment this and make sure all inputs are being used
    # save all the inputs for recursive ccdproc calls
    flags = locals()

    inputs = file_handler(images)
    # can't assume the output files exist
    outputs = file_handler(output, exists=False)

    # if the output isn't empty but doesn't match the input length, we have
    # a problem
    if 0 < len(outputs) != len(inputs):
        raise ValueError("Input and output lists do not match")

    # was given a string or something else, so set up the instrument object
    if not isinstance(instrument, Instrument):
        instrument = Instrument(instrument)

    # start of cal_open
    if ccdtype is None or len(ccdtype.strip().lower()) == 0:
        ccdtype = 'none'
    elif ccdtype.strip().lower() not in _imagetypes:
        ccdtype = 'unknown'
    ccdtype = ccdtype.strip().lower()

    calimages = []
    nscans = []
    caltypes = []
    subsets = []
    scantype = scantype.strip().lower()
    if scantype not in ['shortscan', 'longscan']:
        raise ValueError(f"Unrecognized scantype: '{scantype}'")

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

    # XXX: do we want to explicitly require calibration images be fed
    # via their appropriate argument rather than allowing everything to be
    # sent in as a single input list and relying on ccdproc to be "smart"
    # and sort out calibration files? If we want to be pythonic and
    # explicit, I think we need to remove this line.
    cal_list(inputs, 'unknown', instrument, calimages, nscans,
             caltypes, subsets, scantype, nscan, scancor)

    calibs = (calimages, nscans, caltypes, subsets)
    # end of cal_open

    origccdtype = ccdtype

    # Process each image.
    for imct, image in enumerate(inputs):
        if noproc:
            print(f'{image}:\n')

        # NOTE: Bug in IRAF here (?). If input and output lists are the same
        # size, but an image in the input list doesn't exist, the i+1 input
        # image will then get processed to the i output image name because
        # the output list does not advance to match the step of the input list.
        # For us, we'll give the error that input and output
        # lists don't match, since the input list is created using exists=True.
        # We therefore can't process any images if any one in the list doesn't
        # exist or can't be opened. IRAF will process all but that one, but the
        # output files will not be matched up in names as expected.
        try:
            imin = image_open(image)
        except OSError:
            warnings.warn(f"Could not open image {image}; skipping to next.",
                          category=CCDProcWarning)
            continue

        ccdtype = ccdtypes(imin, instrument)
        if ccdtype != origccdtype:
            imin.close()
            warnings.warn(f'Image {image} has type "{ccdtype}" instead of '
                          f'"{origccdtype}"; skipping to next.',
                          category=CCDProcWarning)
            continue

        # Set output image.
        if not noproc:
            if len(outputs) > 0:
                replace_input = False
                outreal = outputs[imct]
            else:
                outreal = image
                replace_input = True
            outims = tempfile.NamedTemporaryFile(delete=False)
            outims.close()
            outim = outims.name
            file_new_copy(outim, imin, mode='NEW_COPY', overwrite=overwrite,
                          instrument=instrument)
            out = image_open(outim, mode='update')
        else:
            replace_input = False
            outim = None
            out = None
            outreal = None

        if pixeltype is not None and len(pixeltype) > 0:
            otyp = pixeltype.strip().split()[0]

            otypes = "short|ushort|integer|real|double".split('|')
            ndtypes = [np.short, np.ushort, np.int32, np.single, np.double]
            if otyp in otypes:
                outtype = ndtypes[otypes.index(otyp)]
            else:
                raise ValueError(f'Unknown pixeltype: {otyp}')
            # make sure the output type will work given the input
            try:
                outtype = type_max(imin[0].data.dtype, outtype)
            except TypeError:
                # trying to cast a ushort to short or vice versa, which
                # is dangerous and shouldn't be allowed. Or else
                # int to real, which is also dangerous
                outtype = imin[0].data.dtype
                warnings.warn(f"Cannot cast input array of type {outtype} to "
                              f"requested pixeltype '{otyp}'. Output will "
                              f"remain {outtype}.", category=CCDProcWarning)
        else:
            raise ValueError(f'Unknown pixeltype: {pixeltype}')

        # Set processing parameters applicable to all images.
        # Create the ccd structure.
        ccd = CCD()
        ccd.inim = imin
        ccd.outim = out
        ccd.outtype = outtype
        readaxis = readaxis.strip().lower()
        # in IRAF readaxis line == 1, readaxis column == 2
        if readaxis in ['line', 'column']:
            ccd.readaxis = readaxis
        else:
            raise ValueError(f'Invalid readaxis parameter {readaxis}')
        ccd.minreplace = minreplace

        # begin set_sections
        """
        How do I want to handle IRAF/FITS vs numpy axis order? What if we
        want to extend to hdf5 or other formats? C-order vs F-order arrays
        is going to be tough to handle. Same with FITS 1-indexing vs Python
        0-indexing. This ccd_section is designed for 1-indexing and I later
        adjust to 0-indexing, but what about in the future when I want to save
        as hdf5 and its datasec parameter is already 0-indexed?
        
        For now we're just assuming we're only working with FITS files, but will
        need to be very careful in the future and revisit this and probably 
        rewrite chunks to make less assumptions.
        """
        # python and IRAF axis order is backwards. first axis in IRAF is last
        # in python
        nc = ccd.inim[0].data.shape[1]
        nl = ccd.inim[0].data.shape[0]

        # The default data section is the entire image.
        datasec = get_header_value(ccd.inim, instrument, 'datasec')
        c1, c2, cs, l1, l2, ls = ccd_section(datasec, defaults=(1, nc, 1,
                                                                1, nl, 1))
        if c1 < 1 or c2 > nc or cs != 1 or l1 < 1 or l2 > nl or ls != 1:
            raise CCDProcError(f"Error in datasec parameter: {datasec}")

        ccd.inc1 = c1
        ccd.inc2 = c2
        ccd.inl1 = l1
        ccd.inl2 = l2

        # The default trim section is the data section.
        # Defer limit checking until actually used.
        ts = trimsec
        if trimsec == 'image':
            ts = get_header_value(ccd.inim, instrument, 'trimsec')
        c1, c2, cs, l1, l2, ls = ccd_section(ts, defaults=(c1, c2, 1,
                                                           l1, l2, 1))
        if cs != 1 or ls != 1:
            raise CCDProcError(f"Error in trimsec parameter: {ts}")

        ccd.trimc1 = c1
        ccd.trimc2 = c2
        ccd.triml1 = l1
        ccd.triml2 = l2

        # The default bias section is the whole image.
        bs = biassec
        if biassec == 'image':
            bs = get_header_value(ccd.inim, instrument, 'biassec')
        c1, c2, cs, l1, l2, ls = ccd_section(bs, defaults=(1, nc, 1, 1, nl, 1))
        if cs != 1 or ls != 1:
            raise CCDProcError(f"Error in biassec parameter: {bs}")

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
            raise CCDProcError(f"Error in ccdsec parameter: {ccs}")

        ccd.ccdc1 = c1
        ccd.ccdc2 = c2
        ccd.ccdl1 = l1
        ccd.ccdl2 = l2

        if (ccd.inc2 - ccd.inc1 != ccd.ccdc2 - ccd.ccdc1 or
                ccd.inl2 - ccd.inl1 != ccd.ccdl2 - ccd.ccdl1):
            raise CCDProcError("Size of DATASEC and CCDSEC do not agree")

        # The default output data section is the input data section.
        ccd.outc1 = ccd.inc1
        ccd.outc2 = ccd.inc2
        ccd.outl1 = ccd.inl1
        ccd.outl2 = ccd.inl2

        # Set the physical WCS to be CCD coordinates.
        # XXX: need to implement a bunch of WCS stuff here.

        # end set_sections

        # begin set_trim
        if trim and not already_processed(ccd.inim, instrument, 'trim'):
            # Check trim section.
            if (ccd.trimc1 < 1 or ccd.trimc2 > nc or ccd.triml1 < 1 or
                    ccd.triml2 > nl):
                estr = f"Error in trim section: image={ccd.inim.filename()}" \
                       f"[{nc},{nl}], trimsec=[{ccd.trimc1}:{ccd.trimc2}," \
                       f"{ccd.triml1}:{ccd.triml2}]"
                raise CCDProcError(estr)
            # If no processing is desired print trim section and return.
            if noproc:
                ostr = f"  [TO BE DONE] Trim section is [{ccd.trimc1}:" \
                       f"{ccd.trimc2},{ccd.triml1}:{ccd.triml2}]."
                print(ostr)
            else:
                # if trim limits are inside inc1 (datasec), clip inc1 to the
                # trim limits, otherwise leave inc1 alone if trim is wider on
                # either side. Then make sure ccdc1 is adjusted the same so
                # they remain the same size
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

                logstr = logstring(ostr, ccd.inim, verbose, logfile)
                set_header_value(ccd.outim, instrument, 'trim', logstr)

        # end set_trim

        # begin set_fixpix
        if fixpix and not already_processed(ccd.inim, instrument, 'fixpix'):
            """
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
                logstr = logstring(ostr, ccd.inim, verbose, logfile)
                set_header_value(ccd.outim, instrument, 'fixpix', logstr)
            """
            raise NotImplementedError("fixpix not yet implemented.")
        # end set_fixpix

        # begin set_overscan
        if overscan and not already_processed(ccd.inim, instrument, 'overscan'):
            # Check bias section.
            if (ccd.biasc1 < 1 or ccd.biasc2 > nc or ccd.biasl1 < 1 or
                    ccd.biasl2 > nl):
                estr = f"Error in bias section: image={ccd.inim.filename()}" \
                       f"[{nc},{nl}], biassec=[{ccd.biasc1}:{ccd.biasc2}," \
                       f"{ccd.biasl1}:{ccd.biasl2}]"
                raise CCDProcError(estr)
            if (ccd.biasc1 == 1 and ccd.biasc2 == nc and ccd.biasl1 == 1 and
                    ccd.biasl2 == nl):
                estr = "Bias section not specified or given as full image"
                raise CCDProcError(estr)

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
                    raise CCDProcError(f'Could not recognize overscan function '
                                       f'{overscan_function}')
                # Determine the overscan section parameters. The readout axis
                # determines the type of overscan.  The step sizes are ignored.
                # The limits in the long dimension are replaced by the trim
                # limits.
                if overscan_function in ['mean', 'median', 'minmax']:
                    fitoscan = None
                    if ccd.readaxis == 'column':
                        estr = f"Overscan function type {overscan_function}" \
                               " not allowed with readaxis of column"
                        raise ValueError(estr)
                else:
                    if ccd.readaxis == 'line':
                        first = ccd.biasc1
                        last = ccd.biasc2
                        # it's supposed to be the mean in every line between
                        # c1 and c2
                        amean = ccd.inim[0].data[:, first - 1:last].mean(axis=1)
                        # Trim the overscan vector and set the pixel coordinate.
                        veclen = ccd.ccdl2 - ccd.ccdl1 + 1
                        trimoscan = amean[ccd.inl1 - 1:ccd.inl1 - 1 + veclen]
                    else:
                        first = ccd.biasl1
                        last = ccd.biasl2
                        # it's supposed to be the mean in every column between
                        # l1 and l2
                        amean = ccd.inim[0].data[first - 1:last, :].mean(axis=0)
                        # Trim the overscan vector and set the pixel coordinate.
                        veclen = ccd.ccdc2 - ccd.ccdc1 + 1
                        trimoscan = amean[ccd.inc1 - 1:ccd.inc1 - 1 + veclen]
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
                logstr = logstring(ostr, ccd.inim, verbose, logfile)
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

                cal = cal_image(ccd.inim, instrument, 'zero', znscan, calibs,
                                scancor)
                # If no processing is desired print zero correction image
                #  and return.
                if noproc:
                    ot = f"  [TO BE DONE] Zero level correction image is {cal}."
                    print(ot)
                else:
                    # Map the image and return on an error.
                    # Process the zero image if necessary.
                    # If nscan > 1 then the zero may not yet exist so create it
                    # from the unscanned zero.
                    try:
                        zeroim = image_open(cal)
                    except OSError:
                        scancal = cal_image(ccd.inim, instrument, 'zero', 1,
                                            calibs, scancor)
                        scanim = image_open(scancal)
                        if ccdcheck(scanim, instrument, 'zero', flags):
                            scanim.close()
                            calflags = copy.deepcopy(flags)
                            calflags['images'] = scancal
                            calflags['output'] = None
                            calflags['ccdtype'] = 'zero'
                            ccdproc(**calflags)
                        # XXX: calls scancor() Need to write that.
                        # converts from scancal to cal.
                        zeroim = image_open(cal)

                    if ccdcheck(zeroim, instrument, 'zero', flags):
                        zeroim.close()
                        calflags = copy.deepcopy(flags)
                        calflags['images'] = cal
                        calflags['output'] = None
                        calflags['ccdtype'] = 'zero'
                        ccdproc(**calflags)
                        zeroim = image_open(cal)

                    # Set the processing parameters in the CCD structure.
                    znc = zeroim[0].data.shape[1]
                    znl = zeroim[0].data.shape[0]

                    zdatasec = get_header_value(zeroim, instrument, 'datasec')
                    rr = ccd_section(zdatasec, defaults=(1, znc, 1, 1, znl, 1))
                    zc1a, zc2a, zcsa, zl1a, zl2a, zlsa = rr

                    if (zc1a < 1 or zc2a > znc or zcsa != 1 or zl1a < 1 or
                            zl2a > znl or zlsa != 1):
                        estr = f"Data section error: image={cal}[{znc},{znl}]" \
                               f", datasec=[{zc1a}:{zc2a},{zl1a}:{zl2a}]"
                        raise CCDProcError(estr)
                    # save the datasec starting points
                    datac1 = zc1a
                    datal1 = zl1a

                    zcsec = get_header_value(zeroim, instrument, 'ccdsec')
                    rr = ccd_section(zcsec, defaults=(zc1a, zc2a, zcsa,
                                                      zl1a, zl2a, zlsa))
                    zc1, zc2, zcs, zl1, zl2, zls = rr

                    if zc2a - zc1a != zc2 - zc1 or zl2a - zl1a != zl2 - zl1:
                        raise CCDProcError(f"Size of DATASEC and CCDSEC do not"
                                           f" agree in zero image {cal}")

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
                        raise CCDProcError(estr)

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
                    logstr = logstring(ostr, ccd.inim, verbose, logfile)
                    set_header_value(ccd.outim, instrument, 'zerocor', logstr)
            # end set_zero

        if ccdtype not in ['zero', 'dark']:
            # begin set_dark
            if darkcor and not already_processed(ccd.inim, instrument,
                                                 'darkcor'):
                # Get the dark count correction image name.
                if scancor:
                    dnscan = ccdnscan(ccd.inim, instrument, ccdtype, scantype,
                                      nscan, scancor)
                else:
                    dnscan = 1

                cal = cal_image(ccd.inim, instrument, 'dark', dnscan, calibs,
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
                    try:
                        darkim = image_open(cal)
                    except OSError:
                        scancal = cal_image(ccd.inim, instrument, 'dark', 1,
                                            calibs, scancor)
                        scanim = image_open(scancal)
                        if ccdcheck(scanim, instrument, 'dark', flags):
                            scanim.close()
                            calflags = copy.deepcopy(flags)
                            calflags['images'] = scancal
                            calflags['output'] = None
                            calflags['ccdtype'] = 'dark'
                            ccdproc(**calflags)
                        # XXX: calls scancor() Need to write that.
                        # converts from scancal to cal.
                        darkim = image_open(cal)

                    if ccdcheck(darkim, instrument, 'dark', flags):
                        darkim.close()
                        calflags = copy.deepcopy(flags)
                        calflags['images'] = cal
                        calflags['output'] = None
                        calflags['ccdtype'] = 'dark'
                        ccdproc(**calflags)
                        darkim = image_open(cal)

                    # Set the processing parameters in the CCD structure.
                    dnc = darkim[0].data.shape[1]
                    dnl = darkim[0].data.shape[0]

                    ddatasec = get_header_value(darkim, instrument, 'datasec')
                    rr = ccd_section(ddatasec, defaults=(1, dnc, 1, 1, dnl, 1))
                    dc1, dc2, dcs, dl1, dl2, dls = rr

                    if (dc1 < 1 or dc2 > dnc or dcs != 1 or dl1 < 1 or
                            dl2 > dnl or dls != 1):
                        estr = f"Data section error: image={cal}[{dnc},{dnl}]" \
                               f", datasec=[{dc1}:{dc2},{dl1}:{dl2}]"
                        raise CCDProcError(estr)

                    datac1 = dc1
                    datal1 = dl1

                    dcsec = get_header_value(darkim, instrument, 'ccdsec')
                    rr = ccd_section(dcsec, defaults=(dc1, dc2, dcs,
                                                      dl1, dl2, dls))
                    dc1, dc2, dcs, dl1, dl2, dls = rr
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
                        raise CCDProcError(estr)

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
                        raise ValueError(f"Dark time is zero for {cal}")

                    ccd.darkscale = dt1 / dt2
                    ccd.cors['darkcor'] = True
                    ccd.cor = True

                    # Log the operation.
                    ostr = f"Dark count correction image is {cal} with " \
                           f"scale={ccd.darkscale:g}"
                    logstr = logstring(ostr, ccd.inim, verbose, logfile)
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
                    flnscan = ccdnscan(ccd.inim, instrument, ccdtype, scantype,
                                       nscan, scancor)
                else:
                    flnscan = 1

                cal = cal_image(ccd.inim, instrument, 'flat', flnscan, calibs,
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
                    try:
                        flatim = image_open(cal)
                    except OSError:
                        scancal = cal_image(ccd.inim, instrument, 'flat', 1,
                                            calibs, scancor)
                        scanim = image_open(scancal)
                        if ccdcheck(scanim, instrument, 'flat', flags):
                            scanim.close()
                            calflags = copy.deepcopy(flags)
                            calflags['images'] = scancal
                            calflags['output'] = None
                            calflags['ccdtype'] = 'flat'
                            ccdproc(**calflags)
                        # XXX: calls scancor() Need to write that.
                        # converts from scancal to cal.
                        flatim = image_open(cal)

                    if ccdcheck(flatim, instrument, 'flat', flags):
                        flatim.close()
                        calflags = copy.deepcopy(flags)
                        calflags['images'] = cal
                        calflags['output'] = None
                        calflags['ccdtype'] = 'flat'
                        ccdproc(**calflags)
                        flatim = image_open(cal)

                    # Set the processing parameters in the CCD structure.
                    fnc = flatim[0].data.shape[1]
                    fnl = flatim[0].data.shape[0]

                    fdatasec = get_header_value(flatim, instrument, 'datasec')
                    rr = ccd_section(fdatasec, defaults=(1, fnc, 1, 1, fnl, 1))
                    fc1, fc2, fcs, fl1, fl2, fls = rr

                    if (fc1 < 1 or fc2 > fnc or fcs != 1 or fl1 < 1 or
                            fl2 > fnl or fls != 1):
                        estr = f"Data section error: image={cal}[{fnc},{fnl}]" \
                               f", datasec=[{fc1}:{fc2},{fl1}:{fl2}]"
                        raise CCDProcError(estr)

                    datac1 = fc1
                    datal1 = fl1

                    fcsec = get_header_value(flatim, instrument, 'ccdsec')
                    rr = ccd_section(fcsec, defaults=(fc1, fc2, fcs,
                                                      fl1, fl2, fls))
                    fc1, fc2, fcs, fl1, fl2, fls = rr

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
                        raise CCDProcError(estr)

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
                    logstr = logstring(ostr, ccd.inim, verbose, logfile)
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
                        raise CCDProcError(estr)

                    # If no mean value for the scale factor compute it.
                    iscale = get_header_value(illumim, instrument, 'ccdmean')
                    ccd.illumscale = iscale

                    modtime = os.path.getmtime(cal)
                    epoch = datetime(1980, 1, 1, 0, 0).timestamp()
                    modtime -= epoch
                    modtime = int(modtime)

                    # NOTE: ccdmean opens a file under READ_WRITE mode.
                    # This appears to mess with the image modification time.
                    # It sets it as an integer value, but on my computer seems
                    # to set the modification time 3588 seconds earlier than
                    # the system time. ??? I have no idea where this comes from.
                    # has something to do with IM_MTIME/IM_SVMTIME?
                    # I saw somewhere it was adding 4 seconds to times, so maybe
                    # it's 4*3 seconds ahead then an hour behind for some
                    # reason?
                    # The creation date and IRAF-TLM in the header are also
                    # an hour behind what I'm writing in the header, and what
                    # the current UTC time is, so something is wrong there.
                    # the IRAF-TLM is 4 seconds ahead of creation though, so
                    # that part is true. It appears to be a known IRAF issue
                    # dealing with DST.
                    itime = get_header_value(illumim, instrument, 'ccdmeant')
                    if itime is None:
                        itime = modtime

                    # XXX: call ccdmean
                    if iscale is None or itime < modtime:
                        pass

                    iscale = get_header_value(illumim, instrument, 'ccdmean')
                    if iscale is None:
                        iscale = 1.

                    # XXX: this is a hack for now. Need to figure out
                    # what part of the image ccdmean is called on.
                    # when set with ccdproc, only uses the datasec or whatever
                    # part is processed, but when ccdmean() sets it, it uses
                    # the whole image?
                    iscale = illumim[0].data.mean()

                    ccd.illumscale = iscale

                    # Set the processing parameters in the CCD structure.
                    inc = illumim[0].data.shape[1]
                    inl = illumim[0].data.shape[0]

                    idatasec = get_header_value(illumim, instrument, 'datasec')
                    rr = ccd_section(idatasec, defaults=(1, inc, 1, 1, inl, 1))
                    ic1, ic2, ics, il1, il2, ils = rr

                    if (ic1 < 1 or ic2 > inc or ics != 1 or il1 < 1 or
                            il2 > inl or ils != 1):
                        estr = f"Data section error: image={cal}[{inc},{inl}]" \
                               f", datasec=[{ic1}:{ic2},{il1}:{il2}]"
                        raise CCDProcError(estr)

                    datac1 = ic1
                    datal1 = il1

                    icsec = get_header_value(illumim, instrument, 'ccdsec')
                    rr = ccd_section(icsec, defaults=(ic1, ic2, ics,
                                                      il1, il2, ils))
                    ic1, ic2, ics, il1, il2, ils = rr

                    ccdc1 = ic1
                    ccdl1 = il1

                    if (ic1 > ccd.ccdc1 or ic2 < ccd.ccdc2 or
                            il1 > ccd.ccdl1 or il2 < ccd.ccdl2):
                        estr = f'CCD section error: input=[{ccd.ccdc1}:' \
                               f'{ccd.ccdc2},{ccd.ccdl1}:{ccd.ccdl2}], ' \
                               f'{cal}=[{ic1}:{ic2},{il1}:{il2}]'
                        raise CCDProcError(estr)

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
                    logstr = logstring(ostr, ccd.inim, verbose, logfile)
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
                        raise CCDProcError(estr)

                    # Set the processing parameters in the CCD structure.
                    fnc = fringeim[0].data.shape[1]
                    fnl = fringeim[0].data.shape[0]

                    fdatasec = get_header_value(fringeim, instrument, 'datasec')
                    rr = ccd_section(fdatasec, defaults=(1, fnc, 1, 1, fnl, 1))
                    fc1, fc2, fcs, fl1, fl2, fls = rr

                    if (fc1 < 1 or fc2 > fnc or fcs != 1 or fl1 < 1 or
                            fl2 > fnl or fls != 1):
                        estr = f"Data section error: image={cal}[{fnc},{fnl}]" \
                               f", datasec=[{fc1}:{fc2},{fl1}:{fl2}]"
                        raise CCDProcError(estr)

                    datac1 = fc1
                    datal1 = fl1

                    fcsec = get_header_value(fringeim, instrument, 'ccdsec')
                    rr = ccd_section(fcsec, defaults=(fc1, fc2, fcs,
                                                      fl1, fl2, fls))
                    fc1, fc2, fcs, fl1, fl2, fls = rr

                    ccdc1 = fc1
                    ccdl1 = fl1

                    if (fc1 > ccd.ccdc1 or fc2 < ccd.ccdc2 or
                            fl1 > ccd.ccdl1 or fl2 < ccd.ccdl2):
                        estr = f'CCD section error: input=[{ccd.ccdc1}:' \
                               f'{ccd.ccdc2},{ccd.ccdl1}:{ccd.ccdl2}], ' \
                               f'{cal}=[{fc1}:{fc2},{fl1}:{fl2}]'
                        raise CCDProcError(estr)

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
                    logstr = logstring(ostr, ccd.inim, verbose, logfile)
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
                    hstr = f'[{ccd.biasc1}:{ccd.biasc2},' \
                           f'{ccd.biasl1}:{ccd.biasl2}]'
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
            now = datetime.now()
            now = now.strftime('%Y-%m-%d %H:%M:%S')
            ostr = f'{now} CCD processing done'
            set_header_value(ccd.outim, instrument, 'ccdproc', ostr)
            # end set_header

        # there was a complex problem when input was ushort and output
        # was anything else where the bzero from the input file was being
        # applied to the output and adding 2^15 to all results. This appears
        # to fix it.
        if (np.issubdtype(imin[0].data.dtype, np.unsignedinteger) and
                not np.issubdtype(out[0].data.dtype, np.unsignedinteger)):
            out[0].scale(bscale=1, bzero=0)
        out.close()
        imin.close()

        if not noproc:
            if not ccd.cor:
                os.remove(outim)
            elif replace_input:
                os.replace(outim, image)
            else:
                out_exists = os.path.exists(outreal)
                if out_exists and not overwrite:
                    os.remove(outim)
                    raise OSError(f'File {outreal} already exists and overwrite'
                                  f'=False, so it cannot be overwritten.')
                else:
                    os.replace(outim, outreal)

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
            if replace_input:
                readim = image
            else:
                readim = outreal
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
                    rr = ccd_section(rdatasec, defaults=(1, nc, 1, 1, nl, 1))
                    inc1, inc2, incs, inl1, inl2, inls = rr

                    if (inc1 < 1 or inc2 > nc or inl1 < 1 or inl2 > nl or
                            incs != 1 or inl2 != 1):
                        raise CCDProcError('Error in DATASEC parameter')

                    # The default ccd section is the data section.
                    rccdsec = get_header_value(rim, instrument, 'ccdsec')
                    rr = ccd_section(rccdsec, defaults=(inc1, inc2, incs,
                                                        inl1, inl2, inls))
                    ccdc1, ccdc2, ccdcs, ccdl1, ccdl2, ccdls = rr

                    if ccdcs != 1 or ccdls != 1:
                        raise CCDProcError('Error in CCDSEC parameter')
                    if (inc2 - inc1 != ccdc2 - ccdc1 or
                            inl2 - inl1 != ccdl2 - ccdl1):
                        raise CCDProcError('Size of DATASEC and CCDSEC do not '
                                           'agree')

                    # XXX: make sure the output file deals with the
                    # 'pixeltype' correctly as in set_output though
                    outdata = rim[0].data * 1

                    rtmp = tempfile.NamedTemporaryFile(delete=False)
                    rtmp.close()
                    newout = rtmp.name
                    file_new_copy(newout, rim, instrument=instrument,
                                  overwrite=overwrite)
                    newhdr = image_open(newout, mode='update')
                    rim.close()
                    # Average across the readout axis.

                    # zero out the parts not in the data section
                    outdata[:, inc1 - 1] = 0.
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
                    logstr = logstring(ostr, newout, verbose, logfile)
                    # XXX: can't be what is actually put in the header right?
                    set_header_value(newout, instrument, 'readcor', logstr)

                    newout.close()
            else:
                rim.close()
            # end readcor

        if ccdtype == 'flat':
            if replace_input:
                readim = image
            else:
                readim = outreal
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

                    # IRAF uses 1/1/1980 as its epoch
                    diff = datetime.now() - datetime(1980, 1, 1)
                    set_header_value(meanim, instrument, 'ccdmean', mean)
                    set_header_value(meanim, instrument, 'ccdmeant',
                                     int(diff.total_seconds()),
                                     comment="Time mean was computed "
                                             "(seconds since 1/1/1980)")
            meanim.close()
            # end ccdmean

    # XXX: close all calibration images
