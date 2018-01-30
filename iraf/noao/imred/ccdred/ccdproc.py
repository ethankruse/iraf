from iraf.utils import file_handler
from .combine import Instrument, ccdtypes, ccdsubset, get_header_value
from .combine import file_new_copy, type_max, set_header_value
import numpy as np
import os
from iraf.sys import image_open, image_close
import tempfile
import datetime

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

    """
    def __init__(self):
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

        image_close(image)


def ccd_section(section):
    """
    CCD_SECTION -- Parse a 2D image section into its elements.
    1. The default values must be set by the caller.
    2. A null image section is OK.
    3. The first nonwhitespace character must be '['.
    4. The last interpreted character must be ']'.

    If they are not explicitly given, default values for the first index are 0,
    and default step size is 1. However, None is returned for the last index
    and must be filled in by the caller.

    E.g. [:, 1:20] will return (0, None, 1, 1, 20, 1)

    HOWEVER if nothing is given (input is either None or empty string) then
    the first index will also be None.

    E.g. '  ' will return (None, None, 1, None, None, 1)

    Parameters
    ----------
    section

    Returns
    -------
    Tuple of 6 parameters: first dimension start, stop, step and second
    dimension start, stop, step

    """
    if section is None:
        return None, None, 1, None, None, 1

    section = section.strip()
    if len(section) == 0:
        return None, None, 1, None, None, 1

    # XXX: need to check this function against the IRAF version and see how
    # compatible they are

    if section[0] != '[' and section[-1] != ']':
        raise Exception(f"Error in 2D image section specification {section}")
    # remove the brackets
    osection = section
    section = section[1:-1]
    dims = section.split(',')
    if len(dims) != 2:
        raise Exception(f"Error in 2D image section specification {osection}")
    ret = []
    for dim in dims:
        dim = dim.strip()
        split = dim.split(':')
        if len(split) == 1:
            start = 0
            end = None
            step = 1
        elif len(split) == 2 or len(split) == 3:
            step = 1
            try:
                start = int(split[0])
            except ValueError:
                start = 0
            try:
                end = int(split[1])
            except ValueError:
                end = None
            if len(split) == 3:
                try:
                    step = int(split[2])
                except ValueError:
                    step = 1
        else:
            raise Exception(
                f"Error in 2D image section specification {osection}")
        ret.append(start)
        ret.append(end)
        ret.append(step)
    return ret


def ccdproc(images, output, *, ccdtype='object', noproc=False, fixpix=True,
            overscan=True, trim=True, zerocor=True, darkcor=True, flatcor=True,
            illumcor=False, fringecor=False, readcor=False, scancor=False,
            readaxis='line', fixfile=None, biassec=None, trimsec=None,
            zero=None, dark=None, flat=None, illum=None, fringe=None,
            minreplace=1., scantype='shortscan', nscan=1, interactive=False,
            function='legendre', order=1, sample='*', naverage=1, niterate=1,
            low_reject=3., high_reject=3., grow=0., instrument=None,
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
    logfile
    verbose

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
        datasec = get_header_value(ccd.inim, instrument, 'datasec')
        c1, c2, cs, l1, l2, ls = ccd_section(datasec)
        # XXX: should these be nc or nc - 1. also check upper bounds on
        # exception below
        if c1 is None:
            c1 = 0
        if c2 is None:
            c2 = nc
        if l1 is None:
            l1 = 0
        if l2 is None:
            l2 = nl
        if c1 < 0 or c2 > nc or cs != 1 or l1 < 0 or l2 > nl or ls != 1:
            raise Exception(f"Error in datasec parameter: {datasec}")

        ccd.inc1 = c1
        ccd.inc2 = c2
        ccd.inl1 = l1
        ccd.inl2 = l2

        # The default trim section is the data section.
        # Defer limit checking until actually used.
        ts = trimsec
        if trimsec == 'image':
            ts = get_header_value(ccd.inim, instrument, 'trimsec')
        c1, c2, cs, l1, l2, ls = ccd_section(ts)
        if c1 is None:
            c1 = ccd.inc1
        if c2 is None:
            c2 = ccd.inc2
        if l1 is None:
            l1 = ccd.inl1
        if l2 is None:
            l2 = ccd.inl2
        if cs != 1 or ls != 1:
            raise Exception(f"Error in trimsec parameter: {ts}")

        ccd.trimc1 = c1
        ccd.trimc2 = c2
        ccd.triml1 = l1
        ccd.triml2 = l2

        # The default bias section is the whole image.
        bs = biassec
        if biassec == 'image':
            bs = get_header_value(ccd.inim, instrument, 'biassec')
        c1, c2, cs, l1, l2, ls = ccd_section(bs)
        # XXX: should these be nc or nc - 1.
        if c1 is None:
            c1 = 0
        if c2 is None:
            c2 = nc
        if l1 is None:
            l1 = 0
        if l2 is None:
            l2 = nl
        if cs != 1 or ls != 1:
            raise Exception(f"Error in biassec parameter: {bs}")

        ccd.biasc1 = c1
        ccd.biasc2 = c2
        ccd.biasl1 = l1
        ccd.biasl2 = l2

        # The default ccd section is the size of the data section.
        ccs = get_header_value(ccd.inim, instrument, 'ccdsec')
        c1, c2, cs, l1, l2, ls = ccd_section(ccs)
        if c1 is None:
            c1 = 0
        # XXX: check this size
        if c2 is None:
            c2 = ccd.inc2 - ccd.inc1 + 1
        if l1 is None:
            l1 = 0
        if l2 is None:
            l2 = ccd.inl2 - ccd.inl1 + 1
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
        # XXX: ccdflag (IN_IM(ccd), "trim")
        if trim and False:
            # Check trim section.
            # XXX: check these upper limits
            if (ccd.trimc1 < 0 or ccd.trimc2 > nc or ccd.triml1 < 0 or
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

                ccd.outc1 = ccd.inc1 - ccd.trimc1 + 1
                ccd.outc2 = ccd.inc2 - ccd.trimc1 + 1
                ccd.outl1 = ccd.inl1 - ccd.triml1 + 1
                ccd.outl2 = ccd.inl2 - ccd.triml1 + 1

                ccd.cors['trim'] = True
                ccd.cor = True

                ostr = f"Trim data section is [{ccd.trimc1:d}:{ccd.trimc2:d}," \
                       f"{ccd.triml1:d}:{ccd.triml2:d}]"

                now = datetime.datetime.now()
                now = now.strftime('%Y-%m-%d %H:%M:%S')

                ostr = f'{now} {ostr}'
                logstr = f'{ccd.inim.filename()}: {ostr}\n'
                if verbose:
                    print(logstr)
                if logfd is not None:
                    logfd.write(logstr)
                # XXX: this can't be what is actually put in the header right?
                set_header_value(ccd.inim, instrument, 'trim', logstr)

        # end set_trim

        image_close(imin)
        image_close(out)
        if outtmp:
            # XXX: will eventually need to copy to input parameters
            os.remove(outim)
    if logfd is not None:
        logfd.close()
