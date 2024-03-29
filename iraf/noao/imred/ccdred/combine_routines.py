import os

import numpy as np

from iraf.sys import image_open
from iraf.utils import file_handler
from . import Instrument
from .utils import CCDProcError, ccdsubset, ccdtypes, file_new_copy
from .utils import get_header_value, set_header_value, type_max

__all__ = ['combine']


def ic_setout(inputs, output, nimages, project, offsets):
    """
    Set output image size and offsets of input images.

    Parameters
    ----------
    inputs
    output
    nimages
    project
    offsets

    Returns
    -------

    """

    indim = len(inputs[0][0].data.shape)
    outdim = len(output[0][0].data.shape)

    if project:
        outdim = indim - 1
    else:
        for im in inputs:
            if len(im[0].data.shape) != outdim:
                raise CCDProcError("Input image dimensions are not the same")

    # Parse the user offset string.  If "none" then there are no offsets.
    # If "wcs" then set the offsets based on the image WCS.
    # If "grid" then set the offsets based on the input grid parameters.
    # If a file scan it.
    if offsets == 'none':
        offsetsarr = np.zeros((nimages, outdim), dtype=int)
        reloff = True
    # XXX: implement these
    elif offsets.lower() == 'wcs':
        raise NotImplementedError('WCS offsets not implemented yet.')
    elif offsets.lower() == 'grid':
        raise NotImplementedError('grid offsets not implemented yet.')
    else:
        raise NotImplementedError('Manual file offsets not implemented yet.')

    aligned = True

    for jj in np.arange(outdim):
        aa = offsetsarr[0, jj]
        bb = inputs[0][0].data.shape[jj] + aa
        amin = aa
        bmax = bb
        for ii in np.arange(nimages - 1) + 1:
            aa = offsetsarr[ii, jj]
            bb = inputs[ii][0].data.shape[jj] + aa
            if aa != amin or bb != bmax or not reloff:
                aligned = False
            amin = min(aa, amin)
            bmax = max(bb, bmax)
        # XXX: also do? IM_LEN(out[1],j) = bmax
        if reloff or amin < 0:
            offsetsarr[:, jj] -= amin

            # XXX: also do? IM_LEN(out[1],j) = IM_LEN(out[1],j) - amin

    # Update the WCS.
    # XXX: do this
    if project or not aligned or not reloff:
        raise NotImplementedError("WCS updates to output files not "
                                  "implemented yet!")

    return aligned, offsetsarr


def ic_mopen(in_images, mtype, mvalue, instrument):
    """
    Open pixel mask files.

    Parameters
    ----------
    in_images
    mtype
    mvalue
    instrument

    Returns
    -------

    """
    mtype = mtype.strip().lower()

    if mtype not in ['none', 'goodvalue', 'badvalue', 'goodbits', 'badbits']:
        raise ValueError(f'masktype {mtype} not recognized.')

    # Check for special cases.  The BOOLEAN type is used when only
    # zero and nonzero are significant; i.e. the actual mask values are
    # not important.  The invert flag is used to indicate that
    # empty masks are all bad rather than all good.
    if mtype == 'badbits' and mvalue == 0:
        mtype = 'none'
    if mvalue == 0 and mtype in ['goodvalue', 'goodbits']:
        mtype = 'boolean'
    if ((mvalue == 0 and mtype in ['badvalue', 'goodbits']) or
            (mvalue != 0 and mtype == 'goodvalue')):
        invert = True
    else:
        invert = False

    # If mask images are to be used, get the mask name from the image
    # header and open it saving the descriptor in the pms array.
    # Empty masks (all good) are treated as if there was no mask image.
    npms = 0
    pms = []

    for im in in_images:
        if mtype != 'none':
            fname = get_header_value(im, instrument, 'BPM')
            if fname is None:
                pms.append(None)
                continue

            pm = image_open(fname)
            if np.allclose(pm[0].data, 0) and not invert:
                pm.close()
                pm = None
            else:
                npms += 1
            pms.append(pm)
        else:
            pms.append(None)
    # If no mask images are found and the mask parameters imply that
    # good values are 0 then use the special case of no masks.
    if npms == 0 and not invert:
        mtype = 'none'

    return mtype, pms


def combine(images, output, *, plfile=None, sigmafile=None, ccdtype=None,
            subsets=False, delete=False, method='average',
            reject='none', project=False, outtype=None, offsets='none',
            masktype='none', maskvalue=0., blank=0., scale=None, zero=None,
            weight=None, statsec=None, lthreshold=None, hthreshold=None,
            nlow=1, nhigh=1, nkeep=1, mclip=True, lsigma=3.0,
            hsigma=3.0, rdnoise=0., gain=1., snoise=0., sigscale=0.1,
            pclip=-0.5, grow=0, instrument=None, logfile=None, verbose=False):
    """

    Parameters
    ----------
    images :
        List of images to combine
    output : str
        Base output file name
    plfile :
        List of output pixel list files (optional)
    sigmafile :
        List of sigma images (optional)
    ccdtype  :
        CCD image type to combine (optional)
    subsets :
        Combine images by subset parameter?
    delete :
        Delete input images after combining?
    method : "average|median"
        Type of combine operation
    reject : "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip"
        Type of rejection
    project :
        Project highest dimension of input images?
    outtype :
        Output image pixel datatype (e.g. np.float64).
        Default None gives the output image the datatype of the inputs.
        Will also revert to the input type if outtype is too low.
        E.g. inputs have int64 and outtype is int32.
    offsets :
        Input image offsets. 'none, wcs, grid, or file name'
    masktype : "none|goodvalue|badvalue|goodbits|badbits"
        Mask type
    maskvalue :
        Mask value
    blank :
        Value if there are no pixels
    scale :
        Image scaling
    zero :
        Image zero point offset
    weight :
        Image weights
    statsec :
        Image section for computing statistics
    lthreshold :
        Lower threshold
    hthreshold :
        Upper threshold
    nlow :
        minmax: Number of low pixels to reject
    nhigh :
        minmax: Number of high pixels to reject
    nkeep :
        Minimum to keep (pos) or maximum to reject (neg)
    mclip :
        Use median in sigma clipping algorithms?
    lsigma :
        Lower sigma clipping factor
    hsigma :
        Upper sigma clipping factor
    rdnoise :
        ccdclip: CCD readout noise (electrons).
        Either a number or string of the keyword to pull the number
        from the image header file.
    gain :
        ccdclip: CCD gain (electrons/DN)
        Either a number or string of the keyword to pull the number
        from the image header file.
    snoise :
        ccdclip: Sensitivity noise (fraction)
        Either a number or string of the keyword to pull the number
        from the image header file.
    sigscale :
        Tolerance for sigma clipping scaling corrections
    pclip :
        Percentile clipping parameter
    grow :
        Radius (pixels) for 1D neighbor rejection
    instrument
    logfile
    verbose :
        Print log information to the standard output?

    Returns
    -------

    Rejection methods:
    none: all pixels used.

    minmax: nhigh and nlow pixels are rejected no matter what.
    Can be either fractions of number of images (nhigh < 1) or
    integer number of images to remove.

    ccdclip

    """
    # basic error checking before we start manipulating files
    method = method.strip().lower()
    if method not in ['average', 'median']:
        raise ValueError(f'Combine method not recognized: {method}')

    reject = reject.lower().strip()
    rejectopts = "none|ccdclip|crreject|minmax|pclip|sigclip|avsigclip"
    if reject not in rejectopts.split('|'):
        raise ValueError(f'Could not recognize reject parameter {reject}')

    # was given a string or something else, so set up the instrument object
    if not isinstance(instrument, Instrument):
        instrument = Instrument(instrument)

    # start of IRAF cmb_images.

    inputs = file_handler(images)

    if len(inputs) == 0:
        return

    # Go through the input list and eliminate images not satisfying the
    # CCD image type.  Separate into subsets if desired.  Create image
    # and subset lists.

    # lists of images in each subset
    images = []
    # subset names
    subset = []

    for image in inputs:
        hdulist = image_open(image)
        if hdulist is None:
            continue

        # what image type is this one
        thistype = ccdtypes(hdulist, instrument)

        # if this isn't the image type we're looking for, skip it
        if ccdtype is not None and thistype != ccdtype:
            hdulist.close()
            continue

        if subsets:
            subsetstr = ccdsubset(hdulist, instrument)
        else:
            subsetstr = ''

        # add the file to the correct subset list
        if subsetstr not in subset:
            subset.append(subsetstr)
            images.append([image])
        else:
            images[subset.index(subsetstr)].append(image)

        hdulist.close()

    # end of IRAF cmb_images

    if len(images) == 0:
        print("No images to combine.")
        return

    # set threshold flag
    if lthreshold is None and hthreshold is None:
        dothresh = False
    else:
        dothresh = True
        if lthreshold is None:
            lthreshold = -np.inf
        if hthreshold is None:
            hthreshold = np.inf

    # get copies of these before adding to them for each subset
    outroot = output.strip()
    if plfile is not None:
        plroot = plfile.strip()
    else:
        plroot = None
    if sigmafile is not None:
        sigroot = sigmafile.strip()
    else:
        sigroot = None

    # ignore NaN operation warnings
    olderr = np.seterr(invalid='ignore')

    # don't overwrite these values for subsequent subsets
    origoffsets = offsets
    origmasktype = masktype
    origreject = reject
    orignlow = nlow * 1
    orignhigh = nhigh * 1
    orignkeep = nkeep * 1
    origpclip = pclip * 1

    # Combine each input subset.
    for zz, iset in enumerate(subset):
        # restore these values to the inputs
        offsets = origoffsets
        masktype = origmasktype
        reject = origreject
        nlow = orignlow * 1
        nhigh = orignhigh * 1
        nkeep = orignkeep * 1
        pclip = origpclip * 1

        iimages = images[zz]

        # this apparently isn't done in IRAF, but I like separating the two
        # parts with a '.' if there is a subset string we're appending
        if len(iset) > 0:
            comb = '.'
        else:
            comb = ''

        base, ext = os.path.splitext(outroot)
        output = f'{base}{comb}{iset}{ext}'

        if plroot is not None:
            base, ext = os.path.splitext(plroot)
            plfile = f'{base}{comb}{iset}{ext}'
        else:
            plfile = None

        if sigroot is not None:
            base, ext = os.path.splitext(sigroot)
            sigmafile = f'{base}{comb}{iset}{ext}'
        else:
            sigmafile = None

        # icombine starts here

        # Set number of images to combine.
        if project:
            if len(iimages) > 1:
                print("Cannot project combine a list of images")
                return
            hdulist = image_open(iimages[0])
            shp = np.array(hdulist[0].data.shape)
            hdulist.close()
            if shp.size == 1 or shp[shp > 1].size == 1:
                print("Can't project one dimensional images")
                return
            # XXX: which dimension do we project over?
            nimages = shp[-1]
        else:
            nimages = len(iimages)

        # Convert the nkeep parameter if needed.
        if nkeep < 0:
            nkeep = max(0, nimages + nkeep)

        # Convert the pclip parameter to a number of pixels rather than
        # a fraction.  This number stays constant even if pixels are
        # rejected.  The number of low and high pixel rejected, however,
        # are converted to a fraction of the valid pixels.

        # REJECT	"|none|ccdclip|crreject|minmax|pclip|sigclip|avsigclip|"
        if reject == 'pclip':
            if pclip == 0.:
                estr = "pclip parameter can not be zero when reject=='pclip'"
                raise ValueError(estr)

            ii = nimages // 2
            if np.abs(pclip) < 1.:
                pclip *= ii
            # pclip has to be at least 1 image away from the median value
            if pclip < 0.:
                pclip = min(-1, max(-ii, int(pclip)))
            else:
                pclip = max(1, min(ii, int(pclip)))

        if reject == 'minmax':
            if nlow >= 1.:
                nlow /= nimages
            if nhigh >= 1.:
                nhigh /= nimages

            ii = nlow * nimages
            jj = nhigh * nimages
            if ii + jj == 0:
                reject = 'none'
            elif ii + jj >= nimages:
                raise ValueError(f"Bad minmax rejection parameters: {nlow:.2f}"
                                 f" + {nhigh:.2f} > 1")

        imin = []
        # Map the input image(s).
        for im in iimages:
            tmp = image_open(im)
            imin.append(tmp)

        out = []

        # Map the output image and set dimensions and offsets.
        file_new_copy(output, imin[0], mode='NEW_COPY', overwrite=True,
                      instrument=instrument)
        out.append(image_open(output, mode='update'))

        aligned, offarr = ic_setout(imin, out, nimages, project, offsets)

        # Determine the highest precedence datatype and set output datatype.
        intype = imin[0][0].data.dtype
        for im in imin:
            intype = type_max(intype, im[0].data.dtype)

        # don't overwrite input outtype for subsequent subsets
        oouttype = outtype
        if oouttype is None:
            oouttype = intype
        else:
            oouttype = oouttype.strip().lower()
            otypes = "short|ushort|integer|long|real|double".split('|')
            ndtypes = [np.short, np.ushort, np.int_, np.int_,
                       np.single, np.double]
            if oouttype in otypes:
                oouttype = ndtypes[otypes.index(oouttype)]
            else:
                raise ValueError(f'Unrecognized outtype: {oouttype}')
        # make sure the output type will work given the input
        oouttype = type_max(intype, oouttype)

        # Open pixel list file if given.
        if plfile is not None:
            file_new_copy(plfile, out[0], mode='NEW_COPY', overwrite=True,
                          instrument=instrument)
            out.append(image_open(plfile, mode='update'))
        else:
            out.append(None)

        sigmatype = None
        # Open the sigma image if given.
        if sigmafile is not None:
            file_new_copy(sigmafile, out[0], mode='NEW_COPY', overwrite=True,
                          instrument=instrument)
            out.append(image_open(sigmafile, mode='update'))
            # has to be a float
            sigmatype = type_max(np.float, oouttype)
        else:
            out.append(None)

        # Open masks.
        masktype, pms = ic_mopen(imin, masktype, maskvalue, instrument)

        # Open the log file.
        logfd = None
        if logfile is not None:
            logfd = open(logfile, 'a')

        # this is where the icombiner function starts

        # icombiner seems to just be a bunch of memory handling to see if the
        # computer will run out of memory. We're skipping this and assuming
        # memory management will be taken care of behind the scenes. Move
        # to ic_combiner.

        # beginning of ic_scale

        # Set the defaults.
        ncombine = np.ones(nimages).astype(int)
        exptime = np.zeros(nimages)
        means = np.zeros(nimages)
        medians = np.zeros(nimages)
        modes = np.zeros(nimages)
        scales = np.ones(nimages)
        zeros = np.zeros(nimages)
        wts = np.ones(nimages)

        # Get the number of images previously combined and the exposure times.
        # The default combine number is 1 and the default exposure is 0.
        for ii, im in enumerate(imin):
            nc = get_header_value(im, instrument, 'ncombine')
            et = get_header_value(im, instrument, 'exptime')
            if nc is None:
                nc = 1
            if et is None:
                et = 0.
            ncombine[ii] = nc
            exptime[ii] = et
            # all the same image, so repeat these values
            if project:
                ncombine[:] = nc
                exptime[:] = et
                break

        # Set scaling factors.
        # all can also be '@file' or '!keyword' which turns into
        # type == 'file' or 'keyword'
        stypes = "none|mode|median|mean|exposure".split('|')
        ztypes = "none|mode|median|mean".split('|')
        wtypes = "none|mode|median|mean|exposure".split('|')

        stype = ic_gscale(scale, stypes, imin, exptime, scales, nimages,
                          instrument, project)
        ztype = ic_gscale(zero, ztypes, imin, exptime, zeros, nimages,
                          instrument, project)
        wtype = ic_gscale(weight, wtypes, imin, exptime, wts, nimages,
                          instrument, project)

        # Get image statistics only if needed.
        domode = 'mode' in [stype, ztype, wtype]
        domedian = 'median' in [stype, ztype, wtype]
        domean = 'mean' in [stype, ztype, wtype]

        # statsec options: "|input|output|overlap|"
        if statsec is None:
            statsec = ''
        statsec = statsec.strip().lower()

        if domode or domedian or domean:
            # oimref stuff only needed if it needs to be an argument to
            # ic_stat again
            # oimref = None
            section = statsec
            if statsec == 'input':
                section = ''
            elif statsec == 'output':
                section = ''
                # oimref = out[0]
            elif statsec == 'overlap':
                section = '['
                for ii in np.arange(out[0][0].data.ndim):
                    kk = offarr[0, ii]
                    ll = offarr[0, ii] + imin[0][0].data.shape[ii]
                    for jj in np.arange(1, nimages):
                        kk = max(kk, offarr[jj, ii])
                        ll = min(ll,
                                 offarr[jj, ii] + imin[jj][0].data.shape[ii])
                    section += f'{kk:d}:{ll:d},'
                section = section[:-1]
                section += ']'
                # oimref = out[0]

            for ii in np.arange(nimages):
                mn, md, mode = ic_stat(imin[ii], section, offarr,
                                       project, ii, masktype, dothresh,
                                       lthreshold, hthreshold, domode=domode)
                if domode:
                    modes[ii] = mode
                medians[ii] = md
                means[ii] = mn

                dc = {'mean': mn, 'median': md, 'mode': mode}
                if stype in dc:
                    scales[ii] = dc[stype]
                if ztype in dc:
                    zeros[ii] = dc[ztype]
                if wtype in dc:
                    wts[ii] = dc[wtype]

        if (scales <= 0.).any():
            print("WARNING: Negative scale factors -- ignoring scaling.")
            scales = np.ones(nimages)

        snorm = False
        if stype in ['file', 'keyword']:
            snorm = True
        znorm = False
        if ztype in ['file', 'keyword']:
            znorm = True
        wflag = False
        if wtype in ['file', 'keyword']:
            wflag = True

        if snorm:
            scales = 1. / scales
        else:
            scales /= scales.mean()

        zeros /= scales
        zmean = zeros.mean()

        bad = False
        if wtype != 'none':
            if (wts <= 0.).any():
                em = "WARNING: Negative weights"
                em += " -- using only NCOMBINE weights."
                print(em)
                wts = ncombine * 1
                bad = True
            if (zeros <= 0.).any():
                em = "WARNING: Negative zero offsets"
                em += " -- ignoring zero weight adjustments."
                print(em)
                wts = ncombine * wts
                bad = True
            if not bad:
                if ztype == 'none' or znorm or wflag:
                    wts = ncombine * wts
                else:
                    wts = ncombine * wts * zmean / zeros

        if znorm:
            zeros *= -1
        else:
            zeros -= zmean
            # Because of finite arithmetic it is possible for the zero offsets
            # to be nonzero even when they are all equal.  Just for the sake of
            # a nice log message set the zero offsets in this case.
            allclose = True
            for ii in np.arange(zeros.size):
                if not np.isclose([zeros[ii]], [0.]):
                    allclose = False
            if allclose:
                zeros = np.zeros(nimages)

        # normalize the weights to sum to 1
        wts /= wts.sum()

        # Set flags for scaling, zero offsets, sigma scaling, weights.
        # Sigma scaling may be suppressed if the scales or zeros are
        # different by a specified tolerance.
        doscale = False
        dozero = False
        doscale1 = False
        dowts = False
        for ii in np.arange(nimages - 1) + 1:
            if snorm or scales[ii] != scales[0]:
                doscale = True
            if znorm or zeros[ii] != zeros[0]:
                dozero = True
            if wts[ii] != wts[0]:
                dowts = True

        if doscale and sigscale != 0.:
            if (np.abs(scales - 1) > sigscale).any():
                doscale1 = True
            if not doscale1 and zmean > 0.:
                if (np.abs(zeros / zmean) > sigscale).any():
                    doscale1 = True

        # Set the output header parameters.
        nout = ncombine.sum()
        set_header_value(out[0], instrument, 'ncombine', nout)
        exposure = (wts * exptime / scales).sum()
        darktime = 0.
        mean = 0.
        for ii in np.arange(nimages):
            dark = get_header_value(imin[ii], instrument, 'darktime')
            if dark is not None:
                darktime += wts[ii] * dark / scales[ii]
            else:
                darktime += wts[ii] * exptime[ii] / scales[ii]
            mode = get_header_value(imin[ii], instrument, 'ccdmean')
            if mode is not None:
                mean += wts[ii] * mode / scales[ii]

        set_header_value(out[0], instrument, 'exptime', exposure)
        set_header_value(out[0], instrument, 'darktime', darktime)
        exists = get_header_value(out[0], instrument, 'ccdmean')
        if exists is not None:
            set_header_value(out[0], instrument, 'ccdmean', mean)

        if out[1] is not None:
            oname = out[1].filename()
            set_header_value(out[0], instrument, 'BPM', oname)

        # Start the log here since much of the info is only available here.
        if verbose or logfd is not None:
            import datetime

            stack = False
            if project:
                fname = get_header_value(imin[0], instrument, 'stck0001')
                if fname is not None:
                    stack = True

            # Time stamp the log and print parameter information.
            logtxt = ''
            now = datetime.datetime.now()
            now = now.strftime('%Y-%m-%d %H:%M:%S')
            logtxt += f'{now} : IMCOMBINE\n'
            logtxt += f'  combine = {method}, '
            logtxt += f'scale = {scale}, zero = {zero}, '
            logtxt += f'weight = {weight}\n'

            # REJECT "|none|ccdclip|crreject|minmax|pclip|sigclip|avsigclip|"
            if reject == 'minmax':
                ostr = '  reject = minmax, nlow = {0:d}, nhigh = {1:d}\n'
                slow = int(np.round(nlow * nimages))
                shigh = int(np.round(nhigh * nimages))
                logtxt += ostr.format(slow, shigh)
            elif reject == 'ccdclip':
                ostr = '  reject = ccdclip, mclip = {0}, nkeep = {1:d}\n'
                logtxt += ostr.format(mclip, nkeep)
                if isinstance(rdnoise, str):
                    ostr = f'  rdnoise = {rdnoise}'
                else:
                    ostr = f'  rdnoise = {rdnoise:g}'
                if isinstance(gain, str):
                    ostr += f', gain = {gain}, '
                else:
                    ostr += f', gain = {gain:g}, '
                if isinstance(snoise, str):
                    ostr += f'snoise = {snoise}, '
                else:
                    ostr += f'snoise = {snoise:g}, '
                ostr += f'lsigma = {lsigma:g}, hsigma = {hsigma:g}\n'
                logtxt += ostr
            elif reject == 'crreject':
                ostr = '  reject = crreject, mclip = {0}, nkeep = {1:d}\n'
                logtxt += ostr.format(mclip, nkeep)
                if isinstance(rdnoise, str):
                    ostr = f'  rdnoise = {rdnoise}'
                else:
                    ostr = f'  rdnoise = {rdnoise:g}'
                if isinstance(gain, str):
                    ostr += f', gain = {gain}, '
                else:
                    ostr += f', gain = {gain:g}, '
                if isinstance(snoise, str):
                    ostr += f'snoise = {snoise}, '
                else:
                    ostr += f'snoise = {snoise:g}, '
                ostr += f'hsigma = {hsigma:g}\n'
                logtxt += ostr
            elif reject == 'pclip':
                ostr = f'  reject = pclip, nkeep = {nkeep:d}\n'
                logtxt += ostr
                ostr = '  pclip = {0:g}, lsigma = {1:g}, hsigma = {2:g}\n'
                logtxt += ostr.format(pclip, lsigma, hsigma)
            elif reject == 'sigclip':
                ostr = '  reject = sigclip, mclip = {0}, nkeep = {1:d}\n'
                logtxt += ostr.format(mclip, nkeep)
                ostr = '  lsigma = {0:g}, hsigma = {1:g}\n'
                logtxt += ostr.format(lsigma, hsigma)
            elif reject == 'avsigclip':
                ostr = '  reject = avsigclip, mclip = {0}, nkeep = {1:d}\n'
                logtxt += ostr.format(mclip, nkeep)
                ostr = '  lsigma = {0:g}, hsigma = {1:g}\n'
                logtxt += ostr.format(lsigma, hsigma)

            if reject != 'none' and grow > 0:
                logtxt += f'  grow = {grow:d}\n'

            if dothresh:
                if lthreshold > -np.inf and hthreshold < np.inf:
                    ostr = '  lthreshold = {0:g}, hthreshold = {1:g}\n'
                    ostr = ostr.format(lthreshold, hthreshold)
                elif lthreshold > -np.inf:
                    ostr = f'  lthreshold = {lthreshold:g}\n'
                else:
                    ostr = f'  hthreshold = {hthreshold:g}\n'
                logtxt += ostr

            logtxt += f'  blank = {blank:g}\n'
            if len(statsec) > 0:
                logtxt += f'  statsec = {statsec}\n'

            if masktype != 'none':
                raise NotImplementedError('Mask types not yet supported')

            dop = {'ncombine': False, 'exptime': False, 'mode': False,
                   'median': False, 'mean': False, 'mask': False, 'rdn': False,
                   'gain': False, 'sn': False}
            for ii in np.arange(nimages):
                if ncombine[ii] != ncombine[0]:
                    dop['ncombine'] = True
                if exptime[ii] != exptime[0]:
                    dop['exptime'] = True
                if modes[ii] != modes[0]:
                    dop['mode'] = True
                if medians[ii] != medians[0]:
                    dop['median'] = True
                if means[ii] != means[0]:
                    dop['mean'] = True
                if masktype != 'none':
                    raise NotImplementedError('Mask types not yet supported')
            if reject == 'ccdclip' or reject == 'crreject':
                if isinstance(rdnoise, str):
                    dop['rdn'] = True
                if isinstance(gain, str):
                    dop['gain'] = True
                if isinstance(snoise, str):
                    dop['sn'] = True
            hdr = '  {0:20s} '.format('Images')
            if dop['ncombine']:
                hdr += ' {0:6s}'.format('NComb')
            if dop['exptime']:
                hdr += ' {0:6s}'.format('ExpT')
            if dop['mode']:
                hdr += ' {0:7s}'.format('Mode')
            if dop['median']:
                hdr += ' {0:7s}'.format('Median')
            if dop['mean']:
                hdr += ' {0:7s}'.format('Mean')
            if dop['rdn']:
                hdr += ' {0:7s}'.format('RdNoise')
            if dop['gain']:
                hdr += ' {0:6s}'.format('Gain')
            if dop['sn']:
                hdr += ' {0:6s}'.format('SNoise')
            if doscale:
                hdr += ' {0:6s}'.format('Scale')
            if dozero:
                hdr += ' {0:7s}'.format('Zero')
            if dowts:
                hdr += ' {0:6s}'.format('Weight')
            if not aligned:
                hdr += ' {0:9s}'.format('Offsets')
            if dop['mask']:
                hdr += ' {0:s}'.format('Maskfile')
            hdr += '\n'
            logtxt += hdr

            for ii in np.arange(nimages):
                if stack:
                    stc = f'stck{ii:04d}'
                    vl = get_header_value(imin[ii], instrument, stc)
                    if vl is not None:
                        ostr = f'  {vl:21s}'
                    else:
                        ostr = f'  {imin[ii].filename():16s}[{ii:03d}]'
                elif project:
                    ostr = f'  {imin[ii].filename():16s}[{ii:03d}]'
                else:
                    ostr = f'  {imin[ii].filename():21s}'
                if dop['ncombine']:
                    ostr += f' {ncombine[ii]:6d}'
                if dop['exptime']:
                    ostr += f' {exptime[ii]:6.1f}'
                if dop['mode']:
                    ostr += f' {modes[ii]:7.5g}'
                if dop['median']:
                    ostr += f' {medians[ii]:7.5g}'
                if dop['mean']:
                    ostr += f' {means[ii]:7.5g}'
                if dop['rdn']:
                    rval = get_header_value(imin[ii], instrument, rdnoise)
                    ostr += f' {rval:7g}'
                if dop['gain']:
                    rval = get_header_value(imin[ii], instrument, gain)
                    ostr += f' {rval:6g}'
                if dop['sn']:
                    rval = get_header_value(imin[ii], instrument, snoise)
                    ostr += f' {rval:6g}'
                if doscale:
                    ostr += f' {1./scales[ii]:6.3f}'
                if dozero:
                    ostr += f' {-zeros[ii]:7.5g}'
                if dowts:
                    ostr += f' {wts[ii]:6.3f}'
                if not aligned:
                    nd = np.array(out[0][0].data.shape)
                    # number of dimensions in out array
                    nd = len(nd[nd > 1])
                    if nd == 1:
                        ostr += f' {offarr[ii, 0]:9d}'
                    else:
                        use = offarr[ii, :]
                        strs = [f' {ixx:4d}' for ixx in use]
                        for istr in strs:
                            ostr += istr
                if dop['mask']:
                    raise NotImplementedError('Mask types not yet supported')
                ostr += '\n'
                logtxt += ostr

            # Log information about the output images.
            ostr = '\n  Output image = {0}, ncombine = {1:d}\n'
            ostr = ostr.format(out[0].filename(), nout)
            logtxt += ostr

            if out[1] is not None:
                logtxt += f'  Pixel list image = {out[1].filename()}\n'

            if out[2] is not None:
                logtxt += f'  Sigma image = {out[2].filename()}\n'

            if verbose:
                print(logtxt)
            if logfd is not None:
                logfd.write(logtxt)

        doscale = (doscale or dozero)
        # end of the icscale function

        # Set rejection algorithm specific parameters
        # the column order is readnoise, gain, snoise
        nm = np.zeros((nimages, 3))
        if reject in ['ccdclip', 'crreject']:
            if isinstance(rdnoise, str):
                for ii in np.arange(nimages):
                    nm[ii, 0] = get_header_value(imin[ii], instrument, rdnoise)
            else:
                nm[:, 0] = rdnoise

            if isinstance(gain, str):
                for ii in np.arange(nimages):
                    nm[ii, 1] = get_header_value(imin[ii], instrument, gain)
            else:
                nm[:, 1] = gain
            max_real = 0.99e37
            small = 1e4 / max_real
            # adjust the readnoise values to be (readnoise / gain)**2
            # which is what is actually used in the sigma calcluation.
            # don't let it be 0.
            rdsq = (nm[:, 0] / nm[:, 1]) ** 2
            nm[:, 0] = np.where(rdsq > small, rdsq, [small])

            if isinstance(snoise, str):
                for ii in np.arange(nimages):
                    nm[ii, 2] = get_header_value(imin[ii], instrument, snoise)
            else:
                nm[:, 2] = snoise

            if reject == 'crreject':
                lsigma = np.inf
        elif reject == 'none':
            grow = 0

        # start of ic_gdatar
        oshp = out[0][0].data.shape
        oshp += (nimages,)
        # do all calculations at the highest possible precision of
        # output data types, then downgrade later if requested
        data = np.zeros(oshp, dtype=np.double)
        # set bad values
        data[:] = np.nan
        # easy case, no offsets to deal with
        if aligned:
            for ii in np.arange(nimages):
                data[..., ii] = imin[ii][0].data
        else:
            raise NotImplementedError('need to handle offsets in combine')

        # remove the data points with bad mask values (set to nan)
        if masktype != 'none':
            raise NotImplementedError('need to implement mask types')

        # Apply threshold if needed
        if dothresh:
            nothr = np.where((data < lthreshold) | (data > hthreshold))
            data[nothr] = np.nan

        # Apply scaling (avoiding masked pixels which might overflow?)
        if doscale:
            # XXX: this seems to be the wrong order. bug in ic_gdatar?
            data /= scales
            data -= zeros

        # end of ic_gdatar
        if reject in ['ccdclip', 'crreject', 'sigclip', 'avsigclip']:
            minclip = 2
            if reject in ['sigclip', 'avsigclip']:
                minclip = 3
            # number of good points to use
            npts = (~np.isnan(data)).sum(axis=-1)
            minkeep = np.where(npts < nkeep, npts, [nkeep])

            finsig = np.zeros_like(np.sum(npts, axis=-1), dtype=float)
            if reject == 'avsigclip':
                if mclip:
                    meds = np.nanmedian(data, axis=-1)
                else:
                    totals = np.nansum(data, axis=-1)
                    # if there are 3+ data points, use the min/max removed avg
                    adjnpts = npts * 1
                    # this adjustment method will fail if there are infinities
                    totals[adjnpts >= 3] -= np.nanmin(data[adjnpts >= 3],
                                                      axis=-1)
                    totals[adjnpts >= 3] -= np.nanmax(data[adjnpts >= 3],
                                                      axis=-1)
                    adjnpts[adjnpts >= 3] -= 2
                    meds = totals / adjnpts
                # Compute the poisson scaled average sigma about the median.
                # There must be at least three pixels at each point to define
                # the mean sigma.  Corrections for differences in the image
                # scale factors are selected by the doscale1 flag.
                if doscale1:
                    rr = (meds[..., np.newaxis] + zeros) / scales
                    rr[rr < 1.] = 1.
                    ss = (data - meds[..., np.newaxis]) ** 2. / rr
                else:
                    rr = np.where(meds > 1., meds, [1.])
                    ss = ((data - meds[..., np.newaxis]) ** 2. /
                          rr[..., np.newaxis])
                ss = np.nansum(ss, axis=-1)
                n2 = npts >= 3
                # don't contribute pixels with less than 3 data points to the
                # sigma calculations
                clippts = npts * 1
                clippts[~n2] = 0
                ss[~n2] = 0
                # Here is the final sigma.
                # in IRAF this sum is only done along each "row", ie along
                # the first dimension, since that's what ic_gdatar and impnlr
                # deal with.
                # NOTE: the sum is along NAXIS=1 in IRAF, which is column-order,
                # so in numpy which is row-order it is the last dimension.
                if n2.any():
                    finsig = np.sqrt(np.sum(ss, axis=-1) /
                                     (np.sum(clippts, axis=-1) - 1))
                else:
                    # don't do any clipping, so effectively ignore everything
                    finsig = np.ones_like(np.sum(npts, axis=-1), dtype=float)
                    finsig *= np.inf
                # avoid divide by 0 if all images are the same
                finsig[finsig == 0.] = np.inf

            finished = np.zeros_like(npts, dtype=bool)
            finished[(npts < minclip) | (npts <= minkeep)] = True
            while not finished.all():
                if mclip:
                    meds = np.nanmedian(data, axis=-1)
                else:
                    totals = np.nansum(data, axis=-1)
                    # if there are 3+ data points, use the min/max removed avg
                    adjnpts = npts * 1
                    # this adjustment method will fail if there are infinities
                    totals[adjnpts >= 3] -= np.nanmin(data[adjnpts >= 3],
                                                      axis=-1)
                    totals[adjnpts >= 3] -= np.nanmax(data[adjnpts >= 3],
                                                      axis=-1)
                    adjnpts[adjnpts >= 3] -= 2
                    meds = totals / adjnpts

                if reject == 'sigclip':
                    if doscale1:
                        # Compute the sigma with scaling correction.
                        rr = (meds[..., np.newaxis] + zeros) / scales
                        rr[rr < 1.] = 1.
                        rr = np.sqrt(rr)
                        ss = ((data - meds[..., np.newaxis]) / rr) ** 2.
                        # this is just to compute the standard deviation in
                        # each pixel (with 1 degree of freedom)
                        ss = np.sqrt(np.nansum(ss, axis=-1) / (npts - 1))
                        # ignore points where ss == 0 to avoid divide by 0
                        # and because we can't cut them out anyway
                        ss[ss == 0.] = np.inf
                        negresids = ((meds[..., np.newaxis] - data) /
                                     (ss[..., np.newaxis] * rr))
                        bad = ((negresids >= lsigma) |
                               (-1. * negresids >= hsigma))
                    else:
                        # Compute the sigma without scaling correction.
                        ss = (data - meds[..., np.newaxis]) ** 2.
                        # this is just to compute the standard deviation in
                        # each pixel (with 1 degree of freedom)
                        ss = np.sqrt(np.nansum(ss, axis=-1) / (npts - 1))
                        # ignore points where ss == 0 to avoid divide by 0
                        # and because we can't cut them out anyway
                        ss[ss == 0.] = np.inf
                        negresids = ((meds[..., np.newaxis] - data) /
                                     ss[..., np.newaxis])
                        bad = ((negresids >= lsigma) |
                               (-1. * negresids >= hsigma))
                elif reject == 'avsigclip':
                    if doscale1:
                        rr = (meds[..., np.newaxis] + zeros) / scales
                        rr[rr < 1.] = 1.
                        ss = np.sqrt(rr) * finsig[..., np.newaxis]
                        negresids = (meds[..., np.newaxis] - data) / ss
                    else:
                        rr = np.where(meds > 1., meds, [1.])
                        ss = np.sqrt(rr) * finsig[..., np.newaxis]
                        negresids = ((meds[..., np.newaxis] - data) /
                                     ss[..., np.newaxis])
                    bad = ((negresids > lsigma) |
                           (-1. * negresids > hsigma))
                else:
                    if doscale1:
                        rescale = scales * (meds[..., np.newaxis] + zeros)
                        rr = np.where(rescale > 0., rescale, [0.])
                        ss = np.sqrt(nm[:, 0] + rr / nm[:, 1] +
                                     (rr * nm[:, 2]) ** 2) / scales
                    else:
                        rr = np.where(meds > 0., meds, [0.])
                        ss = np.sqrt(nm[:, 0] + rr[..., np.newaxis] / nm[:, 1] +
                                     (rr[..., np.newaxis] * nm[:, 2]) ** 2)

                    negresids = (meds[..., np.newaxis] - data) / ss
                    # these don't pass the cut
                    bad = ((negresids >= lsigma) |
                           (-1. * negresids >= hsigma))

                # number that are bad in a given pixel
                totbad = bad.sum(axis=-1)
                # the easy case of things to update
                toup = (totbad > 0) & (npts - totbad >= minkeep)
                replace = toup[..., np.newaxis] & bad
                data[replace] = np.nan
                # only remove the worst outliers in these cases.
                # not sure there's a pythonic way to do this, so loop it.
                caution = np.where((totbad > 0) & (npts - totbad < minkeep))
                for ii in np.arange(len(caution[0])):
                    inds = tuple(cc[ii] for cc in caution)
                    # get only this pixel's values in every image
                    maxresid = np.abs(negresids[inds])
                    # how many pixels we're allowed to remove
                    nrem = npts[inds] - minkeep[inds]
                    while nrem > 0:
                        # find all values equal to the current max
                        torem = np.where(np.isclose(maxresid,
                                                    maxresid.max()))[0]
                        # only remove all equal values if doing so doesn't
                        # put us below the limit
                        # NOTE: I believe this is a "bug" in IRAF. See this
                        # implementation around line 860 in iccclip.x. Same
                        # logic is applied in icsclip and icaclip. But
                        # IRAF is only wrong for the mclip version,
                        # not the average clip verison. I think
                        # this is the correct way to handle it consistently.
                        if len(torem) <= nrem:
                            data[inds + (torem,)] = np.nan
                        nrem -= len(torem)
                    # mark this one as finished, even if we're keeping
                    # values above minkeep because there are
                    # multiple equal values
                    finished[inds] = True
                # finished is updated where totbad == 0 or npts == minkeep
                # or npts < minclip
                npts = (~np.isnan(data)).sum(axis=-1)
                finished[(totbad == 0) | (npts == minkeep) |
                         (npts < minclip)] = True

        elif reject == 'minmax':
            # number of good points to use
            npts = (~np.isnan(data)).sum(axis=-1)
            npts = npts.astype(float)
            # number of points to reject both high and low
            totlow = nlow * npts
            tothigh = nhigh * npts
            # if we need to remove lower points
            if totlow.max() >= 1.:
                # save the locations of bad points
                rbad = np.isnan(data)
                # while finding mins, replace them with infs
                # necessary because rows of all nans will raise a ValueError in
                # nanargmin or nanargmax
                data[rbad] = np.inf
                # loop through and replace as many as needed
                while totlow.max() >= 1.:
                    mins = np.nanargmin(data, axis=-1)
                    # only the spots we want to remove the min points
                    torep = np.where(totlow >= 1.)
                    torep += (mins[np.where(totlow >= 1.)],)
                    data[torep] = np.inf
                    rbad[torep] = True
                    totlow -= 1.
                data[rbad] = np.nan

            # if we need to remove higher points
            if tothigh.max() >= 1.:
                # save the locations of bad points
                rbad = np.isnan(data)
                # while finding mins, replace them with negative infs
                # necessary because rows of all nans will raise a ValueError in
                # nanargmin or nanargmax
                data[rbad] = -np.inf
                # loop through and replace as many as needed
                while tothigh.max() >= 1.:
                    maxs = np.nanargmax(data, axis=-1)
                    # only the spots we want to remove the max points
                    torep = np.where(tothigh >= 1.)
                    torep += (maxs[np.where(tothigh >= 1.)],)
                    data[torep] = -np.inf
                    rbad[torep] = True
                    tothigh -= 1.
                data[rbad] = np.nan
        elif reject == 'pclip':
            minclip = 3
            npts = (~np.isnan(data)).sum(axis=-1)

            minkeep = np.where(npts < nkeep, npts, [nkeep])

            # Set sign of pclip parameter
            if pclip < 0.:
                tt = -1.
            else:
                tt = 1.

            meds = np.nanmedian(data, axis=-1)
            # the median value or the rightmost of the 2 median values
            n2 = npts // 2

            # pclip should be an integer already, so can probably get
            # rid of the rounding
            if pclip < 0.:
                # whether or not we have an even number of points
                even = (npts % 2) == 0
                # set the clipping index at n2 - abs(pclip), but using the left
                # edge of the median if even
                # can't be below the first index of course
                n3 = np.where(np.round(n2 - even + pclip) > 0,
                              np.round(n2 - even + pclip), [0])
            else:
                # set the clipping index at n2 + pclip, but can't be
                # beyond the max number of images
                n3 = np.where(npts - 1 < np.round(n2 + pclip), npts - 1,
                              np.round(n2 + pclip))
            n3 = n3.astype(int)

            # the data needs to be sorted. NaNs are at the end, beyond npts,
            # so they are ignored.
            data = np.sort(data, axis=-1)

            # Define sigma for clipping
            oinds = np.meshgrid(*[np.arange(dd) for dd in meds.shape],
                                indexing='ij')
            # to fix PyCharm type hinting being wrong for now.
            # meaningless statement in reality
            oinds = list(oinds)
            oinds.append(n3)
            # tuple for a numpy deprecation. Figure out a cleaner way to do
            # this?
            msigma = tt * (data[tuple(oinds)] - meds)
            # skip over ones where sigma is 0 or there's not enough pixels
            msigma[(msigma == 0.) | (npts < minclip)] = np.inf

            negresids = (meds[..., np.newaxis] - data) / msigma[..., np.newaxis]

            # these don't pass the cut
            bad = ((negresids > lsigma) | (-1. * negresids > hsigma))
            # number that are bad in a given pixel
            totbad = bad.sum(axis=-1)
            # the easy case of things to update
            toup = (totbad > 0) & (npts - totbad >= minkeep)
            replace = toup[..., np.newaxis] & bad
            data[replace] = np.nan
            # only remove the worst outliers in these cases.
            # not sure there's a pythonic way to do this, so loop it.
            caution = np.where((totbad > 0) & (npts - totbad < minkeep))
            for ii in np.arange(len(caution[0])):
                inds = tuple(cc[ii] for cc in caution)
                # get only this pixel's values in every image
                maxresid = np.abs(negresids[inds])
                # how many pixels we're allowed to remove
                nrem = npts[inds] - minkeep[inds]
                while nrem > 0:
                    # find all values equal to the current max
                    torem = np.where(np.isclose(maxresid,
                                                maxresid.max()))[0]
                    # only remove all equal values if doing so doesn't
                    #  put us below the limit
                    if len(torem) <= nrem:
                        data[inds + (torem,)] = np.nan
                    nrem -= len(torem)

        if grow > 0.:
            raise NotImplementedError('grow needs to be implemented')
            # grow is only 1-D in IRAF. along the first? dimension?
            # badpts = np.isnan(data)

        if method == 'average':
            # same shape as data, but all ones except for the NaNs
            fullwts = data * 0. + 1.
            # fill in the actual wts values, keeping the NaNs
            fullwts *= wts
            totsum = np.nansum(data * fullwts, axis=-1)
            totwts = np.nansum(fullwts, axis=-1)
            # compute the average
            avg = np.zeros_like(totsum)
            avg[totwts > 0.] = totsum[totwts > 0.] / totwts[totwts > 0.]
            # fill in empty spots with the blank value
            avg[totwts == 0.] = blank
        else:
            npts = (~np.isnan(data)).sum(axis=-1)
            # compute the median
            avg = np.nanmedian(data, axis=-1)
            # fill in empty spots with the blank value
            avg[npts == 0.] = blank
        # save the final result to the output image
        out[0][0].data = avg.astype(oouttype)

        if out[1] is not None:
            npts = (~np.isnan(data)).sum(axis=-1)
            nbad = np.ones_like(npts) * nimages
            # how many bad pixels there were
            nbad -= npts
            out[1][0].data = nbad

        if out[2] is not None:
            # Compute the sigma image line.
            # The estimated sigma includes a correction for the
            # finite population. Weights are used if desired.
            npts = (~np.isnan(data)).sum(axis=-1)
            # same shape as data, but all ones except for the NaNs
            fullwts = data * 0. + 1.
            # fill in the actual wts values, keeping the NaNs
            fullwts *= wts

            sigcor = np.where(npts > 1, npts / (npts - 1), [1])
            wtsum = (data - avg[..., np.newaxis]) ** 2 * fullwts
            wtsum = np.nansum(wtsum, axis=-1)
            totwts = np.nansum(fullwts, axis=-1)

            osigma = np.zeros_like(wtsum, dtype=sigmatype)
            pos = totwts > 0.
            osigma[pos] = np.sqrt(sigcor[pos] * wtsum[pos] /
                                  totwts[pos])
            osigma[~pos] = blank

            out[2][0].data = osigma

        # this is where the icombiner function ends

        # close the input images
        for ifile in imin:
            ifile.close()
        # close the output images
        for ifile in out:
            if ifile is not None:
                ifile.close()
        if logfd is not None:
            logfd.close()
        # delete the input images if requested
        if delete:
            for ifile in iimages:
                os.remove(ifile)
        # close the mask images
        for pm in pms:
            if pm is not None:
                pm.close()

    # return numpy error reporting to its previous state
    np.seterr(invalid=olderr['invalid'])
    return


def ic_stat(imin, section, offarr, project, nim, masktype,
            dothresh, lower, upper, domode=False):
    # Determine the image section parameters.  This must be in terms of
    # the data image pixel coordinates though the section may be specified
    # in terms of the reference image coordinates.  Limit the number of
    # pixels in each dimension to a maximum.

    ndim = imin[0].data.ndim
    if project:
        ndim -= 1

    data = imin[0].data

    starts = np.zeros(ndim)
    ends = np.array(data.shape)[:ndim]
    starts = starts.astype(int)
    ends = ends.astype(int)

    # XXX: implement this ic_section function
    # it's only active when statsec == 'overlap'
    # seems to set va, vb, dv to be image from va to vb with step size dv
    # in each dimension. sticks with default of va = 0 (1 in IRAF),
    # vb = shape (length) of reference image, and dv = 1.
    # call ic_section (section, Memi[va], Memi[vb], Memi[dv], ndim)
    if len(section) > 0:
        raise NotImplementedError('statsec: overlap not yet implemented.')

    oends = ends * 1
    # adjust based on the offsets, but be wary of the bounds
    # equivalent to v1 and v2 in IRAF in the 1->10 loop
    starts -= offarr[nim, :]
    ends -= offarr[nim, :]

    starts[starts < 0] = 0
    ends[ends > oends] = oends[ends > oends]

    # Accumulate the pixel values within the section.  Masked pixels and
    # thresholded pixels are ignored.

    if masktype != 'none':
        # need to carry in the pixel masks to this function,
        # then check that pixels aren't included in the masks before
        # adding them to the stats
        raise NotImplementedError(
            "need to implement bad pixel masks in ic_stat")

    # NOTE: this may not work if project is True
    # get the subgrid of points we want
    idx = tuple(slice(x, y, 1) for x, y in zip(starts, ends))
    data = data[idx]

    # if we're dealing with masks, now would be the time to remove
    # points with bad pixel masks

    # now flatten down the data we want
    data = data.flatten()

    # remove points outside the threshold
    if dothresh:
        data = data[(data >= lower) & (data <= upper)]

    if data.size < 1:
        # this is going to break for things other than fits files
        raise CCDProcError(
            f"Image section contains no pixels: {imin.filename()}")

    mean = data.mean()
    median = np.median(data)
    # weird IRAF mode calculation is brute force slow even on modern
    # computers, so only do it if necessary
    if domode:
        mode = ic_mode(data)
    else:
        mode = None

    return mean, median, mode


def ic_mode(data, zrange=0.8, zstep=0.01, zbin=0.1, nmin=10, maxsize=10000):
    """
    Compute mode of an array.  The mode is found by binning
    with a bin size based on the data range over a fraction of the
    pixels about the median and a bin step which may be smaller than the
    bin size.  If there are too few points the median is returned.

    Parameters
    ----------
    data : ndarray
        Input array to calculate the mode
    zrange : float
        Fraction of pixels about median to use
    zstep : float
        Step size for search for mode
    zbin : float
        Bin size for mode.
    nmin : int
        Minimum number of pixels for mode calculation
    maxsize : int
        This function is crazy slow for limited use, so only use up to
        maxsize random points in this calculation. IRAF use 10k as the
        maxsize for the entire ic_stat function instead of just silly
        mode calculations.
    Returns
    -------
    float : mode of the input array
    """
    if data.size > maxsize:
        data = np.random.choice(data, size=maxsize, replace=False)

    data.sort()

    nn = data.size
    if nn < nmin:
        return np.median(data)

    ii = int(np.floor(nn * (1. - zrange) / 2.))
    jj = int(np.ceil(nn * (1. + zrange) / 2.))

    z1 = data[ii]
    z2 = data[jj]
    if np.isclose(z1, z2):
        return z1

    zstep = zstep * (z2 - z1)
    zbin = zbin * (z2 - z1)

    z1 = z1 - zstep
    kk = ii * 1
    nmax = 0
    mode = np.median(data)
    while kk < jj:
        z1 += zstep
        z2 = z1 + zbin
        while ii < jj and data[ii] < z1:
            ii += 1
        while kk < jj and data[kk] < z2:
            kk += 1
        if kk - ii > nmax:
            nmax = kk - ii
            mode = data[(ii + kk) // 2]

    return mode


def ic_gscale(param, dic, inp, exptime, values, nimages, instrument, project):
    if param is None:
        stype = 'none'
    elif param[0] == '@':
        stype = 'file'
        tmp = np.loadtxt(param[1:])
        if len(tmp.shape) != 1:
            raise ValueError(f"Could not understand values in {param[1:]}")
        if tmp.size < nimages:
            raise ValueError(f"Insufficient values in {param[1:]}")
        if tmp.size > nimages:
            print(f"Warning: Ignoring additional values in {param[1:]}")
        values[:] = tmp[:nimages]
    elif param[0] == '!':
        stype = 'keyword'
        for ii, im in enumerate(inp):
            hv = get_header_value(im, instrument, param[1:])
            if hv is not None:
                values[ii] = hv
                if project:
                    values[:] = hv
                    break
    else:
        param = param.strip().lower()
        if param in dic:
            stype = param
            if stype == 'exposure':
                tmp = np.where(exptime < 0.001)[0]
                values[tmp] = 0.001
        else:
            raise ValueError(f"Unknown scale, zero, or weight type: {param}")
    return stype
