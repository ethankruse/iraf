from .imstat import loadparams, file_handler
from . import instrument, logfile, ssfile
import numpy as np
from astropy.io import fits
import os
import csv
from astropy.wcs import WCS
import sys


def make_fits(path):
    splits = path.split('.')
    if splits[-1] not in ['fits', 'fit']:
        path += '.fits'
    return path

# XXX: how do I "search up the tree" for cl.par mode etc?


def ccdtypes(header):
    # XXX: this needs to have the instrument file header conversions
    options = "object|zero|dark|flat|illum|fringe|other|comp".split('|')
    try:
        typ = header['imagetyp'].lower()
        if typ not in options:
            typ = 'unknown'
    except KeyError:
        typ = 'none'
    return typ


def ccdsubset(im):
    # XXX: this needs to have the instrument file header conversions
    try:
        # XXX: change back to subset
        subsetstr = im[0].header['subset']
        # The default subset identifier is the first word of the subset string.
        subset1 = subsetstr.strip().split()
        if len(subset1) > 0:
            subset1 = subset1[0]
        else:
            subset1 = None
    except KeyError:
        subsetstr = None
        subset1 = None

    # A null subset string is ok.  If not null check for conflict
    # with previous subset IDs.
    if subset1 is not None:
        orig = subset1
        # whether or not to append this to the subsets file
        append = True
        # Search the subset record file for the same subset string.
        # If found use the ID string.  If the subset ID has been
        # used for another subset string then increment an integer
        # suffix to the default ID and check the list again.
        if os.path.exists(ssfile):
            subsetstrs = []
            subsetids = []
            with open(ssfile, 'r') as ff:
                reader = csv.reader(ff, delimiter='\t')
                for row in reader:
                    # skip over blank lines and comment lines
                    if (len(row) == 0 or len(row[0].strip()) == 0 or
                                row[0].strip()[0] == '#'):
                        continue
                    # make sure we have a complete row
                    assert len(row) > 1
                    subsetstrs.append(row[0])
                    subsetids.append(row[1])

            ii = 1
            while True:
                if subsetstr in subsetstrs:
                    append = False
                    subset1 = subsetids[subsetstrs.index(subsetstr)]
                    break
                else:
                    if subset1 in subsetids:
                        subset1 = '{0}{1:d}'.format(orig, ii)
                        ii += 1
                    else:
                        break

        if append:
            with open(ssfile, 'a') as ff:
                writer = csv.writer(ff, delimiter='\t')
                writer.writerow([subsetstr, subset1])

        # Set the subset ID string and replace magic characters by '_'
        # since the subset ID is used in forming image names.
        for ii, ichar in enumerate(subset1):
            if not ichar.isalnum() and ichar != '.':
                subset1[ii] = '_'

    return subset1


def ic_setout(inputs, output, nimages, project, offsets):
    indim = len(inputs[0][0].data.shape)
    outdim = len(output[0][0].data.shape)

    if project:
        outdim = indim - 1
        # also do? IM_NDIM(out[1]) = outdim
    else:
        for im in inputs:
            if len(im[0].data.shape) != outdim:
                print "Image dimensions are not the same"
                sys.exit(1)

    """
    # Set the reference point to that of the first image.
    mw = mw_openim (in[1])
    mwdim = mw_stati (mw, MW_NPHYSDIM)
    call mw_gwtermd (mw, Memd[lref], Memd[wref], Memd[cd], mwdim)
    ct = mw_sctran (mw, "world", "logical", 0)
    call mw_ctrand (ct, Memd[wref], Memd[lref], mwdim)
    call mw_ctfree (ct)
    if (project)
        Memd[lref+outdim] = 1
    """

    # Parse the user offset string.  If "none" then there are no offsets.
    # If "wcs" then set the offsets based on the image WCS.
    # If "grid" then set the offsets based on the input grid parameters.
    # If a file scan it.
    if offsets is None or offsets.lower() == 'none':
        offsetsarr = np.zeros((nimages, outdim), dtype=int)
        reloff = True
    # XXX: implement these
    elif offsets.lower() == 'wcs':
        print 'WCS offsets not implemented yet.'
        sys.exit(1)
    elif offsets.lower() == 'grid':
        print 'grid offsets not implemented yet.'
        sys.exit(1)
    else:
        print 'Manual file offsets not implemented yet.'
        sys.exit(1)

    aligned = True

    for jj in np.arange(outdim):
        aa = offsetsarr[0, jj]
        bb = inputs[0][0].data.shape[jj] + aa
        amin = aa
        bmax = bb
        for ii in np.arange(nimages-1) + 1:
            aa = offsetsarr[ii, jj]
            bb = inputs[ii][0].data.shape[jj] + aa
            if aa != amin or bb != bmax or not reloff:
                aligned = False
            amin = min(aa, amin)
            bmax = max(bb, bmax)
        # also do? IM_LEN(out[1],j) = bmax
        if reloff or amin < 0:
            offsetsarr[:, jj] -= amin

        # also d? IM_LEN(out[1],j) = IM_LEN(out[1],j) - amin

    # Update the WCS.
    # XXX: do this
    if project or not aligned or not reloff:
        print "WCS updates to output files not implemented yet!"

    return offsetsarr, aligned


# this is meant to be equivalent to immap (output, NEW_COPY, Memi[in]) in IRAF
def file_new_copy(outstr, in_header, mode='NEW_COPY', clobber=True):
    if mode != 'NEW_COPY':
        print 'other modes of file_new_copy not supported.'
        return

    # XXX: check this part. Need to add lines to history files, etc.
    in_header.writeto(outstr, clobber=clobber)
    return


def ic_mopen(in_images, out_images, nimages, mtype, mvalue):
    # MASKTYPES	"|none|goodvalue|badvalue|goodbits|badbits|"
    npix = out_images[0][0].data.shape[0]
    """
    # pointer to pms and bufs for each input image
    # pms are the pointers to each image's pixel masks
    call calloc (pms, nimages, TY_POINTER)
    # bufs is just an array of 1s of length npix?
    call calloc (bufs, nimages, TY_POINTER)
    # for every input image, create an array of ones of length npix
    do i = 1, nimages {
        call malloc (Memi[bufs+i-1], npix, TY_INT)
        call amovki (1, Memi[Memi[bufs+i-1]], npix)
    }
    """

    # Check for special cases.  The BOOLEAN type is used when only
    # zero and nonzero are significant; i.e. the actual mask values are
    # not important.  The invert flag is used to indicate that
    # empty masks are all bad rather the all good.
    if mtype.lower() == 'badbits' and mvalue == 0:
        mtype = 'none'
    if mvalue == 0 and mtype.lower() in ['goodvalue', 'goodbits']:
        mtype = 'boolean'
    if ((mvalue == 0 and mtype.lower() in ['badvalue', 'goodbits']) or
            (mvalue != 0 and mtype.lower() == 'goodvalue')):
        invert = True
    else:
        invert = False

    # If mask images are to be used, get the mask name from the image
    # header and open it saving the descriptor in the pms array.
    # Empty masks (all good) are treated as if there was no mask image.
    npms = 0
    if mtype.lower() != 'none':
        for im in in_images:
            try:
                fname = im[0].header['bpm']
                # XXX: implement this
                print 'pixel maps not yet implemented.'
                sys.exit(1)
            except KeyError:
                continue

    # If no mask images are found and the mask parameters imply that
    # good values are 0 then use the special case of no masks.
    if npms == 0 and not invert:
        mtype = 'none'

    """
    # Set up mask structure.
    call calloc (icm, ICM_LEN, TY_STRUCT)
    ICM_TYPE(icm) = mtype
    ICM_VALUE(icm) = mvalue
    ICM_BUFS(icm) = bufs
    ICM_PMS(icm) = pms

    """
    return mtype


def type_max(type1, type2):
    right = np.can_cast(type1, type2, casting='safe')
    left = np.can_cast(type2, type1, casting='safe')

    if left:
        return type1
    if right:
        return type2

    # likely case of an unsigned int and signed int of same size
    ints = [np.int8, np.int16, np.int32, np.int64]
    if (np.issubdtype(type1.type, np.unsignedinteger) and
            np.issubdtype(type2.type, np.integer)):
        for iint in ints:
            if np.can_cast(type1, iint, casting='safe'):
                return np.dtype(iint)

    elif (np.issubdtype(type2.type, np.unsignedinteger) and
              np.issubdtype(type1.type, np.integer)):
        for iint in ints:
            if np.can_cast(type2, iint, casting='safe'):
                return np.dtype(iint)

    print "Unrecognized dtype or cannot safely cast between {0} and {1}.".format(type1, type2)
    sys.exit(1)


def combine(*args, **kwargs):
    params = loadparams(*args, **kwargs)

    inputs = file_handler(params['input'].value)

    if len(inputs) == 0:
        return

    if instrument is not None:
        print "Instrument translation files not yet supported."
        # XXX: need to implement this part

    # Determine whether to divide images into subsets and append extensions.
    dosubsets = params['subsets'].value

    # Go through the input list and eliminate images not satisfying the
    # CCD image type.  Separate into subsets if desired.  Create image
    # and subset lists.

    ccdtypestr = params['ccdtype'].value
    """
    pointer	images		# Pointer to lists of subsets (allocated)
    pointer	extns		# Image extensions for each subset (allocated)
    pointer	subsets		# Subset names (allocated)
    pointer	nimages		# Number of images in subset (allocated)
    int	nsubsets	# Number of subsets
    """
    # lists of images in each subset
    images = []
    # subset names
    subset = []

    for image in inputs:
        # open the image
        try:
            hdulist = fits.open(image)
        except IOError:
            print "Error reading image {0} ...".format(image)
            continue

        thistype = ccdtypes(hdulist[0].header)

        if ccdtypestr is not None and thistype != ccdtypestr:
            hdulist.close()
            continue

        if dosubsets:
            subsetstr = ccdsubset(hdulist)
        else:
            subsetstr = None

        if subsetstr not in subset:
            subset.append(subsetstr)
            images.append([image])
        else:
            images[subset.index(subsetstr)].append(image)

        hdulist.close()

    if len(images) == 0:
        print "No images to combine."
        return

    # Get task parameters.  Some additional parameters are obtained later.
    outroot = params['output'].value
    if len(outroot) == 0:
        print "Must give an output base name"
        return
    plroot = params['plfile'].value
    sigroot = params['sigma'].value

    project = params['project'].value
    combine = params['combine'].value
    reject = params['reject'].value
    blank = params['blank'].value
    gain = params['gain'].value
    rdnoise = params['rdnoise'].value
    snoise = params['snoise'].value
    lthresh = params['lthreshold'].value
    hthresh = params['hthreshold'].value
    lsigma = params['lsigma'].value
    hsigma = params['hsigma'].value
    offsets = params['offsets'].value
    masktype = params['masktype'].value
    maskvalue = params['maskvalue'].value
    scale = params['scale'].value
    zero = params['zero'].value
    weight = params['weight'].value
    statsec = params['statsec'].value

    grow = params['grow'].value
    mclip = params['mclip'].value
    sigscale = params['sigscale'].value
    delete = params['delete'].value
    nkeep = params['nkeep'].value
    pclip = params['pclip'].value
    nlow = params['nlow'].value
    nhigh = params['nhigh'].value
    otype = params['outtype'].value
    outtype = None
    # "short|ushort|integer|long|real|double"
    if otype.lower() == 'short':
        outtype = np.short
    elif otype.lower() == 'ushort':
        outtype = np.ushort
    elif otype.lower() == 'integer':
        outtype = np.intc
    elif otype.lower() == 'long':
        outtype = np.long
    elif otype.lower() == 'real':
        outtype = np.float
    elif otype.lower() == 'double':
        outtype = np.double

    # Check parameters, map INDEFs, and set threshold flag
    if blank is None:
        blank = 0.
    if lsigma is None:
        # XXX: typo? This should probably be -np.inf
        lsigma = np.inf
    if hsigma is None:
        hsigma = np.inf
    if grow is None:
        grow = 0
    if sigscale is None:
        sigscale = 0.

    if lthresh is None and hthresh is None:
        dothresh = False
    else:
        if lthresh is None:
            lthresh = -np.inf
        if hthresh is None:
            hthresh = np.inf

    # Combine each input subset.
    for zz, iset in enumerate(subset):
        iimages = images[zz]

        sstring = subset[zz]
        if sstring is None:
            sstring = ''

        output = '{0}{1}'.format(outroot, sstring)
        output = make_fits(output)

        plfile = None
        sigma = None

        if plroot is not None:
            plfile = '{0}{1}'.format(plroot, sstring)

        if sigroot is not None:
            sigma = '{0}{1}'.format(sigroot, sstring)

        """
            # Combine all images from the (subset) list.
            iferr (call icombine (Memc[Memi[images+i-1]], Memi[nimages+i-1],
            Memc[output], Memc[plfile], Memc[sigma],
            Memc[logfile], NO, delete))
        """
        # Set number of images to combine.
        if project:
            if len(iimages) > 1:
                print "Cannot project combine a list of images"
                return
            hdulist = fits.open(iimages[0])
            shp = hdulist[0].data.shape
            hdulist.close()
            if len(shp) == 1 or shp[-1] == 1:
                print "Can't project one dimensional images"
                return
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

        if reject.lower() == 'pclip':
            if pclip == 0.:
                print "Pclip parameter may not be zero"
                return
            if pclip is None:
                pclip = -0.5

            ii = nimages / 2.
            if np.abs(pclip) < 1.:
                pclip *= ii
            if pclip < 0.:
                pclip = min(-1, max(-ii, int(pclip)))
            else:
                pclip = max(1, min(ii, int(pclip)))

        if reject.lower() == 'minmax':
            if nlow is None:
                nlow = 0.
            if nhigh is None:
                nhigh = 0.

            if nlow >= 1.:
                nlow /= nimages
            if nhigh >= 1.:
                nhigh /= nimages

            ii = nlow * nimages
            jj = nhigh * nimages
            if ii + jj == 0:
                reject = 'none'
            elif ii + jj >= nimages:
                print "Bad minmax rejection parameters"
                return

        out = []
        imin = []
        if project:
            tmp = fits.open(images[0])
            imin.append(tmp)
        else:
            for im in iimages:
                tmp = fits.open(im)
                imin.append(tmp)

        # Map the output image and set dimensions and offsets.
        file_new_copy(output, imin[0], mode='NEW_COPY', clobber=True)
        out.append(fits.open(output))

        offarr, aligned = ic_setout(imin, out, nimages, project, offsets)

        # Determine the highest precedence datatype and set output datatype.
        intype = imin[0][0].data.dtype
        for im in imin:
            intype = type_max(intype, im[0].data.dtype)

        if outtype is None:
            outtype = intype
        # set this? IM_PIXTYPE(out[1]) = getdatatype (clgetc ("outtype"))

        # XXX: this won't work if we're introducing 'sections' of files too
        # need to have a 'return the file without the sections' function
        # e.g. imgimage in IRAF
        # Open pixel list file if given.
        if plfile is not None:
            # make sure it is a .pl file
            base, tail = os.path.split(plfile)
            if len(tail) > 3:
                if len(tail.split('.')) > 1:
                    tail = '.'.join(tail.split('.')[:-1])
                tail += '.pl'
            plfile = os.path.join(base, tail)
            file_new_copy(plfile, out[0], mode='NEW_COPY', clobber=True)
            out.append(fits.open(plfile))
        else:
            out.append(None)

        # Open the sigma image if given.
        if sigma is not None:
            file_new_copy(sigma, out[0], mode='NEW_COPY', clobber=True)
            out.append(fits.open(sigma))
            # XXX: add this?
            # IM_PIXTYPE(out[3]) = ty_max (TY_REAL, IM_PIXTYPE(out[1]))
            # call sprintf (IM_TITLE(out[3]), SZ_IMTITLE,
            # "Combine sigma images for %s")
            # call pargstr (output)
        else:
            out.append(None)

        # XXX: this currently is useless except to make sure masktype == 'none'
        # Open masks.
        masktype = ic_mopen(imin, out, nimages, masktype, maskvalue)

        # Open the log file.
        if logfile is not None:
            logfd = open(logfile, 'a')

        # Memi[in], out, Memi[offsets], nimages
        # icombine(imin, out, offsets, nimages)
        # XXX: this is where the icombine function starts

        # icombine seems to just be a bunch of memory handling to see if the
        # computer will run out of memory. We're skipping this and assuming
        # memory management will be taken care of behind the scenes

        """
        # If aligned use the IMIO buffer otherwise we need vectors of
        # output length.

        if (!aligned) {
            call salloc (dbuf, nimages, TY_POINTER)
            do i = 1, nimages
            call salloc (Memi[dbuf+i-1], npts, TY_REAL)
        }
        """

        # thus, we skip ahead and are now in the ic_combine function
        # call ic_scale (in, out, offsets, scales, zeros, wts, nimages)

        # Set the defaults.
        ncombine = np.ones(nimages)
        exptime = np.zeros(nimages)
        means = np.zeros(nimages)
        medians = np.zeros(nimages)
        modes = np.zeros(nimages)
        means.fill(np.nan)
        medians.fill(np.nan)
        modes.fill(np.nan)
        scales = np.ones(nimages)
        zeros = np.zeros(nimages)
        wts = np.ones(nimages)

        # Get the number of images previously combined and the exposure times.
        # The default combine number is 1 and the default exposure is 0.
        for ii, im in enumerate(imin):
            try:
                ncombine[ii] = im[0].header['ncombine']
            except KeyError:
                pass
            try:
                exptime[ii] = im[0].header['exptime']
            except KeyError:
                pass

        # Set scaling factors.
        stypes = "none|mode|median|mean|exposure".split('|')
        ztypes = "none|mode|median|mean".split('|')
        wtypes = "none|mode|median|mean|exposure".split('|')

        stype = ic_gscale(params['scale'], stypes, imin, exptime, scales, nimages)
        ztype = ic_gscale(params['zero'], ztypes, imin, exptime, zeros, nimages)
        wtype = ic_gscale(params['weight'], wtypes, imin, exptime, wts, nimages)

        # Get image statistics only if needed.
        domode = 'mode' in [stype, ztype, wtype]
        domedian = 'median' in [stype, ztype, wtype]
        domean = 'mean' in [stype, ztype, wtype]

        if domode or domedian or domean:
            # statsec options: "|input|output|overlap|"
            section = None
            oref = None
            if statsec == 'input':
                pass
            elif statsec == 'output':
                oref = out[0]
            elif statsec == 'overlap':
                section = '['
                for ii in np.arange(out[0][0].data.ndim):
                    kk = offarr[0,ii]
                    ll = offarr[0,ii] + imin[0][0].data.shape[ii]
                    for jj in np.arange(1, nimages):
                        kk = max(kk, offarr[jj, ii])
                        ll = min(ll, offarr[jj, ii] + imin[jj][0].data.shape[ii])
                    section += '{0:d}:{1:d},'.format(kk, ll)
                section = section[:-1]
                section += ']'
                oref = out[0]
            else:
                pass

            for ii in np.arange(nimages):
                if oref is None:
                    imref = imin[ii]
                else:
                    imref = oref

                # ic_stat(imin[ii], imref, section, offarr)
        """
        do i = 1, nimages {
        if (imref != out[1])
            imref = in[i]
        call ic_statr (in[i], imref, Memc[section], offsets,
            i, nimages, domode, domedian, domean, mode, median, mean)
        if (domode) {
            Memr[modes+i-1] = mode
            if (stype == S_MODE)
            scales[i] = mode
            if (ztype == S_MODE)
            zeros[i] = mode
            if (wtype == S_MODE)
            wts[i] = mode
        }
        if (domedian) {
            Memr[medians+i-1] = median
            if (stype == S_MEDIAN)
            scales[i] = median
            if (ztype == S_MEDIAN)
            zeros[i] = median
            if (wtype == S_MEDIAN)
            wts[i] = median
        }
        if (domean) {
            Memr[means+i-1] = mean
            if (stype == S_MEAN)
            scales[i] = mean
            if (ztype == S_MEAN)
            zeros[i] = mean
            if (wtype == S_MEAN)
            wts[i] = mean
        }
        }
        """

        # XXX: this is where the icombine function ends

        # close the input images
        for ifile in imin:
            ifile.close()
        # close the output images
        for ifile in out:
            if ifile is not None:
                ifile.close()
        logfd.close()

    return params


def ic_stat(imin, imref, section, offarr):
    pass


def ic_gscale(param, dic, inp, exptime, values, nimages):
    if param.value is None:
        stype = 'none'
    elif param.value[0] == '@':
        stype = 'file'
        tmp = np.loadtxt(param.value[1:])
        if len(tmp.shape) != 1:
            print "Could not understand {0} values in {1}".format(param.name, param.value[1:])
            sys.exit(1)
        if tmp.size < nimages:
            print "Insufficient {0} values in {1}".format(param.name, param.value[1:])
            sys.exit(1)
        if tmp.size > nimages:
            print "Warning: Ignoring additional {0} values in {1}".format(param.name, param.value[1:])
        values[:] = tmp[:nimages]
    elif param.value[0] == '!':
        stype = 'keyword'
        for ii, im in enumerate(inp):
            values[ii] = im[0].header[param.value[1:]]
    else:
        if param.value in dic:
            stype = param.value
            if stype == 'exposure':
                tmp = np.where(exptime > 0.001)[0]
                values[tmp] = 0.001
        else:
            print "Unknown {0} type".format(param.name)
            sys.exit(1)
    return stype


def icombine(ins, out, offsets, nimages):
    return
