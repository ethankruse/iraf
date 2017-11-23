from iraf.utils import file_handler
import numpy as np
import os
import sys
from iraf.sys import image_open, image_close
import re


__all__ = ['combine', 'Instrument', 'ccdtypes', 'get_header_value']

# XXX: deal with hdulist assuming we have a fits file. trace all image_open and
# see if what I do is fits specific.


class Instrument(object):
    """
    Instrument object.
    """
    """
    To learn about instruments, look into
    'iraf-src/noao/imred/ccdred/src/hdrmap.x' for the functions that create
    a symbol table pointer and do the translations.
    Instrument file tables are located in
    'iraf-src/noao/imred/ccdred/ccddb'.

    values in fits headers/etc actually used so far:
    imagetyp 
    and its header values: object|zero|dark|flat|illum|fringe|other|comp
    and also none, unknown.
    subset
    
    """
    # XXX: need to handle allowing 'default' values for things
    def __init__(self, name=None):
        self.definitions = {}

        if name is None or name.strip().lower() == 'default':
            # insert how to interpret things
            # XXX: filter or filters?
            self.definitions['subset'] = 'filter'
            self.definitions['bias'] = 'zero'
            self.definitions['dome flat'] = 'flat'
            self.definitions['projector flat'] = 'flat'
            self.definitions['comparison'] = 'comp'
            self.definitions['sky flat'] = 'object'

        else:
            print(f'Instrument {name} not implemented yet.')
            sys.exit(1)

    def translate(self, key):
        # if it's a string, strip the value and put it in lowercase.
        try:
            key = key.strip().lower()
        except AttributeError:
            pass

        if key in self.definitions:
            return self.definitions[key]
        else:
            return key

    """
    Set up an Instrument class with a translate() function or something.
    Maybe more like to_default() or from_default() 
    """


def make_fits(path):
    splits = path.split('.')
    if splits[-1] not in ['fits', 'fit']:
        path += '.fits'
    return path


def set_header_value(hdulist, instrument, key, value):
    """
    The equivalent of IRAF's hdmput*. (e.g. hdmputi, hdmputr)

    Parameters
    ----------
    hdulist
    instrument
    key
    value

    Returns
    -------

    """
    # what is the header key in the instrument's language
    key = instrument.translate(key)

    # go through the list of hdus
    found = False
    for hdu in hdulist:
        if key in hdu.header:
            hdu.header.set(key, value)
            found = True
            break
    # if it's not in any headers, put it in the first one
    if not found:
        hdulist[0].header.set(key, value)
    return


def get_header_value(hdulist, instrument, key):
    """
    The equivalent of IRAF's hdmg*. (e.g. hdmgstr)

    Parameters
    ----------
    hdulist
    instrument
    key

    Returns
    -------

    """
    # what is the header key in the instrument's language
    key = instrument.translate(key)
    val = None

    # go reverse order because we'll typically want the first instance
    # if something is in multiple headers
    for hdu in hdulist[::-1]:
        try:
            # get the instrument's value in the header
            val = hdu.header[key]
            # translate that back to our normalized language
            val = instrument.translate(val)
        except KeyError:
            pass

    if val is None:
        # XXX: we need a way to use a default value here if given in the
        # instrument object
        pass

    return val


def ccdtypes(hdulist, instrument):
    """
    Get the header value of 'imagetyp' (or instrument equivalent).
    If that (instrument converted) value is one of
    "object|zero|dark|flat|illum|fringe|other|comp" return that.
    Otherwise return 'unknown' unless the header does not contain any
    'imagetyp' and the instrument doesn't have a default value for it. Then
    return 'none'.

    Parameters
    ----------
    hdulist
    instrument : Instrument

    Returns
    -------

    """

    if not isinstance(instrument, Instrument):
        print('ccdtypes not given an Instrument object.')
        sys.exit(1)

    options = "object|zero|dark|flat|illum|fringe|other|comp".split('|')
    typ = get_header_value(hdulist, instrument, 'imagetyp')

    if typ is None:
        typ = 'none'
    elif typ not in options:
        typ = 'unknown'

    return typ


def ccdsubset(hdulist, instrument, ssfile):

    if not isinstance(instrument, Instrument):
        print('ccdsubset not given an Instrument object.')
        sys.exit(1)

    subsetstr = get_header_value(hdulist, instrument, 'subset')

    if subsetstr is None:
        subsetstr = ''

    subsetstr = subsetstr.strip()

    # Replace non alphanumeric or '.' characters by '_'
    # since the subset ID is used in forming image names.
    subsetstr = re.sub(r'[^\w.]', '_', subsetstr)

    """
    # This bit was a translation of the original IRAF ccdsubset function
    # that uses the subsets file and turns things into shorter subset strings
    # while also making sure there aren't overlaps. This all seems unnecessary
    # and we should just use the full subset string as the identifier.

    import shlex
    
    if ssfile is None:
        print('ssfile must be defined to use subsets')
        sys.exit(1)
        
    # The default subset identifier is the first
    # word/string of the subset string.
    subset1 = shlex.split(subsetstr.strip())
    if len(subset1) > 0:
        subset1 = subset1[0]
    else:
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
                lines = ff.readlines()
                for row in lines:
                    # skip over blank lines and comment lines
                    if (len(row) == 0 or len(row.strip()) == 0 or
                            row[0].strip()[0] == '#'):
                        continue

                    groups = shlex.strip(row)
                    assert len(groups) == 2
                    subsetstrs.append(groups[0])
                    subsetids.append(groups[1])

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
                ff.write('{0}\t{1}\n'.format(subsetstr, subset1))

        # Set the subset ID string and replace magic characters by '_'
        # since the subset ID is used in forming image names.
        for ii, ichar in enumerate(subset1):
            if not ichar.isalnum() and ichar != '.':
                subset1[ii] = '_'

    return subset1
    """

    return subsetstr


def ic_setout(inputs, output, nimages, project, offsets):

    indim = len(inputs[0][0].data.shape)
    outdim = len(output[0][0].data.shape)

    if project:
        outdim = indim - 1
        # XXX: IM_NDIM(out[1]) = outdim
    else:
        for im in inputs:
            if len(im[0].data.shape) != outdim:
                print("Image dimensions are not the same")
                sys.exit(1)

    """
    # Set the reference point to that of the first image.
    # Open an MWCS descriptor on an image
    mw = mw_openim (in[1])
    # Get the value of a MWCS interface parameter.
    mwdim = mw_stati (mw, MW_NPHYSDIM)
    # Get the linear part of the Wterm, i.e., the physical and world
    # coordinates of the reference point and the CD matrix.
    call mw_gwtermd (mw, Memd[lref], Memd[wref], Memd[cd], mwdim)
    # Set up a coordinate transformation (CTRAN) descriptor.
    ct = mw_sctran (mw, "world", "logical", 0)
    # Transform a single N-dimensional point
    call mw_ctrand (ct, Memd[wref], Memd[lref], mwdim)
    call mw_ctfree (ct)
    if (project)
        Memd[lref+outdim] = 1
    """

    # Parse the user offset string.  If "none" then there are no offsets.
    # If "wcs" then set the offsets based on the image WCS.
    # If "grid" then set the offsets based on the input grid parameters.
    # If a file scan it.
    if offsets == 'none':
        offsetsarr = np.zeros((nimages, outdim), dtype=int)
        reloff = True
    # XXX: implement these
    elif offsets.lower() == 'wcs':
        print('WCS offsets not implemented yet.')
        sys.exit(1)
    elif offsets.lower() == 'grid':
        print('grid offsets not implemented yet.')
        sys.exit(1)
    else:
        print('Manual file offsets not implemented yet.')
        sys.exit(1)

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
        print("WCS updates to output files not implemented yet!")
        sys.exit(1)

    return aligned, offsetsarr


# this is meant to be equivalent to immap (output, NEW_COPY, Memi[in]) in IRAF
def file_new_copy(outstr, in_header, mode='NEW_COPY', overwrite=True):
    if mode != 'NEW_COPY':
        print('other modes of file_new_copy not supported.')
        return

    # XXX: check this part. Need to add lines to history files, etc.
    # Also this only works for FITS files.
    if in_header.__filetype__ == 'fits':
        in_header.writeto(outstr, overwrite=overwrite)
        # XXX: are some header parameters added/changed at this stage?
        # see sys/imio/immapz.x

    else:
        err = 'file_new_copy of file type {0} not yet implemented.'
        print(err.format(in_header.__filetype__))
        sys.exit(1)
    return


def ic_mopen(in_images, out_images, nimages, mtype, mvalue, instrument):
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

    mtype = mtype.strip().lower()

    if mtype not in ['none', 'goodvalue', 'badvalue', 'goodbits', 'badbits']:
        print('masktype {0} not recognized. Assuming "none".'.format(mtype))
        mtype = 'none'

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
    if mtype != 'none':
        for im in in_images:
            fname = get_header_value(im, instrument, 'BPM')
            if fname is None:
                continue

            # XXX: implement this
            print('pixel maps not yet implemented.')
            sys.exit(1)

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

    """
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
    """
    errstr = "Unrecognized dtype or cannot safely cast between {0} and {1}."
    print(errstr.format(type1, type2))
    sys.exit(1)


def combine(images, output, *, plfile=None, sigma=None, ccdtype=None,
            subsets=False, delete=False, method='average',
            reject='none', project=False, outtype=None, offsets='none',
            masktype='none', maskvalue=0., blank=0., scale=None, zero=None,
            weight=None, statsec=None, lthreshold=None, hthreshold=None,
            nlow=1, nhigh=1, nkeep=1, mclip=True, lsigma=3.0,
            hsigma=3.0, rdnoise=0., gain=1., snoise=0., sigscale=0.1,
            pclip=-0.5, grow=0, instrument=None, logfile=None, verbose=False,
            ssfile=None):
    """

    Parameters
    ----------
    images :
        List of images to combine
    output :
        List of output images
    plfile :
        List of output pixel list files (optional)
    sigma :
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
    ssfile

    Returns
    -------

    """
    offsets = offsets.strip().lower()

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

    method = method.strip().lower()
    if method not in ['average', 'median']:
        print('Combine method not recognized: {0}'.format(method))
        sys.exit(1)

    for image in inputs:
        # open the image
        hdulist = image_open(image)

        if hdulist is None:
            continue

        # what image type is this one
        thistype = ccdtypes(hdulist, instrument)

        # if this isn't the image type we're looking for, skip it
        if ccdtype is not None and thistype != ccdtype:
            image_close(hdulist)
            continue

        if subsets:
            subsetstr = ccdsubset(hdulist, instrument, ssfile)
            # As far as I can tell, the subset and extn list is the same.
            # extn = subsetstr
        else:
            subsetstr = ''

        if subsetstr not in subset:
            subset.append(subsetstr)
            images.append([image])
        else:
            images[subset.index(subsetstr)].append(image)

        image_close(hdulist)

    # end of cmb_images code.

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
    if sigma is not None:
        sigroot = sigma.strip()
    else:
        sigroot = None

    # XXX: make sure we're not changing input values in each loop.
    # Combine each input subset.
    for zz, iset in enumerate(subset):

        iimages = images[zz]

        # this apparently isn't done in IRAF, but I like separating the two
        # parts with a '.' if they both exist
        if len(iset) > 0:
            comb = '.'
        else:
            comb = ''

        base, ext = os.path.splitext(outroot)
        output = '{0}{1}{2}{3}'.format(base, comb, iset, ext)

        if plroot is not None:
            base, ext = os.path.splitext(plroot)
            plfile = '{0}{1}{2}{3}'.format(base, comb, iset, ext)
        else:
            plfile = None

        if sigroot is not None:
            base, ext = os.path.splitext(sigroot)
            sigma = '{0}{1}{2}{3}'.format(base, comb, iset, ext)
        else:
            sigma = None

        # icombine starts here
        """
        # Combine all images from the (subset) list.
        call icombine (Memc[Memi[images+i-1]], Memi[nimages+i-1],
        Memc[output], Memc[plfile], Memc[sigma],
        Memc[logfile], NO, delete))
        """
        # Set number of images to combine.
        if project:
            if len(iimages) > 1:
                print("Cannot project combine a list of images")
                return
            hdulist = image_open(iimages[0])
            shp = np.array(hdulist[0].data.shape)
            image_close(hdulist)
            if shp.size == 1 or shp[shp > 1].size == 1:
                print("Can't project one dimensional images")
                return
            # XXX: which dimension do we project over?
            nimages = shp[-1]
        else:
            nimages = len(iimages)

        # Convert the nkeep parameter if needed.
        # XXX: do we want to do this here or wait until later?
        if nkeep < 0:
            nkeep = max(0, nimages + nkeep)

        # Convert the pclip parameter to a number of pixels rather than
        # a fraction.  This number stays constant even if pixels are
        # rejected.  The number of low and high pixel rejected, however,
        # are converted to a fraction of the valid pixels.

        reject = reject.lower().strip()
        rejectopts = "none|ccdclip|crreject|minmax|pclip|sigclip|avsigclip"
        if reject not in rejectopts.split('|'):
            print('Could not recognize reject parameter {0}'.format(reject))
            sys.exit(1)
        # define	REJECT	"|none|ccdclip|crreject|minmax|pclip|sigclip|avsigclip|"
        if reject == 'pclip':
            if pclip == 0.:
                print("Pclip parameter can not be zero when reject=='pclip'")
                return

            ii = nimages // 2
            if np.abs(pclip) < 1.:
                pclip *= ii
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
                print("Bad minmax rejection parameters")
                return

        imin = []
        # Map the input image(s).
        for im in iimages:
            tmp = image_open(im)
            imin.append(tmp)

        print(output)

        out = []

        # Map the output image and set dimensions and offsets.
        file_new_copy(output, imin[0], mode='NEW_COPY', overwrite=True)
        out.append(image_open(output))

        # start of ic_setout
        aligned, offarr = ic_setout(imin, out, nimages, project, offsets)

        # Determine the highest precedence datatype and set output datatype.
        intype = imin[0][0].data.dtype
        for im in imin:
            intype = type_max(intype, im[0].data.dtype)

        if outtype is None:
            outtype = intype
        # make sure the output type will work given the input
        outtype = type_max(intype, outtype)

        # Open pixel list file if given.
        if plfile is not None:
            # XXX: this won't work if we're introducing 'sections' of files too
            # need to have a 'return the file without the sections' function
            # e.g. imgimage in IRAF

            # make sure it is a .pl file
            base, tail = os.path.split(plfile)
            if len(tail) > 3:
                # remove the suffix if it has one
                if len(tail.split('.')) > 1:
                    tail = '.'.join(tail.split('.')[:-1])
            tail += '.pl'
            plfile = os.path.join(base, tail)
            file_new_copy(plfile, out[0], mode='NEW_COPY', overwrite=True)
            out.append(image_open(plfile))
        else:
            out.append(None)

        sigmatype = None
        # Open the sigma image if given.
        if sigma is not None:
            file_new_copy(sigma, out[0], mode='NEW_COPY', overwrite=True)
            out.append(image_open(sigma))
            # has to be a float
            sigmatype = type_max(np.float, outtype)
        else:
            out.append(None)

        # XXX: this currently is useless except to make sure masktype == 'none'
        # Open masks.
        masktype = ic_mopen(imin, out, nimages, masktype, maskvalue, instrument)

        # Open the log file.
        logfd = None
        if logfile is not None:
            logfd = open(logfile, 'a')

        # Memi[in], out, Memi[offsets], nimages
        # icombiner(imin, out, offsets, nimages)
        # this is where the icombiner function starts

        # icombiner seems to just be a bunch of memory handling to see if the
        # computer will run out of memory. We're skipping this and assuming
        # memory management will be taken care of behind the scenes. Move
        # to ic_combiner.

        # beginning of ic_scale
        # call ic_scale (in, out, offsets, scales, zeros, wts, nimages)

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

        if domode or domedian or domean:
            # statsec options: "|input|output|overlap|"
            if statsec is None:
                statsec = ''
            statsec = statsec.strip().lower()

            oimref = None
            section = statsec
            if statsec == 'input':
                section = ''
            elif statsec == 'output':
                section = ''
                oimref = out[0]
            elif statsec == 'overlap':
                section = '['
                for ii in np.arange(out[0][0].data.ndim):
                    kk = offarr[0, ii]
                    ll = offarr[0, ii] + imin[0][0].data.shape[ii]
                    for jj in np.arange(1, nimages):
                        kk = max(kk, offarr[jj, ii])
                        ll = min(ll,
                                 offarr[jj, ii] + imin[jj][0].data.shape[ii])
                    section += '{0:d}:{1:d},'.format(kk, ll)
                section = section[:-1]
                section += ']'
                oimref = out[0]

            for ii in np.arange(nimages):
                mn, md, mode = ic_stat(imin[ii], oimref, section, offarr,
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
            scales = 1./scales
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
            # a nice log set the zero offsets in this case.
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
        for ii in np.arange(nimages-1) + 1:
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
            logtxt += '{0} : IMCOMBINE\n'.format(now)
            logtxt += '  combine = {0}, '.format(method)
            logtxt += 'scale = {0}, zero = {1}, '.format(scale, zero)
            logtxt += 'weight = {0}\n'.format(weight)

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
                    ostr = '  rdnoise = {0}'.format(rdnoise)
                else:
                    ostr = '  rdnoise = {0:g}'.format(rdnoise)
                if isinstance(gain, str):
                    ostr += ', gain = {0}, '.format(gain)
                else:
                    ostr += ', gain = {0:g}, '.format(gain)
                if isinstance(snoise, str):
                    ostr += 'snoise = {0}, '.format(snoise)
                else:
                    ostr += 'snoise = {0:g}, '.format(snoise)
                ostr += 'lsigma = {0:g}, hsigma = {1:g}\n'.format(lsigma,
                                                                  hsigma)
                logtxt += ostr
            elif reject == 'crreject':
                ostr = '  reject = crreject, mclip = {0}, nkeep = {1:d}\n'
                logtxt += ostr.format(mclip, nkeep)
                if isinstance(rdnoise, str):
                    ostr = '  rdnoise = {0}'.format(rdnoise)
                else:
                    ostr = '  rdnoise = {0:g}'.format(rdnoise)
                if isinstance(gain, str):
                    ostr += ', gain = {0}, '.format(gain)
                else:
                    ostr += ', gain = {0:g}, '.format(gain)
                if isinstance(snoise, str):
                    ostr += 'snoise = {0}, '.format(snoise)
                else:
                    ostr += 'snoise = {0:g}, '.format(snoise)
                ostr += 'hsigma = {0:g}\n'.format(hsigma)
                logtxt += ostr
            elif reject == 'pclip':
                ostr = '  reject = pclip, nkeep = {0:d}\n'.format(nkeep)
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
                logtxt += '  grow = {0:d}\n'.format(grow)

            if dothresh:
                if lthreshold > -np.inf and hthreshold < np.inf:
                    ostr = '  lthreshold = {0:g}, hthreshold = {1:g}\n'
                    ostr = ostr.format(lthreshold, hthreshold)
                elif lthreshold > -np.inf:
                    ostr = '  lthreshold = {0:g}\n'.format(lthreshold)
                else:
                    ostr = '  hthreshold = {0:g}\n'.format(hthreshold)
                logtxt += ostr

            logtxt += '  blank = {0:g}\n'.format(blank)
            if len(statsec) > 0:
                logtxt += '  statsec = {0}\n'.format(statsec)

            if masktype != 'none':
                print('Mask types not yet supported')
                sys.exit(1)

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
                    print('Mask types not yet supported')
                    sys.exit(1)
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
                    stc = 'stck{0:04d}'.format(ii)
                    vl = get_header_value(imin[ii], instrument, stc)
                    if vl is not None:
                        ostr = '  {0:21s}'.format(vl)
                    else:
                        ostr = '  {0:16s}[{1:03d}]'
                        ostr = ostr.format(imin[ii].filename(), ii)
                elif project:
                    ostr = '  {0:16s}[{1:03d}]'
                    ostr = ostr.format(imin[ii].filename(), ii)
                else:
                    ostr = '  {0:21s}'.format(imin[ii].filename())
                if dop['ncombine']:
                    ostr += ' {0:6d}'.format(ncombine[ii])
                if dop['exptime']:
                    ostr += ' {0:6.1f}'.format(exptime[ii])
                if dop['mode']:
                    ostr += ' {0:7.5g}'.format(modes[ii])
                if dop['median']:
                    ostr += ' {0:7.5g}'.format(medians[ii])
                if dop['mean']:
                    ostr += ' {0:7.5g}'.format(means[ii])
                if dop['rdn']:
                    rval = get_header_value(imin[ii], instrument, rdnoise)
                    ostr += ' {0:7g}'.format(rval)
                if dop['gain']:
                    rval = get_header_value(imin[ii], instrument, gain)
                    ostr += ' {0:6g}'.format(rval)
                if dop['sn']:
                    rval = get_header_value(imin[ii], instrument, snoise)
                    ostr += ' {0:6g}'.format(rval)
                if doscale:
                    ostr += ' {0:6.3f}'.format(1./scales[ii])
                if dozero:
                    ostr += ' {0:7.5g}'.format(-zeros[ii])
                if dowts:
                    ostr += ' {0:6.3f}'.format(wts[ii])
                if not aligned:
                    nd = np.array(out[0][0].data.shape)
                    # number of dimensions in out array
                    nd = len(nd[nd > 1])
                    if nd == 1:
                        ostr += ' {0:9d}'.format(offarr[ii, 0])
                    else:
                        use = offarr[ii, :]
                        strs = [' {0:4d}'.format(ixx) for ixx in use]
                        for istr in strs:
                            ostr += istr
                if dop['mask']:
                    print('Mask types not yet supported')
                    sys.exit(1)
                ostr += '\n'
                logtxt += ostr

            # Log information about the output images.
            ostr = '\n  Output image = {0}, ncombine = {1:d}\n'
            ostr = ostr.format(out[0].filename(), nout)
            logtxt += ostr

            if out[1] is not None:
                logtxt += '  Pixel list image = {0}\n'.format(out[1].filename())

            if out[2] is not None:
                logtxt += '  Sigma image = {0}\n'.format(out[1].filename())

            if verbose:
                print(logtxt)
            if logfd is not None:
                logfd.write(logtxt)

        doscale = (doscale or dozero)
        # end of the icscale function
        keepids = False
        if method == 'average' and dowts:
            keepids = True
        elif method == 'median':
            dowts = False

        docombine = True

        # Set rejection algorithm specific parameters
        if reject in ['ccdclip', 'crreject']:
            # the column order is readnoise, gain, snoise
            nm = np.zeros((nimages,3))
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
            # adjust the readnoise values? don't let it be actually 0?
            rdsq = (nm[:, 0]/nm[:, 1])**2
            nm[:, 0] = np.where(rdsq > small, rdsq, [small])

            if isinstance(snoise, str):
                for ii in np.arange(nimages):
                    nm[ii, 2] = get_header_value(imin[ii], instrument, snoise)
            else:
                nm[:, 2] = snoise

            if not keepids:
                if doscale1 or grow > 0:
                    keepids = True
                else:
                    # see if every row is the same as the first row
                    for ii in np.arange(nimages-1)+1:
                        if (nm[ii, :] != nm[0, :]).any():
                            keepids = True
                            break
            if reject == 'crreject':
                lsigma = np.inf
        elif reject == 'minmax':
            mclip = False
            if grow > 0:
                keepids = True
        elif reject == 'pclip':
            mclip = True
            if grow > 0:
                keepids = True
        elif reject in ['sigclip', 'avsigclip']:
            if doscale1 or grow > 0:
                keepids = True
        else:
            # case 'none'
            mclip = False
            grow = 0

        # start of ic_gdatar
        oshp = out[0][0].data.shape
        oshp += (nimages,)
        data = np.zeros(oshp)
        # set bad values
        data[:] = np.nan
        # easy case, no offsets to deal with
        if aligned:
            for ii in np.arange(nimages):
                data[..., ii] = imin[ii][0].data
        else:
            print('need to handle offsets in combine')
            sys.exit(1)

        # remove the data points with bad mask values (set to nan)
        if masktype != 'none':
            print('need to implement mask types')
            sys.exit(1)

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
        if reject in ['ccdclip', 'crreject']:
            minclip = 2
            if mclip:
                # XXX: set docombine somewhere if that is actually useful?
                # number of good points to use
                npts = (~np.isnan(data)).sum(axis=-1)
                if nkeep < 0:
                    minkeep = np.where(npts + nkeep > 0, npts + nkeep, [0])
                else:
                    minkeep = np.where(npts < nkeep, npts, [nkeep])

                finished = np.zeros_like(npts, dtype=bool)
                finished[(npts < minclip) | (npts <= minkeep)] = True
                while not finished.all():
                    # XXX: catch RuntimeWarning for all nan rows
                    meds = np.nanmedian(data, axis=-1)

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
                    bad = ((negresids >= lthreshold) |
                           (-1. * negresids >= hthreshold))
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
                    data[torep] = np.nan
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
                    data[torep] = np.nan
                    tothigh -= 1.
                data[rbad] = np.nan
        elif reject == 'pclip':
            pass
        elif reject == 'sigclip':
            pass
        elif reject == 'avsigclip':
            pass


        # XXX: this is where the icombiner function ends

        # close the input images
        for ifile in imin:
            image_close(ifile)
        # close the output images
        for ifile in out:
            if ifile is not None:
                image_close(ifile)
        if logfd is not None:
            logfd.close()

    return


def ic_stat(imin, imref, section, offarr, project, nim, masktype,
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
        print('statsec: overlap not yet implemented.')
        sys.exit(1)

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
        print('need to implement bad pixel masks in ic_stat')
        # need to carry in the pixel masks to this function,
        # then check that pixels aren't included in the masks before
        # adding them to the stats
        sys.exit(1)

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
        # this is going to break for things other than fits
        print('Image section contains no pixels: {0}'.format(imin.filename()))
        sys.exit(1)

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
            mode = data[(ii + kk)//2]

    return mode


def ic_gscale(param, dic, inp, exptime, values, nimages, instrument, project):
    if param is None:
        stype = 'none'
    elif param[0] == '@':
        stype = 'file'
        tmp = np.loadtxt(param[1:])
        if len(tmp.shape) != 1:
            print("Could not understand values in {0}".format(param[1:]))
            sys.exit(1)
        if tmp.size < nimages:
            print("Insufficient values in {0}".format(param[1:]))
            sys.exit(1)
        if tmp.size > nimages:
            print("Warning: Ignoring additional values in {0}".format(param[1:]))
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
            print("Unknown scale, zero, or weight type: {0}".format(param))
            sys.exit(1)
    return stype
