from .imstat import loadparams, file_handler
from . import instrument, logfile, ssfile
import numpy as np
from astropy.io import fits
import os
import csv


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
    print dosubsets
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

    grow = params['grow'].value
    mclip = params['mclip'].value
    sigscale = params['sigscale'].value
    delete = params['delete'].value
    nkeep = params['nkeep'].value
    pclip = params['pclip'].value
    nlow = params['nlow'].value
    nhigh = params['nhigh'].value

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
        if project:
            tmp = fits.open(images[0])
            out.append(tmp)
        else:
            for im in iimages:
                tmp = fits.open(im)
                out.append(tmp)

        # Map the output image and set dimensions and offsets.
        imout = out[0]
        # XXX: check this function. Need to add lines to history files, etc.
        imout.writeto(output)




        # stack1 = stack = NO = False
        """
        # Map the output image and set dimensions and offsets.
        tmp = immap (output, NEW_COPY, Memi[in]); out[1] = tmp
        if (stack1 == YES) {
        call salloc (key, SZ_FNAME, TY_CHAR)
        do i = 1, nimages {
            call sprintf (Memc[key], SZ_FNAME, "stck%04d")
            call pargi (i)
            call imdelf (out[1], Memc[key])
        }
        }
        call salloc (offsets, nimages*IM_NDIM(out[1]), TY_INT)
        call ic_setout (Memi[in], out, Memi[offsets], nimages)
        """

        # close the input images
        for ifile in out:
            ifile.close()

    return params
