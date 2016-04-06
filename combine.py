from .imstat import loadparams, file_handler
from . import instrument, logfile
import numpy as np
from astropy.io import fits


def ccdtypes(header):
    options = "object|zero|dark|flat|illum|fringe|other|comp".split('|')
    try:
        typ = header['imagetyp'].lower()
        if typ not in options:
            typ = 'unknown'
    except KeyError:
        typ = 'none'
    return typ


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

    ccdtype = 0
    ccdtypestr = params['ccdtype'].value

    images = []
    for image in inputs:
        # open the image
        try:
            hdulist = fits.open(image)
        except IOError:
            print "Error reading image {0} ...".format(image)
            continue

        thistype = ccdtypes(hdulist[0].header)
        if ccdtypestr is not None and thistype != ccdtypestr:
            continue

        if dosubsets:
            # XXX: implement this
            print "Subsets not implemented yet."

        images.append(image)

        hdulist.close()

    if len(images) == 0:
        print "No images to combine."
        return

    # Get task parameters.  Some additional parameters are obtained later.
    output = file_handler(params['output'].value)
    plfile = file_handler(params['plfile'].value)
    sigma = file_handler(params['sigma'].value)

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

    if project:
        if len(images) > 1:
            print "Cannot project combine a list of images"
            return
        hdulist = fits.open(images[0])
        shp = hdulist[0].data.shape
        if len(shp) == 1 or shp[-1] == 1:
            print "Can't project one dimensional images"
            return
        nimages = shp[-1]
    else:
        nimages = len(images)

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


    """

    """







    return params
