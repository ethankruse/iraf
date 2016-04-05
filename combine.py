from .imstat import loadparams, file_handler
from . import instrument, logfile
import numpy as np
from astropy.io import fits


def ccdtypes(header):
    options = "object|zero|dark|flat|illum|fringe|other|comp".split('|')
    try:
        typ = header['imagetyp']
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

    """


    """







    return params
