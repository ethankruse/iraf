from iraf.utils import file_handler
import numpy as np
import scipy.stats
from iraf.sys import image_open


__all__ = ['imstatistics']


def imstatistics(images, *, fields="image,npix,mean,stddev,min,max",
                 lower=None, upper=None, nclip=0, lsigma=3.0, usigma=3.0,
                 print_format=True):
    """
    Get general statistics about the data in a file's (or list of files')
    primary FITS HDU.

    Parameters
    ----------
    images : IRAF file list
        Any valid IRAF file descriptors to run statistics for.
    fields : string
        Which statistics to calculate and print. CSV list chosen amongst
        "image,npix,min,max,mean,midpt,mode,stddev,skew,kurtosis".
    lower : float
        Lower limit for pixel values to include in the statistics
    upper : float
        Upper limit for pixel values to include in the statistics.
    nclip : int
        Number of times to run sigma clipping before calculating
        the statistics.
    lsigma : float
        Lower clipping factor (in sigma). e.g. 3.0 would remove any pixels
        with values 3.0 sigma below the mean.
    usigma : float
        Upper clipping factor (in sigma). e.g. 3.0 would remove any pixels
        with values 3.0 sigma above the mean.
    print_format : boolean
        Whether or not to format the output and print column labels
    """

    images = file_handler(images)

    if lower is None:
        lower = -np.inf
    if upper is None:
        upper = np.inf

    possible_fields = "image|npix|min|max|mean|midpt|mode|stddev|skew|kurtosis"
    possible_fields = possible_fields.split('|')

    # user input fields
    in_fields = [x.strip().lower() for x in fields.split(',')]
    # XXX: print warning if a user field isn't in the possible list?
    # retain the same order as in_fields, but only the valid ones
    use_fields = [x for x in in_fields if x in possible_fields]

    if len(use_fields) == 0:
        return

    # how wide the images list column will be
    filecollen = 20
    # how wide all other output columns will be
    collen = 10
    # floating point printed precision
    floatprec = 4

    # print the column labels if desired
    if print_format:
        headerstrings = {'image': 'IMAGE', 'npix': 'NPIX', 'min': 'MIN',
                         'max': 'MAX', 'mean': 'MEAN', 'midpt': 'MIDPT',
                         'mode': 'MODE', 'stddev': 'STDDEV', 'skew': 'SKEW',
                         'kurtosis': 'KURTOSIS'}
        outstring = '#'
        for ifield in use_fields:
            if ifield == 'image':
                slen = filecollen
            else:
                slen = collen
            outstring += headerstrings[ifield].rjust(slen)
        print(outstring)

    for image in images:
        # open the image
        hdulist = image_open(image)
        if hdulist is None:
            continue

        results = {'npix': 0, 'min': None,
                   'max': None, 'mean': None, 'midpt': None,
                   'mode': None, 'stddev': None, 'skew': None,
                   'kurtosis': None}

        # only use the first (primary) header for the statistics
        data = hdulist[0].data
        npix = 0
        valid = data * 1
        # if there's valid data to look at
        if data is not None:
            data = data.flatten()
            # only grab the points we want to look at
            valid = data[(data >= lower) & (data <= upper)]
            npix = valid.size

            # do the sigma clipping
            for _ in np.arange(nclip):
                if npix > 0:
                    if lsigma > 0.:
                        lowlim = valid.mean() - lsigma * valid.std()
                    else:
                        lowlim = -np.inf
                    if usigma > 0.:
                        upperlim = valid.mean() + usigma * valid.std()
                    else:
                        upperlim = np.inf
                    # don't go below or above previously defined limits
                    lower = max(lower, lowlim)
                    upper = min(upper, upperlim)

                    valid = data[(data >= lower) & (data <= upper)]
                    # no more changes, so stop the loop
                    if valid.size == npix:
                        break
                    npix = valid.size

        # calculate all the statistics for this file
        if npix > 0:
            results['npix'] = npix

            if 'min' in use_fields:
                results['min'] = valid.min()

            if 'max' in use_fields:
                results['max'] = valid.max()

            if 'mean' in use_fields:
                results['mean'] = valid.mean()

            if 'midpt' in use_fields:
                results['midpt'] = np.median(valid)
                # NOTE: the original IRAF definition of median was convoluted
                # and inaccurate. It essentially involved finding the center
                # of a histogram of pixel values. This is way better.

            if 'mode' in use_fields:
                results['mode'] = scipy.stats.mstats.mode(valid)[0][0]
                """
                # the original IRAF calculation of mode. no clue what exactly
                # it is doing.
                # it does seem to have the benefit of returning something near
                # a peak of a distribution in the case of floating point pixels,
                # which are unlikely to have a mode of more than 1.
                # if used, need to include binwidth=0.1 as a keyword argument
                
                mode = None
                hwidth = binwidth * valid.std()
                nbins = (valid.max() - valid.min()) / hwidth + 1
                # does not allow bins < 3 (see ist_ihist)
                if nbins >= 3:
                    # create the histogram based on the binwidth
                    bins = np.arange(valid.min(), valid.max(), hwidth)
                    histo, _ = np.histogram(valid, bins=bins)
                    # Find the bin containing the histogram maximum.
                    bpeak = histo.argmax()
                    # If the maximum is in the first bin return the midpoint
                    # of the bin.
                    if bpeak == 0:
                        mode = valid.min() + 0.5 * hwidth
                    # If the maximum is in the last bin return the midpoint
                    # of the bin.
                    elif bpeak == histo.size - 1:
                        mode = valid.min() + (nbins - 0.5) * hwidth
                    else:
                        # Compute the lower limit of bpeak.
                        bpeak -= 1
                        # Do a parabolic interpolation to find the peak.
                        dh1 = histo[bpeak + 1] - histo[bpeak]
                        dh2 = histo[bpeak + 1] - histo[bpeak + 2]
                        denom = dh1 + dh2
                        if np.isclose(denom, 0.):
                            mode = valid.min() + (bpeak + 0.5) * hwidth
                        else:
                            mode = bpeak + 1 + 0.5 * (dh1 - dh2) / denom
                            mode = valid.min() + (mode - 0.5) * hwidth
                results['mode'] = mode
                """

            if 'stddev' in use_fields:
                results['stddev'] = valid.std()

            if 'skew' in use_fields:
                results['skew'] = scipy.stats.skew(valid)

            if 'kurtosis' in use_fields:
                results['kurtosis'] = scipy.stats.kurtosis(valid)

        # for some reason IRAF does this.
        if print_format:
            outstring = ' '
        else:
            outstring = ''

        # print the statistics requested
        for ii, key in enumerate(use_fields):
            # not enough data to calculate these
            if key != 'image' and results[key] is None:
                tmpstr = 'INDEF'
                if print_format:
                    outstring += tmpstr.rjust(collen)
                else:
                    outstring += tmpstr
                    if ii < len(use_fields) - 1:
                        outstring += '  '
                continue

            if key == 'image':
                if print_format:
                    outstring += image.rjust(filecollen)
                else:
                    outstring += image
            elif key == 'npix':
                if print_format:
                    tmpstr = '{0:{1}d}'.format(results[key], collen)
                    outstring += tmpstr.rjust(collen)
                else:
                    outstring += '{0:d}'.format(results[key])
            else:
                if print_format:
                    tmpstr = '{0:{1}.{2}g}'.format(results[key], collen,
                                                   floatprec)
                    outstring += tmpstr.rjust(collen)
                else:
                    outstring += '{0:g}'.format(results[key])

            if ii < len(use_fields) - 1:
                outstring += '  '
        print(outstring)

        hdulist.close()
    return
