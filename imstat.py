import inspect
import os
import csv
from . import uparam_dir
from glob import glob
from astropy.io import fits
import numpy as np
import scipy.stats


class Parameter(object):
    def __init__(self, value, learn):
        self.value = value
        self.learn = learn


def loadparams(*args, **kwargs):
    # name of the function that is loading parameters
    myname = inspect.stack()[1][3]
    # and its source file
    myfile = inspect.stack()[1][1]
    # base of the parameter file name
    parname = '{0}.par'.format(myname)
    # assume for now that the parameter file is in the same directory
    # as the function's source file
    myparamfile = os.path.join(os.path.dirname(myfile), parname)
    # if a user parameter file exists, use that instead
    if os.path.exists(os.path.join(uparam_dir, parname)):
        myparamfile = os.path.join(uparam_dir, parname)

    # read in the default parameters for the function
    defaultparams = []
    reader = csv.reader(open(myparamfile, 'r'))
    for row in reader:
        # skip over blank lines and comment lines
        if (len(row) == 0 or len(row[0].strip()) == 0 or
                    row[0].strip()[0] == '#'):
            continue
        # make sure we have a complete row
        assert len(row) == 7
        defaultparams.append([x.strip() for x in row])

    automode = 'h'
    for param in defaultparams:
        if param[0] == 'mode':
            mode = param[3]
            # get modes into the single letter categories
            if mode == 'hidden':
                mode = 'h'
            elif mode == 'learn':
                mode = 'l'
            elif mode == 'query':
                mode = 'q'
            if mode in ['h', 'l', 'q', 'ql', 'hl', 'lh', 'lq']:
                automode = mode
            else:
                automode = 'h'

    params = {}
    for ii, param in enumerate(defaultparams):
        name, dtype, mode, default, dmin, dmax, prompt_str = param

        # first switch auto into the appropriate bin
        if mode == 'auto' or mode == 'a':
            mode = automode
        # get modes into the single letter categories
        if mode == 'hidden':
            mode = 'h'
        elif mode == 'learn':
            mode = 'l'
        elif mode == 'query':
            mode = 'q'

        learn = False
        if mode in ['l', 'ql', 'hl', 'lh', 'lq']:
            learn = True

        prompt = False
        if mode in ['q', 'ql', 'lq']:
            prompt = True

        # kwargs take precedent over positional args
        if name in kwargs:
            value = kwargs[name]
            prompt = False
        elif ii < len(args):
            value = args[ii]
            prompt = False
        else:
            value = default

        if len(prompt_str) == 0:
            prompt_str = name

        allowrange = ''
        if len(dmax) and len(dmin):
            allowrange = '<{0} to {1}>'.format(dmin, dmax)
        elif len(dmin) and len(dmin.split('|')) > 1:
            allowrange = '<{0}>'.format(dmin)

        if len(default) > 0:
            default_str = '[{0}]'.format(default)
        else:
            default_str = ''

        while True:
            # try to convert to the appropriate type
            try:
                if value is None:
                    pass
                # boolean data types
                elif dtype == 'b':
                    if value == 'y' or value == 'yes' or value is True:
                        value = True
                    elif value == 'n' or value == 'no' or value is False:
                        value = False
                    elif len(str(value)) == 0:
                        value = None
                    else:
                        print 'Boolean input must be [y/n]. Try again.'
                        prompt = True
                # int data types
                elif dtype == 'i':
                    if len(str(value)) == 0 or str(value).upper() == 'INDEF':
                        value = None
                    else:
                        value = int(value)
                    # constrain by min/max values
                    if len(dmax) and len(dmin):
                        if not (int(dmin) <= value <= int(dmax)):
                            print 'Input outside the bounds. Try again.'
                            prompt = True
                    # check for enumerated list of values
                    elif len(dmin):
                        # more than one option available
                        if len(dmin.split('|')) > 1:
                            # make sure it's one of the options
                            enums = [int(x) for x in dmin.split('|')]
                            if value not in enums:
                                print 'Input not one of the available options. Try again.'
                                prompt = True
                # float data types
                elif dtype == 'r':
                    if len(str(value)) == 0 or str(value).upper() == 'INDEF':
                        value = None
                    else:
                        value = float(value)

                    if len(dmax) and len(dmin):
                        if not (float(dmin) <= value <= float(dmax)):
                            print 'Input outside the bounds. Try again.'
                            prompt = True
                    # check for enumerated list of values
                    elif len(dmin):
                        # more than one option available
                        if len(dmin.split('|')) > 1:
                            # make sure it's one of the options
                            enums = [float(x) for x in dmin.split('|')]
                            if value not in enums:
                                print 'Input not one of the available options. Try again.'
                                prompt = True
                # string data types
                elif dtype == 's':
                    value = str(value)
                    # check for enumerated list of values
                    if len(dmin):
                        # more than one option available
                        if len(dmin.split('|')) > 1:
                            # make sure it's one of the options
                            enums = [x.strip().upper() for x in dmin.split('|')]
                            if value.upper() not in enums:
                                print 'Input not one of the available options. Try again.'
                                prompt = True
                # filename data types
                elif dtype[0] == 'f':
                    # XXX: do we want to call file_handler on this one?
                    value = str(value)
                    expanded = os.path.expanduser(value)

                    # check for enumerated list of values
                    if len(dmin):
                        # more than one option available
                        if len(dmin.split('|')) > 1:
                            # make sure it's one of the options
                            enums = [os.path.expanduser(x.strip()) for x in
                                     dmin.split('|')]
                            if expanded not in enums:
                                print 'Input not one of the available options. Try again.'
                                prompt = True

                    # XXX: documentation says that min/max field is valid
                    #  for files?

                    readacc = False
                    writeacc = False
                    nonexistent = False
                    exists = False
                    # should we check for these categories
                    if len(dtype) > 1:
                        if 'r' in dtype[1:]:
                            readacc = True
                        if 'w' in dtype[1:]:
                            writeacc = True
                        if 'n' in dtype[1:]:
                            nonexistent = True
                        if 'e' in dtype[1:]:
                            exists = True

                    if exists and not os.path.exists(expanded):
                        print 'Input file must exist. Try again.'
                        prompt = True
                    if nonexistent and os.path.exists(expanded):
                        print 'Input file can not already exist. Try again.'
                        prompt = True
                    if readacc and not os.access(expanded, os.R_OK):
                        print 'Input file must have read access. Try again.'
                        prompt = True
                    if writeacc and not os.access(expanded, os.W_OK):
                        print 'Input file must have write access. Try again.'
                        prompt = True

            except ValueError:
                print 'Could not interpret input. Try again.'
                prompt = True

            if not prompt:
                break

            value = raw_input('{0} {1}{2}: '.format(prompt_str, allowrange,
                                                    default_str))

            if len(value.strip()) == 0:
                value = default
            prompt = False

        params[name] = Parameter(value, learn)

    return params


def is_iterable(obj):
    """
    Returns True if the input is iterable, but not a string. Returns
    False otherwise.
    """
    from collections import Iterable
    return not isinstance(obj, str) and isinstance(obj, Iterable)


def file_handler(filelist):
    """
    General purpose file type interpreter. Takes IRAF file inputs and returns
    list with actual paths for all matching files.

    Allows for a file with a list of file patterns by starting with '@', e.g.
    '@filelist.txt'. Input string can also be comma separated list of filenames.
    All file names can also accept ~ and wildcard (* ?) expansions as well.

    Parameters
    ----------
    filelist : string or iterable
        String or list of input file strings to expand and create the final
        list of files.

    Returns
    -------
    outlist : list
        List of file names matching all patterns.
    """

    # XXX: this does not allow for image subsections,
    #  e.g. imagename[x1:x2,y1:y2]
    # How should that be handled?

    # see if the input is already a list
    is_list = is_iterable(filelist)
    # make the input string a list
    if not is_list:
        filelist = [filelist]

    outlist = []
    # go through every input and create the output list of files
    for istr in filelist:
        istr = istr.strip()
        # every file here to look for
        files = []

        # we have an input file with a list of files to use.
        if istr[0] == '@':
            fname = os.path.expanduser(istr[1:])

            reader = csv.reader(open(fname, 'r'))
            for row in reader:
                # skip over blank lines and comment lines
                if (len(row) == 0 or len(row[0].strip()) == 0 or
                            row[0].strip()[0] == '#'):
                    continue

                files.append(os.path.expanduser(row[0].strip()))
        else:
            # this is either a single file or a CSV list of files
            tmp = istr.split(',')
            for ifile in tmp:
                files.append(os.path.expanduser(ifile.strip()))

        # now that we have separated this up into individual path names or
        # loaded them from the input list, deal with wildcards
        for ifile in files:
            # search for the matching files
            results = glob(ifile)
            # add any matches to the output
            for jfile in results:
                outlist.append(jfile)

    return outlist


def imstatistics(*args, **kwargs):
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
    binwidth : float
        Bin width of histogram (in sigma). Currently obsolete; was only used
        in the original IRAF implementation to calculate the mode.
    format : boolean
        Whether or not to format the output and print column labels
    cache : boolean
        Obsolete. Was previously used to ask about caching the images
        in memory.
    """
    params = loadparams(*args, **kwargs)
    # list of images to run on
    images = params['images'].value
    images = file_handler(images)
    # which statistics we want to calculate
    fields = params['fields'].value
    # lower/upper limits for pixel values
    lower = params['lower'].value
    if lower is None:
        lower = -np.inf
    upper = params['upper'].value
    if upper is None:
        upper = np.inf
    # number of clipping iterations
    nclip = params['nclip'].value
    # lower and upper sigma clipping level
    lsigma = params['lsigma'].value
    usigma = params['usigma'].value
    # bin width of histogram in sigma
    # this is only used in the original IRAF method of defining mode.
    # could now be eliminated.
    binwidth = params['binwidth'].value
    # should the output be pretty formatted
    print_format = params['format'].value
    # cache image in memory (obsolete)
    cache = params['cache'].value

    possible_fields = "image|npix|min|max|mean|midpt|mode|stddev|skew|kurtosis"
    possible_fields = possible_fields.split('|')

    in_fields = [x.strip().lower() for x in fields.split(',')]
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
        print outstring

    for image in images:
        # open the image
        try:
            hdulist = fits.open(image)
        except IOError:
            print "Error reading image {0} ...".format(image)
            continue

        results = {'npix': 0, 'min': None,
                   'max': None, 'mean': None, 'midpt': None,
                   'mode': None, 'stddev': None, 'skew': None,
                   'kurtosis': None}

        # only use the first (primary) header for the statistics
        data = hdulist[0].data
        npix = 0
        # if there's valid data to look at
        if data is not None:
            data = data.flatten()
            # only grab the points we want to look at
            valid = data[(data >= lower) & (data <= upper)]
            npix = valid.size

            # do the sigma clipping
            for ii in np.arange(nclip):
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
        print outstring

        hdulist.close()

    # XXX: need to deal with 'learning' some parameters
    
    return
