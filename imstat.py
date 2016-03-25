import inspect
import os
import csv
from . import uparam_dir


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
        if len(row) == 0 or row[0].strip()[0] == '#':
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
                    if len(value) == 0:
                        value = None
                    elif value == 'y' or value == 'yes' or value is True:
                        value = True
                    elif value == 'n' or value == 'no' or value is False:
                        value = False
                    else:
                        print 'Boolean input must be [y/n]. Try again.\n'
                        prompt = True
                # int data types
                elif dtype == 'i':
                    if len(value) == 0 or str(value).upper() == 'INDEF':
                        value = None
                    else:
                        value = int(value)
                    # constrain by min/max values
                    if len(dmax) and len(dmin):
                        if not (int(dmin) <= value <= int(dmax)):
                            print 'Input outside the bounds. Try again.\n'
                            prompt = True
                    # check for enumerated list of values
                    elif len(dmin):
                        # more than one option available
                        if len(dmin.split('|')) > 1:
                            # make sure it's one of the options
                            enums = [int(x) for x in dmin.split('|')]
                            if value not in enums:
                                print 'Input not one of the available options. Try again.\n'
                                prompt = True
                # float data types
                elif dtype == 'r':
                    if len(value) == 0 or str(value).upper() == 'INDEF':
                        value = None
                    else:
                        value = float(value)

                    if len(dmax) and len(dmin):
                        if not (float(dmin) <= value <= float(dmax)):
                            print 'Input outside the bounds. Try again.\n'
                            prompt = True
                    # check for enumerated list of values
                    elif len(dmin):
                        # more than one option available
                        if len(dmin.split('|')) > 1:
                            # make sure it's one of the options
                            enums = [float(x) for x in dmin.split('|')]
                            if value not in enums:
                                print 'Input not one of the available options. Try again.\n'
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
                                print 'Input not one of the available options. Try again.\n'
                                prompt = True
                # filename data types
                elif dtype[0] == 'f':
                    value = str(value)
                    expanded = os.path.expanduser(value)

                    # check for enumerated list of values
                    if len(dmin):
                        # more than one option available
                        if len(dmin.split('|')) > 1:
                            # make sure it's one of the options
                            enums = [os.path.expanduser(x.strip()) for x in dmin.split('|')]
                            if expanded not in enums:
                                print 'Input not one of the available options. Try again.\n'
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
                        print 'Input file must exist. Try again.\n'
                        prompt = True
                    if nonexistent and os.path.exists(expanded):
                        print 'Input file can not already exist. Try again.\n'
                        prompt = True
                    if readacc and not os.access(expanded, os.R_OK):
                        print 'Input file must have read access. Try again.\n'
                        prompt = True
                    if writeacc and not os.access(expanded, os.W_OK):
                        print 'Input file must have write access. Try again.\n'
                        prompt = True

            except ValueError:
                print 'Could not interpret input. Try again.\n'
                prompt = True

            if not prompt:
                break

            value = raw_input('{0} {1}{2}: '.format(prompt_str, allowrange, default_str))

            if len(value.strip()) == 0:
                value = default
            prompt = False

        params[name] = Parameter(value, learn)

    return params


# XXX: make a relatively simple file interpretation system.
# so it can handle iraf things like @lists, wildcards, etc.
# probably use it in loadparams() for the filetype thing (see imslice.par)
# also allow for CSV separated lists of images? Also iterable Pythonic lists
# as inputs

def imstatistics(*args, **kwargs):
    params = loadparams(*args, **kwargs)
    images = params['images'].value
    fields = params['fields'].value
    lower = params['lower'].value
    upper = params['upper'].value
    nclip = params['nclip'].value
    lsigma = params['lsigma'].value
    usigma = params['usigma'].value
    binwidth = params['binwidth'].value
    print_format = params['format'].value
    cache = params['cache'].value

    possible_fields = "image|npix|min|max|mean|midpt|mode|stddev|skew|kurtosis".split('|')

    fcol = '%10d'
    fint = '%10d'
    fflt = '%10.4f'
    fstr = '%20s'

    fields = [x.strip().lower() for x in fields.split(',')]
    # retain the same order as possible_fields, but only the ones requested
    combined = [x for x in possible_fields if x in fields]

    if len(combined) == 0:
        return

    if print_format:
        headerstrings = {'image': 'IMAGE', 'npix': 'NPIX', 'min': 'MIN',
                         'max': 'MAX', 'mean': 'MEAN', 'midpt': 'MIDPT',
                         'mode': 'MODE', 'stddev': 'STDDEV', 'skew': 'SKEW',
                         'kurtosis': 'KURTOSIS'}
        outstring = '#'
        for ifield in combined:
            if ifield == 'image':
                slen = 20
            else:
                slen = 10
            outstring += headerstrings[ifield].ljust(slen)
        outstring += '\n'
        print outstring





    print combined



    return params


