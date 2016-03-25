import inspect
import os
import csv
from . import uparam_dir


class Parameter(object):
    def __init__(self, name, value, learn):
        self.name = name
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

        while True:

            print value
            # move on if this is supposed to be undefined
            if str(value).upper() == 'INDEF':
                value = None

            # XXX: what if there is no input so default == ''

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
                    else:
                        print 'Boolean input must be [y/n]. Try again.\n'
                        prompt = True
                # int data types
                elif dtype == 'i':
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

                    # XXX: finish this

                    # XXX: documentation says that min/max field is valid for files?

            except ValueError:
                print 'Could not interpret input. Try again.\n'
                prompt = True

            if not prompt:
                break

            value = raw_input('{0} {1}{2}: '.format(prompt_str, allowrange, default_str))

            if len(value.strip()) == 0:
                value = default
            prompt = False

        print param

    return params


def imstatistics(*args, **kwargs):
    return loadparams(*args, **kwargs)


