from __future__ import print_function
import inspect
import os
import csv
from . import uparam_dir
from glob import glob
import sys
import shutil

__all__ = ['Package', 'Parameter', 'loadparams', 'is_iterable', 'file_handler',
           'cl', 'clget', 'startfunc', 'endfunc']


class Package(dict):
    def __init__(self, name):
        dict.__init__(self)
        self.__name__ = name
        self.__parent__ = None
        self.__nparam__ = 0

    def __dir__(self):
        return self.keys()

    def __getattr__(self, attr):
        return self[attr]

        # XXX: implement a print tree type function


# XXX: set it up so that setting iraf.cl.parameter = 5 will set the value to 5?
# this could be done in the Package __setattr__ function!
class Parameter(object):
    def __init__(self, value, learn, name, order, mode, default, dmin,
                 dmax, prompt_str, dtype, default_str):
        self.default = default
        self.dtype = dtype
        self.min = dmin
        self.max = dmax
        self.value = value
        self.learn = learn
        self.name = name
        self.order = order
        self.mode = mode
        self.default_str = default_str
        self.prompt = prompt_str
        self.changed = False

    def __setattr__(self, key, value):
        if key == 'value':
            value, valid = convert_value(value, self.dtype, self.min, self.max)
            if self.default != value:
                self.changed = True
            # if we're changing it back to the default, no need to rewrite
            else:
                self.changed = False
        object.__setattr__(self, key, value)

    def __repr__(self):
        return 'Parameter(name={0!r}, value={1!r}, default={2!r})'.format(
            self.name, self.value, self.default)

    def __str__(self):
        allowrange = ''
        if len(self.max) and len(self.min):
            allowrange = '<{0} to {1}>'.format(self.min, self.max)
        elif len(self.min) and len(self.min.split('|')) > 1:
            allowrange = '<{0}>'.format(self.min)
        if self.changed:
            default_str = '[{0}]'.format(self.value)
        elif len(self.default_str) > 0:
            default_str = '[{0}]'.format(self.default_str)
        else:
            default_str = ''
        return '{0}|{1} {2}{3}'.format(self.name, self.prompt, allowrange,
                                       default_str)


def convert_value(value, dtype, dmin, dmax):
    valid = True
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
                print('Boolean input must be [y/n]. Try again.')
                valid = False
        # int data types
        elif dtype == 'i':
            if len(str(value)) == 0 or str(value).upper() == 'INDEF':
                value = None
            else:
                value = int(value)
            # constrain by min/max values
            if len(dmax) and len(dmin):
                if not (int(dmin) <= value <= int(dmax)):
                    print('Input outside the bounds. Try again.')
                    valid = False
            # check for enumerated list of values
            elif len(dmin):
                # more than one option available
                if len(dmin.split('|')) > 1:
                    # make sure it's one of the options
                    enums = [int(x) for x in dmin.split('|')]
                    if value not in enums:
                        print(
                            'Input not one of the available options. Try again.')
                        valid = False
        # float data types
        elif dtype == 'r':
            if len(str(value)) == 0 or str(value).upper() == 'INDEF':
                value = None
            else:
                value = float(value)

            if len(dmax) and len(dmin):
                if not (float(dmin) <= value <= float(dmax)):
                    print('Input outside the bounds. Try again.')
                    valid = False
            # check for enumerated list of values
            elif len(dmin):
                # more than one option available
                if len(dmin.split('|')) > 1:
                    # make sure it's one of the options
                    enums = [float(x) for x in dmin.split('|')]
                    if value not in enums:
                        print(
                            'Input not one of the available options. Try again.')
                        valid = False
        # string data types
        elif dtype == 's':
            value = str(value).strip()
            # check for enumerated list of values
            if len(dmin):
                # more than one option available
                if len(dmin.split('|')) > 1:
                    # make sure it's one of the options
                    enums = [x.strip().upper() for x in dmin.split('|')]
                    if value.upper() not in enums:
                        print(
                            'Input not one of the available options. Try again.')
                        valid = False
            # XXX: should we just return an empty string?
            if len(value) == 0:
                value = None

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
                        print(
                            'Input not one of the available options. Try again.')
                        valid = False

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
                print('Input file must exist. Try again.')
                valid = False
            if nonexistent and os.path.exists(expanded):
                print('Input file can not already exist. Try again.')
                valid = False
            if readacc and not os.access(expanded, os.R_OK):
                print('Input file must have read access. Try again.')
                valid = False
            if writeacc and not os.access(expanded, os.W_OK):
                print('Input file must have write access. Try again.')
                valid = False
        # some other kind of parameter type?
        else:
            print('unrecognized parameter type: {0}'.format(dtype))
            value = None
    except ValueError:
        print('Could not interpret input. Try again.')
        valid = False

    return value, valid


def startfunc(func, *args, **kwargs):
    modname = func.__module__
    # remove the iraf. and the name of the .py file the code is in
    modname = '.'.join(modname.split('.')[1:-1])
    if len(modname) > 0:
        modname += '.'
    modname += func.__name__

    curpack = cl
    tree = modname.split('.')
    try:
        for branch in tree:
            curpack = curpack[branch]
    except KeyError:
        print('Function/package {0} not found.'.format(modname))
        sys.exit(1)

    # get the appropriate automode parameter for this function
    automode = clget(func, 'mode').value
    # XXX: need to deal with 'menu' mode somehow.

    # set up the parameter list
    plist = [None for _ in range(curpack.__nparam__)]

    # fill in the empty list with the parameters in order
    keylist = curpack.keys()
    for key in keylist:
        if isinstance(curpack[key], Parameter):
            plist[curpack[key].order] = curpack[key]

    for ii, iparam in enumerate(plist):
        name = iparam.name
        dmin = iparam.min
        dmax = iparam.max
        dtype = iparam.dtype
        value = iparam.value
        default = iparam.value
        mode = ''

        for char in iparam.mode:
            if char == 'a':
                mode += automode
            else:
                mode += char

        prompt = False
        if 'q' in mode:
            prompt = True

        learn = False
        if 'l' in mode:
            learn = True

        # kwargs take precedent over positional args
        if name in kwargs:
            value = kwargs[name]
            prompt = False
        elif ii < len(args):
            value = args[ii]
            prompt = False

        while True:
            value, valid = convert_value(value, dtype, dmin, dmax)

            if not valid:
                prompt = True

            if not prompt:
                break

            value = raw_input('{0}: '.format(iparam))

            if len(value.strip()) == 0:
                value = default
            prompt = False

        iparam.value = value
        iparam.learn = learn

    return


def endfunc(func):
    # path to the function's implementation
    myfile = inspect.getabsfile(func)
    # base of the parameter file name
    parname = '{0}.par'.format(func.__name__)
    userparname = func.__module__
    userparname = '.'.join(userparname.split('.')[1:-1])
    if len(userparname) > 0:
        userparname += '.'
    funcpath = userparname + func.__name__
    userparname += func.__name__ + '.par'
    # assume for now that the parameter file is in the same directory
    # as the function's source file
    myparamfile = os.path.join(os.path.dirname(myfile), parname)
    uparamfile = os.path.join(uparam_dir, userparname)

    try:
        curpack = cl
        for ifunc in funcpath.split('.'):
            curpack = curpack[ifunc]
    except KeyError:
        print('Could not find function {0} in endfunc'.format(funcpath))
        sys.exit(1)

    possible = dir(curpack)
    changed = []
    for key in possible:
        if isinstance(curpack[key], Parameter) and curpack[key].changed:
            changed.append(key)

    # need to check for the 'learn' parameter
    learned = []
    for key in changed:
        if curpack[key].learn:
            learned.append(key)
        else:
            curpack[key].value = curpack[key].default

    if len(learned):
        if not os.path.exists(uparamfile):
            if not os.path.exists(uparam_dir):
                os.mkdir(uparam_dir)
            shutil.copy2(myparamfile, uparamfile)

        tmp = uparamfile + 'tmp'
        shutil.copy2(uparamfile, tmp)

        wf = open(uparamfile, 'w')
        writer = csv.writer(wf)
        with open(tmp, 'r') as ff:
            reader = csv.reader(ff)
            for row in reader:
                if (len(row) == 0 or len(row[0].strip()) == 0 or
                        row[0].strip()[0] == '#'):
                    pass
                elif row[0].strip() in learned:
                    row[3] = curpack[row[0].strip()].value
                writer.writerow(row)
        wf.close()
        os.remove(tmp)
    return


def clget(func, param):
    modname = func.__module__
    # remove the iraf. and the name of the .py file the code is in
    modname = '.'.join(modname.split('.')[1:-1])
    if len(modname) > 0:
        modname += '.'
    modname += func.__name__

    tree = modname.split('.')
    pkg = cl
    try:
        for branch in tree:
            pkg = pkg[branch]
    except KeyError:
        print('Function/package {0} not found.'.format(modname))
        sys.exit(1)

    obj = None
    while pkg is not None:
        if param in pkg.keys() and isinstance(pkg[param], Parameter):
            obj = pkg[param]
            break
        pkg = pkg.__parent__

    if obj is None:
        print('Could not find parameter {0} in {1}'.format(param, modname))
        sys.exit(1)

    return obj


def loadparams(func):
    # path to the function's implementation
    myfile = inspect.getabsfile(func)
    # base of the parameter file name
    parname = '{0}.par'.format(func.__name__)
    userparname = func.__module__
    userparname = '.'.join(userparname.split('.')[1:-1])
    if len(userparname) > 0:
        userparname += '.'
    funcpath = userparname + func.__name__
    userparname += func.__name__ + '.par'
    # assume for now that the parameter file is in the same directory
    # as the function's source file
    uparamfile = os.path.join(uparam_dir, userparname)
    myparamfile = os.path.join(os.path.dirname(myfile), parname)
    # if a user parameter file exists, use that instead
    if os.path.exists(uparamfile):
        myparamfile = uparamfile

    # read in the default parameters for the function
    defaultparams = []
    if os.path.exists(myparamfile):
        with open(myparamfile, 'r') as ff:
            reader = csv.reader(ff)
            for row in reader:
                # skip over blank lines and comment lines
                if (len(row) == 0 or len(row[0].strip()) == 0 or
                            row[0].strip()[0] == '#'):
                    continue
                # make sure we have a complete row
                assert len(row) == 7
                defaultparams.append([x.strip() for x in row])

    # if we're setting the mode, get it in the right format
    for param in defaultparams:
        if param[0] == 'mode':
            mode = param[3]
            # get modes into the single letter categories
            if mode == 'auto':
                mode = 'a'
            elif mode == 'hidden':
                mode = 'h'
            elif mode == 'learn':
                mode = 'l'
            elif mode == 'query':
                mode = 'q'
            elif mode == 'menu':
                mode = 'm'
            param[3] = mode

    newpack = False
    # if cl has already been initialized
    if isinstance(cl, Package):
        try:
            curpack = cl
            for ifunc in funcpath.split('.'):
                curpack = curpack[ifunc]
        except KeyError:
            newpack = True
    else:
        newpack = True

    if newpack:
        curpack = Package(func.__name__)
        addloc = cl
        try:
            for ifunc in funcpath.split('.')[:-1]:
                addloc = addloc[ifunc]
            addloc[func.__name__] = curpack
            curpack.__parent__ = addloc
        except KeyError:
            print('Trying to import function or package out of order.')
            sys.exit(1)
        except TypeError:
            pass

    for ii, param in enumerate(defaultparams):
        name, dtype, mode, default_str, dmin, dmax, prompt_str = param

        default, valid = convert_value(default_str, dtype, dmin, dmax)
        # get modes into the single letter categories
        if mode == 'auto':
            mode = 'a'
        elif mode == 'hidden':
            mode = 'h'
        elif mode == 'learn':
            mode = 'l'
        elif mode == 'query':
            mode = 'q'
        elif mode == 'menu':
            mode = 'm'

        learn = False
        if 'l' in mode:
            learn = True

        curpack[name] = Parameter(default, learn, name, ii, mode,
                                  default, dmin, dmax, prompt_str, dtype,
                                  default_str)

        curpack.__nparam__ += 1

    return curpack


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

    # XXX: should we return None or an empty list?
    if filelist is None:
        return []

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

            with open(fname, 'r') as ff:
                reader = csv.reader(ff)
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


def cl():
    return


cl = loadparams(cl)
