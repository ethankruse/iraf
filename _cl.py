from __future__ import print_function
import inspect
import os
import csv
from . import uparam_dir
from glob import glob
import sys
import shutil

__all__ = ['Package', 'Parameter', 'is_iterable', 'file_handler',
           'clget', 'convert_value']


class Package(dict):
    def __init__(self, func, name):
        # necessary to allow packages/functions to have arbitrary parameters
        # and be accessed via dot autocomplete
        # See http://stackoverflow.com/questions/4984647/
        # accessing-dict-keys-like-an-attribute-in-python
        super(Package, self).__init__()
        self.__dict__ = self
        self.__name__ = name
        self.__parent__ = None
        self.__nparam__ = 0
        self.__function__ = func
        self.__doc__ = func.__doc__

    # again to allow parameters to be accessed via dot autocomplete
    def __dir__(self):
        return self.keys()

    def __setattr__(self, key, value):
        # allow for self.param = 5 to set the parameter's value and not
        # erase the parameter object with an int
        if hasattr(self, key) and isinstance(getattr(self, key), Parameter):
            getattr(self, key).value = value
        else:
            # default action if we're not changing a parameter's value
            super(Package, self).__setattr__(key, value)

    def __repr__(self, level=0):
        # print a tree of all packages below this one showing the structure
        # of everything that is loaded
        indent = '  ' * level
        ret = indent + self.__name__ + '\n'
        attribs = dir(self)
        for attrib in attribs:
            if attrib != '__parent__' and isinstance(getattr(self, attrib),
                                                     Package):
                ret += getattr(self, attrib).__repr__(level=level + 1)
        return ret

    def __call__(self, *args, **kwargs):
        startfunc(self.__function__, *args, **kwargs)
        retval = self.__function__()
        endfunc(self.__function__, self.__name__)
        return retval


class Parameter(object):
    def __init__(self, value, learn, name, order, mode, default, dmin,
                 dmax, prompt_str, dtype, default_str, parent):
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
        self.__package__ = parent

    def __setattr__(self, key, value):
        valid = True
        if key == 'value':
            value, valid = convert_value(value, self.dtype, self.min, self.max)
            if self.default != value and valid:
                self.changed = True
            # if we're changing it back to the default, no need to rewrite
            else:
                self.changed = False
        if valid:
            # actually set the attribute as normal
            super(Parameter, self).__setattr__(key, value)

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

    def update(self, value=None):
        """
        Update this parameter's value, prompting the user for input
        if necessary.

        Parameters
        ----------
        value : any type, optional
            What to update the parameter to. If not supplied, the parameter
            prompt will be displayed asking for user input.

        Returns
        -------

        """
        update_parameter(self, automode='q', value=value)


def convert_value(value, dtype, dmin, dmax):
    """
    Convert an input value (typically a string) to the appropriate type for
    the parameter.

    Also check to see if the input is valid (e.g. falls within min/max values).

    Parameters
    ----------
    value : any type
        Input value to set. Will be converted to the appropriate type.
    dtype : {'b', 'i', 'r', 's', 'f'}
        What type of parameter this is. Boolean, integer, real (float), string,
        or file name (also a string, but expected to be a file path).
    dmin : str
        Can either be an enumerated list of legal values (in which case dmax
        must be None) or the lower bounds for int/real data types. Can also be
        None.
    dmax : str
        Upper bound for int/real data types.

    Returns
    -------
    value : data type
        The input value converted to the appropriate type.
    valid : bool
        Whether or not the input could be interpreted correctly and falls
        within all bounds.

    """
    valid = True
    # try to convert to the appropriate type
    try:
        # None breaks things, so pass it on through
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
                    print(
                        'Input outside the bounds ({0}-{1}). Try again.'.format(
                            int(dmin), int(dmax)))
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
            # constrain by min/max values
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

            # check for these flags
            readacc = False
            writeacc = False
            nonexistent = False
            exists = False
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


def update_parameter(param, automode=None, value=None):
    """
    Update the value of a parameter, asking the user for input if necessary.

    Parameters
    ----------
    param : Parameter
        Which parameter to update.
    automode : string, optional
        What to replace mode 'a' with. If None, search to find the appropriate
        mode.
    value : any type, optional
        What the assumed user input is by default.

    Returns
    -------

    """
    if automode is None:
        # get the appropriate automode parameter for this function
        automode = clget(param.__package__.__function__, 'mode').value

    # fill in automode if necessary
    mode = ''
    for char in param.mode:
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

    # if we're given a starting user input value, turn off the prompting
    if value is not None:
        prompt = False
    else:
        value = param.value

    # prompt the user until valid input is given
    while True:
        value, valid = convert_value(value, param.dtype, param.min, param.max)

        if not valid:
            prompt = True

        if not prompt:
            break

        value = raw_input('{0}: '.format(param))

        # if the user just hit enter, accept the default value.
        if len(value.strip()) == 0:
            value = param.value
        prompt = False

    param.value = value
    param.learn = learn


def startfunc(func, *args, **kwargs):
    """
    Initialize an IRAF task. Must be called at the beginning of any
    IRAF package task's implementation.

    Run the checks to make sure all parameters have been properly initalized,
    asking the user for input for any missing parameters.

    Parameters
    ----------
    func : function
        Function object of the task we're initalizing.
    args : list
        Arguments given to the function to be interpreted as parameter values.
    kwargs : dict
        Keyword arguments given to the function to be interpreted as
         parameter values.

    Returns
    -------

    """
    curpack = func.__package__

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
        value = None
        # kwargs take precedent over positional args
        if iparam.name in kwargs:
            value = kwargs[iparam.name]
        elif ii < len(args):
            value = args[ii]

        update_parameter(iparam, automode=automode, value=value)

    return


def endfunc(func, name):
    """
    Clean up parameter values after an IRAF task and learn any new parameter
    values.

    This must be the last function called in an IRAF task. Learn all changed
    parameters with the learn flag set, and reset all those without the learn
    flag. Save the parameter file if necessary.

    Parameters
    ----------
    func : function
        Function object of the task we're ending.
    name : str
        Name of the function we're ending.
    Returns
    -------

    """
    # base of the parameter file name
    parname = '{0}.par'.format(name)
    # where the function is implemented
    userparname = func.__module__
    # ignore the iraf. and name of the file the function is in
    userparname = '.'.join(userparname.split('.')[1:-1])
    if len(userparname) > 0:
        userparname += '.'
    # now have the full parameter file name, need to add the directories next
    userparname += parname
    # assume for now that the parameter file is in the same directory
    # as the function's source file
    fileloc = inspect.getabsfile(func)
    # default parameter file location
    paramfile = os.path.join(os.path.dirname(fileloc), parname)
    # user parameter file location
    uparamfile = os.path.join(uparam_dir, userparname)

    curpack = func.__package__

    possible = dir(curpack)
    # figure out which parameters have changed their values
    changed = []
    for key in possible:
        if isinstance(curpack[key], Parameter) and curpack[key].changed:
            changed.append(key)

    # which changed parameters do we need to learn?
    learned = []
    for key in changed:
        if curpack[key].learn:
            learned.append(key)
            # change the default value in memory
            curpack[key].default = curpack[key].value
        else:
            # reset values we don't want to learn
            curpack[key].value = curpack[key].default

    # change the default value on disk
    if len(learned):
        # copy over the parameter file to the user directory if necessary
        if not os.path.exists(uparamfile):
            if not os.path.exists(uparam_dir):
                os.mkdir(uparam_dir)
            shutil.copy2(paramfile, uparamfile)

        tmp = uparamfile + 'tmp'
        shutil.copy2(uparamfile, tmp)

        wf = open(uparamfile, 'w')
        writer = csv.writer(wf)
        # read in the old version of the file
        with open(tmp, 'r') as ff:
            reader = csv.reader(ff)
            for row in reader:
                # skip comment or blank rows
                if (len(row) == 0 or len(row[0].strip()) == 0 or
                            row[0].strip()[0] == '#'):
                    pass
                # if we are learning this parameter, replace the default value.
                elif row[0].strip() in learned:
                    row[3] = curpack[row[0].strip()].value
                writer.writerow(row)
        wf.close()
        os.remove(tmp)
    return


def clget(func, param):
    """
    Retrieve a parameter from within the scope of the calling function.

    Search up the IRAF package tree until a package with the
    parameter in question is found.

    Parameters
    ----------
    func : function
        The function object asking for a parameter.
    param : str
        Name of the parameter.

    Returns
    -------
    Parameter
        Parameter object with the name in question.

    """
    pkg = func.__package__

    # search up the tree for a parameter with the right name
    obj = None
    while pkg is not None:
        if param in pkg.keys() and isinstance(pkg[param], Parameter):
            obj = pkg[param]
            break
        pkg = pkg.__parent__

    # XXX: replace all sys.exit() calls with exceptions
    if obj is None:
        print('Could not find parameter {0} in {1}:{2}'.format(param,
                                                               func.__module__,
                                                               pkg.__name__))
        sys.exit(1)

    return obj


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
        # XXX: deal with wildcards here too?
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
            # XXX: deal with lists of '@' files?
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
