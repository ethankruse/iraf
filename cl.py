import inspect
import os
import csv
from . import uparam_dir
from glob import glob

__all__ = ['Package', 'Parameter', 'loadparams', 'is_iterable', 'file_handler']

# XXX: handle package level CL handling by creating some kind of CL object.
# it'll look like CL (param_list; sub_tasks; parent?). Then loading a package
# will have the same structure with its own sub_tasks. So just work your way
# up the chain for each task. Somehow need to work out how the parent thing
# will work. Just pass in the parent object to init() when creating a new
# object? That should work.
# That could even let me create my own __print__ so it'll print out the CL
# parameters and loaded sub-tasks for a package.
# need to figure out how to be able to index though. e.g.
# iraf.cl.noao.imred.ccdred should print out that CL object, but you need
# to go through the list somehow instead and have tab complete do things right.
# which means adding a __dir__ function? see:
# http://stackoverflow.com/questions/13870241/ipython-tab-completion-for-custom-dict-class
# then I can just wrap everything in my version of the IRAF
# clgstr clgwrd clgeti clgetr clgetb etc

# __init__ makes the original CL Package. Each package needs to have a path
# to it. Then when you import a new package you can figure out where to
# add it along the tree based on the paths.
# also assume that the .par files are in the package's directory.

# need a file with the list of functions that should be contained within a
# package and then use eval() to get them.


def cl():
    return


class Package(dict):
    def __init__(self):
        dict.__init__(self)

    def __dir__(self):
        return self.keys()

    def __getattr__(self, attr):
        return self[attr]


class Parameter(object):
    def __init__(self, value, learn, name):
        self.value = value
        self.learn = learn
        self.name = name


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
    with open(myparamfile, 'r') as ff:
        reader = csv.reader(ff)
        for row in reader:
            # skip over blank lines and comment lines
            if (len(row) == 0 or len(row[0].strip()) == 0 or
                        row[0].strip()[0] == '#'):
                continue
            # make sure we have a complete row
            if len(row) != 7:
                print row
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
                    value = str(value).strip()
                    # check for enumerated list of values
                    if len(dmin):
                        # more than one option available
                        if len(dmin.split('|')) > 1:
                            # make sure it's one of the options
                            enums = [x.strip().upper() for x in dmin.split('|')]
                            if value.upper() not in enums:
                                print 'Input not one of the available options. Try again.'
                                prompt = True
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
                # some other kind of parameter type?
                else:
                    print 'unrecognized parameter type: {0}'.format(dtype)
                    value = None
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

        params[name] = Parameter(value, learn, name)

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
