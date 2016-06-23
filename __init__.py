import os
uparam_dir = os.path.join(os.getcwd(), 'iraf_uparam')
del os

# the subdirectory we put all IRAF packages in
_package_dir = '_packages'

# XXX: somehow give first level packages direct access, e.g. iraf.images, and
# have their package point to the CL one.

# make a 'packages' directory and put everything in there. Then import from
# there using a load_task() function.


def load_package(package):
    """

    Parameters
    ----------
    package

    Returns
    -------

    """
    import re
    # Remove invalid characters
    package = re.sub('[^0-9a-zA-Z_]', '', package)

    # Remove leading characters until we find a letter or underscore
    cleaned = re.sub('^[^a-zA-Z_]+', '', package)
    exec('import {0}.{1}'.format(_package_dir, cleaned), globals(), locals())


def load_task(func, name):
    """
    Load a new function or package, including getting its parameters from
    either the default parameter file or a user file.

    Parameters
    ----------
    func : function
        The function object that will be executed when the Package is called.
    name : str
        What the name of the package or function will be.

    Returns
    -------
    Package
        The Package object for the function or package being loaded.
    """
    import sys
    import inspect
    import csv
    import os
    from _cl import Parameter, Package, convert_value

    # XXX: what happens if this is accidentally called twice for something?

    # path to the function's implementation
    myfile = inspect.getabsfile(func)
    # base of the parameter file name
    parname = '{0}.par'.format(name)
    # full package path, ignoring the leading iraf.
    userparname = func.__module__
    tmp = userparname.split('.')[1:-1]
    if len(tmp) and tmp[0] == _package_dir:
        tmp = tmp[1:]
    userparname = '.'.join(tmp)
    if len(userparname) > 0:
        userparname += '.'
    funcpath = userparname + name
    userparname += name + '.par'
    uparamfile = os.path.join(uparam_dir, userparname)
    # assume for now that the parameter file is in the same directory
    # as the function's source file
    paramfile = os.path.join(os.path.dirname(myfile), parname)
    # if a user parameter file exists, use that instead
    if os.path.exists(uparamfile):
        paramfile = uparamfile

    # read in the default parameters for the function
    defaultparams = []
    if os.path.exists(paramfile):
        with open(paramfile, 'r') as ff:
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
            mode = param[3].lower()
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
    # see if this package has already been loaded
    try:
        curpack = cl
        for ifunc in funcpath.split('.'):
            curpack = curpack[ifunc]
    except KeyError:
        newpack = True
    except NameError:
        newpack = True

    if newpack:
        curpack = Package(func, name)
        # put the new package in the right place in the tree
        try:
            addloc = cl
            for ifunc in funcpath.split('.')[:-1]:
                addloc = addloc[ifunc]
            addloc[name] = curpack
            curpack.__parent__ = addloc
            # if this is a top-level package, put it in the iraf namespace.
            if addloc.__parent__ is None:
                setattr(sys.modules[__name__], name, curpack)
                # globals()[name] = curpack
        except KeyError:
            print('Trying to import function or package out of order.')
            sys.exit(1)
        except NameError:
            # this is the cl package right at the beginning
            setattr(sys.modules[__name__], name, curpack)

    for ii, param in enumerate(defaultparams):
        pname, dtype, mode, default_str, dmin, dmax, prompt_str = param
        mode = mode.lower()
        # get the parameter into the correct data type
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
        # add the parameter to the package
        curpack[pname] = Parameter(default, learn, pname, ii, mode,
                                   default, dmin, dmax, prompt_str, dtype,
                                   default_str, curpack)

        curpack.__nparam__ += 1

    # provide a shortcut access to these parameters from the function itself
    func.__package__ = curpack
    return curpack


# dummy function to start the first import of IRAF
def _packages():
    return


# start the initial IRAF package tree
cl = load_task(_packages, 'cl')

from images import *
from noao import *
from plot import *
# default packages
#load_package('images')
#load_package('noao')
#load_package('plot')
