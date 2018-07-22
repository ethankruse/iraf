import os
import csv
from glob import glob

__all__ = ['is_iterable', 'file_handler']


def is_iterable(obj):
    """
    Returns True if the input is iterable, but not a string. Returns
    False otherwise.
    """
    from collections import Iterable
    return not isinstance(obj, str) and isinstance(obj, Iterable)


def file_handler(filelist):
    """
    A limited version of IRAF's imtopen.

    General purpose file type interpreter. Takes IRAF file inputs and returns
    list with actual paths for all matching files.

    Allows for a file with a list of file patterns by starting with '@', e.g.
    '@filelist.txt'. Will recursively search each line in an @list so e.g.
    it can handle @lists of @lists.

    All file names can also accept ~ and wildcard (* ?) expansions as well.

    Currently does not support IRAF's image section format.
    e.g. image.fits[1:50, 1:50] will fail rather than allowing only a portion
    of an image to be opened.

    Parameters
    ----------
    filelist : string or iterable
        String or list of input file strings to expand and create the final
        list of files.

    Returns
    -------
    list[str]
        List of file names matching all patterns.
    """

    if filelist is None:
        return []
    # make the input string a list
    if not is_iterable(filelist):
        filelist = [filelist]

    outlist = []
    # go through every input and add to the output list of files
    for istr in filelist:
        istr = istr.strip()
        # we have an input file with a list of files to use.
        if istr[0] == '@':
            # remove the @ and deal with home directories
            fnames = os.path.expanduser(istr[1:])
            # deal with wildcards
            if '*' in fnames or '?' in fnames:
                fnames = glob(fnames)
            else:
                fnames = [fnames]
            # load each requested list
            for fname in fnames:
                with open(fname, 'r') as ff:
                    reader = csv.reader(ff)
                    for row in reader:
                        # skip over blank lines and comment lines
                        if (len(row) == 0 or len(row[0].strip()) == 0 or
                                row[0].strip()[0] == '#'):
                            continue
                        # recursively handle these files
                        newfiles = file_handler(row)
                        for newfile in newfiles:
                            outlist.append(newfile)
        # just a single file name, no list
        else:
            fnames = os.path.expanduser(istr)
            # deal with wildcards
            if '*' in fnames or '?' in fnames:
                fnames = glob(fnames)
            else:
                fnames = [fnames]
            # add the results to the output list
            for fname in fnames:
                outlist.append(fname)

    return outlist
