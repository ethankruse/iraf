"""
Implementation of some utility functions.
"""

import csv
import os
from glob import glob
from typing import List, Union

__all__ = ['is_iterable', 'file_handler']


def is_iterable(obj):
    """
    Return True if the input is iterable, but not a string. Return
    False otherwise.
    """
    from collections.abc import Iterable
    return not isinstance(obj, str) and isinstance(obj, Iterable)


def file_handler(filelist: Union[List[str], str], exists: bool = True,
                 recursive: bool = False) -> List[str]:
    """
    Convert an input string or list to a Python list of filenames.

    Essentially a combination of the functionality of os.path.expanduser, glob,
    and IRAF's @list structure. A modified, limited version of IRAF's imtopen.
    Takes IRAF file inputs and returns a list with actual paths for all
    matching files.

    Allows for a file with a list of file patterns by starting with '@', e.g.
    '@filelist.txt' will recursively find all matches for each line in the
    list, and can handle @lists of @lists. Comment lines (starting with '#')
    are skipped.

    All file names also accept ~ and wildcard (*,?,[a-f]) expansions as well.
    Wildcard and character ranges are handled by `glob.glob`, and the recursive
    option is passed to `glob.glob` for handling of '**'.

    Currently does not support IRAF's image section format.
    e.g. image.fits[1:50, 1:50] will fail rather than allowing only a portion
    of an image to be opened.

    Parameters
    ----------
    filelist : str or list[str]
        String or list of input file strings to expand and create the final
        list of files.
    exists : bool, optional
        Whether files must exist to be returned. If True, will run the pattern
        matching (glob) to look for existing files that match (default).
        If False, only does ~ and @list expansions. Useful for e.g. an
        @output_files.txt list of requested output file names that do not
        yet exist.
    recursive : bool, optional
        Passed to `glob.glob` to control how '**' is handled.

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
        # skip blank lines
        istr = istr.strip()
        if len(istr) == 0:
            continue
        # we have an input file with a list of files to use.
        if istr[0] == '@':
            # remove the @ and deal with home directories
            fnames = os.path.expanduser(istr[1:])
            # deal with wildcards
            if exists:
                fnames = glob(fnames, recursive=recursive)
            else:
                fnames = [fnames]
            # load each requested list
            for fname in fnames:
                with open(fname) as ff:
                    reader = csv.reader(ff)
                    for row in reader:
                        # skip over blank lines and comment lines
                        if (len(row) == 0 or len(row[0].strip()) == 0 or
                                row[0].strip()[0] == '#'):
                            continue
                        # recursively handle these files
                        newfiles = file_handler(row, recursive=recursive,
                                                exists=exists)
                        for newfile in newfiles:
                            outlist.append(newfile)
        # just a single file name, no list
        else:
            fnames = os.path.expanduser(istr)
            # deal with wildcards
            if exists:
                fnames = glob(fnames, recursive=recursive)
            else:
                fnames = [fnames]
            # add the results to the output list
            for fname in fnames:
                outlist.append(fname)

    return outlist
