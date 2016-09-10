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
