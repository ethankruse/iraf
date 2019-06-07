"""
Utility functions needed by the main tasks in ccdred, but not intended to be
used directly by the user.
"""

# keep the namespace clean
import re as _re
import os as _os
import numpy as _np
from datetime import datetime as _datetime
from iraf.sys import image_open as _image_open
from . import Instrument as _Instrument

__all__ = ['ccdtypes', 'set_header_value', 'get_header_value', 'CCDProcError',
           'delete_header_value', 'ccdsubset', 'file_new_copy', 'type_max']


class CCDProcError(Exception):
    pass


def ccdtypes(hdulist, instrument):
    """
    Get the header value of 'imagetyp' (or instrument equivalent).

    If that (instrument converted) value is one of
    "object|zero|dark|flat|illum|fringe|other|comp" return that.
    Otherwise return 'unknown' unless the header does not contain any
    'imagetyp' and the instrument doesn't have a default value for it. Then
    return 'none'.

    Parameters
    ----------
    hdulist : IRAF image
        An open IRAF image we want the ccd type of.
    instrument : Instrument
        The instrument specific translations.

    Returns
    -------
    str, one of the 10 possible image types:
    "|object|zero|dark|flat|illum|fringe|other|comp|none|unknown|"
    """
    if not isinstance(instrument, _Instrument):
        raise TypeError('ccdtypes not given an Instrument object.')

    typ = get_header_value(hdulist, instrument, 'imagetyp')
    return instrument.get_image_type(typ)


def set_header_value(hdulist, instrument, key, value, comment=None):
    """
    Add a header key/value/comment to an image, or modify the existing value
    if the key already exists. First appropriately translate the IRAF default
    parameter 'key' into the instrument equivalent.

    If comment is None, the existing comment will remain. To clear
    the comment, set comment to the empty string ''.

    The equivalent of IRAF's hdmput*. (e.g. hdmputi, hdmputr)

    Parameters
    ----------
    hdulist : IRAF image
    instrument : Instrument
    key : str
    value
    comment : str, optional
    """
    # what is the header key in the instrument's language
    key = instrument.translate(key)

    # for some reason comment='' fails in pytest, but ' ' works and enters
    # an empty comment.
    if comment is not None and len(comment) == 0:
        comment = ' '

    # go through the list of hdus and put it in the first one we find
    found = False
    for hdu in hdulist:
        if key in hdu.header:
            hdu.header.set(key, value=value, comment=comment)
            found = True
            break
    # if it's not in any headers, put it in the first one
    # XXX: is this what IRAF does? what if the same key is in two headers?
    if not found:
        hdulist[0].header.set(key, value=value, comment=comment)
    return


def get_header_value(hdulist, instrument, key, default=False):
    """
    Retrieve the value of the IRAF parameter 'key' from the instrument
    specific header keyword.

    The equivalent of IRAF's hdmg*. (e.g. hdmgstr) when default is False.
    When default is True, it returns the instrument default value for the
    key, equivalent to IRAF's hdmgdef.

    NOTE: in IRAF there's a distinction between int/real and strings.
    hdmgstr/hdmgdef returns the empty string while hdmgeti/hdmgetr
    raises an error if there is no value or default.
    Here we always just return None and expect that to be handled by the
    calling function depending on what return type it is expecting.

    Parameters
    ----------
    hdulist : IRAF image
    instrument : Instrument
    key : str
    default : bool, optional

    Returns
    -------
    value

    """
    # what is the header key in the instrument's language
    userkey = instrument.translate(key)
    val = None

    if not default:
        # grab the first instance
        # XXX: is this what IRAF does? what if the same key is in two headers?
        for hdu in hdulist:
            if userkey in hdu.header:
                val = hdu.header[userkey]
                break
    if val is None:
        val = instrument.get_default(key)

    return val


def delete_header_value(hdulist, instrument, key):
    """
    Delete the translated IRAF parameter 'key' from the header.

    The equivalent of IRAF's hdmdelf.

    Parameters
    ----------
    hdulist : IRAF image
    instrument : Instrument
    key : str

    """
    # what is the header key in the instrument's language
    key = instrument.translate(key)

    for hdu in hdulist:
        # XXX: does IRAF really remove from all headers?
        if key in hdu.header:
            hdu.header.remove(key, remove_all=True)
    return


def ccdsubset(hdulist, instrument):
    """
    Get the subset identifier as the header value of 'subset' (or instrument
    equivalent). Replace any characters not suitable for file names with '_'.
    If no subset is found, return the empty string ''.

    Parameters
    ----------
    hdulist : IRAF image
    instrument : Instrument

    Returns
    -------
    str
    """

    if not isinstance(instrument, _Instrument):
        raise TypeError('ccdsubset not given an Instrument object.')

    subsetstr = get_header_value(hdulist, instrument, 'subset')

    if subsetstr is None:
        subsetstr = ''

    subsetstr = subsetstr.strip()

    # Replace any non alphanumeric or '.' or '_' characters by '_'
    # since the subset ID is used in forming image names.
    subsetstr = _re.sub(r'[^\w.]', '_', subsetstr)

    """
    # This bit was a translation of the original IRAF ccdsubset function
    # that uses the subsets file and turns things into shorter subset strings
    # while also making sure there aren't overlaps. This all seems unnecessary
    # and we should just use the full subset string as the identifier. No
    # point in maintaining subsets files.

    import shlex

    if ssfile is None:
        print('ssfile must be defined to use subsets')
        sys.exit(1)

    # The default subset identifier is the first
    # word/string of the subset string.
    subset1 = shlex.split(subsetstr.strip())
    if len(subset1) > 0:
        subset1 = subset1[0]
    else:
        subset1 = None

    # A null subset string is ok.  If not null check for conflict
    # with previous subset IDs.

    if subset1 is not None:
        orig = subset1
        # whether or not to append this to the subsets file
        append = True
        # Search the subset record file for the same subset string.
        # If found use the ID string.  If the subset ID has been
        # used for another subset string then increment an integer
        # suffix to the default ID and check the list again.
        if os.path.exists(ssfile):
            subsetstrs = []
            subsetids = []
            with open(ssfile, 'r') as ff:
                lines = ff.readlines()
                for row in lines:
                    # skip over blank lines and comment lines
                    if (len(row) == 0 or len(row.strip()) == 0 or
                            row[0].strip()[0] == '#'):
                        continue

                    groups = shlex.strip(row)
                    assert len(groups) == 2
                    subsetstrs.append(groups[0])
                    subsetids.append(groups[1])

            ii = 1
            while True:
                if subsetstr in subsetstrs:
                    append = False
                    subset1 = subsetids[subsetstrs.index(subsetstr)]
                    break
                else:
                    if subset1 in subsetids:
                        subset1 = '{0}{1:d}'.format(orig, ii)
                        ii += 1
                    else:
                        break

        if append:
            with open(ssfile, 'a') as ff:
                ff.write('{0}\t{1}\n'.format(subsetstr, subset1))

        # Set the subset ID string and replace magic characters by '_'
        # since the subset ID is used in forming image names.
        for ii, ichar in enumerate(subset1):
            if not ichar.isalnum() and ichar != '.':
                subset1[ii] = '_'

    return subset1
    """
    return subsetstr


def file_new_copy(outpath, in_header, mode='NEW_COPY', overwrite=True,
                  instrument=None):
    """
    Copy a file, adding appropriate header info about IRAF sourcing.

    Meant to be equivalent to immap (output, NEW_COPY, Memi[in]) in IRAF.

    Parameters
    ----------
    outpath : str
        Path to new file. Performs ~ expansions first.
    in_header : IRAF image
        Open image to be copied
    mode : {'NEW_COPY'}
        What mode to copy the file
    overwrite : bool, optional
        If True, overwrites any previously existing file. If False, raises an
        error if the output file already exists.
    instrument : Instrument, optional
        Translations for the header values added as part of the copy.
    """
    if mode != 'NEW_COPY':
        raise ValueError(f'Mode "{mode}" of file_new_copy not supported.')

    outpath = _os.path.expanduser(outpath)
    # NOTE: IRAF appears to put these three header values (origin, date,
    # iraf-tlm) at the beginning of the header, right after the 'extend'
    # value. We currently append at the end. Not sure if this matters for
    # anyone.
    ftype = in_header.__filetype__
    if ftype == 'fits':
        # do the actual writing of the copy
        in_header.writeto(outpath, overwrite=overwrite)
        if instrument is None:
            instrument = _Instrument()
        # update header parameters in the new file without altering the old
        hdulist = _image_open(outpath, mode='update')
        orcomm = 'FITS file originator'
        from iraf import __hdrstring__ as fitsname
        set_header_value(hdulist, instrument, 'origin', fitsname, orcomm)
        # don't need microsecond precision in our dates, so leave it out
        dtval = _datetime.utcnow().replace(microsecond=0).isoformat()
        dtcomm = 'Date FITS file was generated (UTC)'
        set_header_value(hdulist, instrument, 'date', dtval, dtcomm)
        lmcomm = 'Time of last modification (UTC)'
        set_header_value(hdulist, instrument, 'iraf-tlm', dtval, lmcomm)
        hdulist.close()
    else:
        raise NotImplementedError(f'file_new_copy of file type {ftype} not yet '
                                  f'implemented.')
    return


def type_max(type1, type2):
    """
    Return the datatype of highest precedence of the two input data types.

    Input two data types, and return the data type of highest precedence, such
    that the lower precedence type can be safely converted to the higher type.

    For example, inputs of np.float64 and np.float32 will return np.float64
    because all float32s can be safely converted to float64s without losing
    any information, but the reverse is not true.

    Parameters
    ----------
    type1, type2 : data type

    Returns
    -------
    data type

    Raises
    ------
    TypeError
        When the two inputs cannot be understood or safely converted in either
        direction.
    """
    right = _np.can_cast(type1, type2, casting='safe')
    left = _np.can_cast(type2, type1, casting='safe')

    if left:
        return type1
    if right:
        return type2

    """
    # if we're here, it's likely case of an unsigned int and signed int
    # of same size
    ints = [np.int8, np.int16, np.int32, np.int64]
    if (np.issubdtype(type1.type, np.unsignedinteger) and
            np.issubdtype(type2.type, np.integer)):
        for iint in ints:
            if np.can_cast(type1, iint, casting='safe'):
                return np.dtype(iint)

    elif (np.issubdtype(type2.type, np.unsignedinteger) and
              np.issubdtype(type1.type, np.integer)):
        for iint in ints:
            if np.can_cast(type2, iint, casting='safe'):
                return np.dtype(iint)
    """
    errstr = "Unrecognized dtype or cannot safely cast between {0} and {1}."
    raise TypeError(errstr.format(type1, type2))
