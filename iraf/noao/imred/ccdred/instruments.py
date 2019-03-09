from iraf.sys import image_open
from ..ccdred import imagetypes
import re
import copy
import shlex
import os
from datetime import datetime

__all__ = ['Instrument', 'get_header_value', 'set_header_value',
           'delete_header_value', 'ccdtypes', 'file_new_copy', 'ccdsubset']


class Instrument(object):
    """
    Provide translations between IRAF standards and parameter names or image
    types used by specific telescopes or instruments.

    Internally, IRAF has a standard set of parameter names, but these may not
    be followed by all observatories. For example, IRAF always assumes the
    total exposure time for an image is found in the header with keyword
    'exptime'. If your image uses 'texpose' instead, the Instrument object
    can provide the necessary translation. Similarly, an image labeled as
    type 'sky flat' must be translated to the IRAF standard image type 'flat'.

    Translation files must be white space separated and use quotes to
    mark one field that contains spaces.

    Parameters
    ----------
    translation_file : str, optional
        File with two or three whitespace separated values. For parameter names,
        the first value is the IRAF standard parameter name, the second value
        is the custom instrument specific header key, and the optional third
        value is a default value if the parameter is not found in the header.
        Ex:
        exptime customexp
        biassec  custombias    [411:431,2:573]

        For image types, the first value is the custom image type and the
        second value is the IRAF default.
        Wrap any values with spaces in quotes, e.g.:
        "sky flat" flat

    Attributes
    ----------
    translation_file : str
        A reference to the translation file used to initialize the object.
    parameters : dict
        Each internal IRAF parameter name is a key in this dictionary, and the
        values are the custom translations used (None if using the default).
    defaults : dict
        Same keys as parameters, this time listing the default values for each
        IRAF parameter
    image_types : dict
        Any custom image types and their translations to the IRAF standard.

    Raises
    ------
    Exception
        When a line in the translation file cannot be interpreted.

    """
    """
    To learn about instruments, look into
    'iraf-src/noao/imred/ccdred/src/hdrmap.x' for the functions that create
    a symbol table pointer and do the translations.
    Instrument file tables are located in
    'iraf-src/noao/imred/ccdred/ccddb'.
    """

    def __init__(self, translation_file=None):
        # this is the full list of parameters ccdred references in image
        # headers. Anything not None says to use that value as the header
        # keyword instead of this default parameter name.
        self.parameters = {'BPM': None, 'biassec': None, 'ccdmean': None,
                           'ccdmeant': None, 'ccdproc': None, 'ccdsec': None,
                           'darkcor': None, 'darktime': None, 'datasec': None,
                           'exptime': None, 'fixfile': None, 'fixpix': None,
                           'flatcor': None, 'fringcor': None, 'gain': None,
                           'illumcor': None, 'imagetyp': None, 'mkfringe': None,
                           'mkillum': None, 'ncombine': None, 'nscanrow': None,
                           'overscan': None, 'rdnoise': None, 'readcor': None,
                           'snoise': None, 'subset': None, 'trim': None,
                           'trimsec': None, 'zerocor': None, 'origin': None,
                           'date': None, 'iraf-tlm': None}
        # Each of these parameters can have a default value, but we
        # start off with all None.
        self.defaults = copy.deepcopy(self.parameters)

        self._default_imagetyps = imagetypes
        # any custom image type conversions go here
        self.image_types = {}

        # store the reference to the file
        self.translation_file = translation_file
        # create the instrument translations, if given
        if translation_file is not None:
            with open(os.path.expanduser(translation_file), 'r') as ff:
                # allow for multiple spaces between fields for pretty alignment
                # and use quotes to mark a field with spaces.
                for row in ff.readlines():
                    row = shlex.split(row)
                    # allow comments at the end of the line, but ignore them
                    for ii in range(len(row)):
                        if row[ii][0] == '#':
                            row = row[:ii]
                            break
                    # ignore blank or fully commented lines
                    if len(row) == 0:
                        continue
                    # translations can only be 2 columns or 3 with a default
                    # value
                    if len(row) < 2 or len(row) > 3:
                        raise Exception(f"Row '{row}' in instrument "
                                        f"translation file {translation_file} "
                                        f"does not have 2 or 3 values.")
                    # change the value of this parameter and its defaults
                    # if present
                    if row[0] in self.parameters:
                        self.parameters[row[0]] = row[1]
                        if len(row) == 3:
                            self.defaults[row[0]] = row[2]
                    # see if this is a custom image type that needs to map
                    # to one of our recognized types
                    elif row[1] in self._default_imagetyps and len(row) == 2:
                        self.image_types[row[0]] = row[1]
                    # not sure what to do with this row
                    else:
                        raise Exception(f"Unrecognized row in translation file "
                                        f"{translation_file}.\nNot sure what "
                                        f"to do with '{row}'.")

    def translate(self, key):
        """
        Return the instrument specific translation of the IRAF parameter 'key'.
        """
        if self.parameters[key] is not None:
            return self.parameters[key]
        else:
            return key

    def get_default(self, key):
        """
        Return the default value of the IRAF parameter 'key'.
        """
        return self.defaults[key]

    def get_image_type(self, key):
        """
        Translate an input image type into one of the 10 possible IRAF image
        types: "|object|zero|dark|flat|illum|fringe|other|comp|none|unknown|".

        Parameters
        ----------
        key : str or None

        Returns
        -------
        str, one of the 10 possible image types:
        "|object|zero|dark|flat|illum|fringe|other|comp|none|unknown|"
        """
        if key is None:
            return 'none'
        key = key.strip()
        if len(key) == 0:
            return 'none'
        if key in self.image_types:
            key = self.image_types[key]
        # make sure no one messed with self.image_types to give a nonvalid
        # image type
        if key in imagetypes:
            return key
        else:
            return 'unknown'


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
    if not isinstance(instrument, Instrument):
        raise Exception('ccdtypes not given an Instrument object.')

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

    if not isinstance(instrument, Instrument):
        raise Exception('ccdsubset not given an Instrument object.')

    subsetstr = get_header_value(hdulist, instrument, 'subset')

    if subsetstr is None:
        subsetstr = ''

    subsetstr = subsetstr.strip()

    # Replace any non alphanumeric or '.' or '_' characters by '_'
    # since the subset ID is used in forming image names.
    subsetstr = re.sub(r'[^\w.]', '_', subsetstr)

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
        raise Exception(f'Mode "{mode}" of file_new_copy not supported.')

    outpath = os.path.expanduser(outpath)
    # NOTE: IRAF appears to put these three header values (origin, date,
    # iraf-tlm) at the beginning of the header, right after the 'extend'
    # value. We currently append at the end. Not sure if this matters for
    # anyone.
    ftype = in_header.__filetype__
    if ftype == 'fits':
        # do the actual writing of the copy
        in_header.writeto(outpath, overwrite=overwrite)
        if instrument is None:
            instrument = Instrument()
        # update header parameters in the new file without altering the old
        hdulist = image_open(outpath, mode='update')
        orcomm = 'FITS file originator'
        from iraf import __hdrstring__ as fitsname
        set_header_value(hdulist, instrument, 'origin', fitsname, orcomm)
        # don't need microsecond precision in our dates, so leave it out
        dtval = datetime.utcnow().replace(microsecond=0).isoformat()
        dtcomm = 'Date FITS file was generated (UTC)'
        set_header_value(hdulist, instrument, 'date', dtval, dtcomm)
        lmcomm = 'Time of last modification (UTC)'
        set_header_value(hdulist, instrument, 'iraf-tlm', dtval, lmcomm)
        hdulist.close()
    else:
        raise Exception(f'file_new_copy of file type {ftype} not yet '
                        f'implemented.')
    return
