from iraf.sys import image_open
import re
import copy
import csv

__all__ = ['Instrument', 'get_header_value', 'set_header_value',
           'delete_header_value', 'ccdtypes', 'file_new_copy', 'ccdsubset']


class Instrument(object):
    """
    Instrument object.

    Provides header translations between IRAF standards and the values
    used by specific telescopes or instruments.

    Translation files must be space separated and use single quotes '' to
    mark one field that contains spaces.
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
                           'darkcor': None,
                           'darktime': None, 'datasec': None, 'exptime': None,
                           'fixfile': None, 'fixpix': None, 'flatcor': None,
                           'fringcor': None, 'gain': None, 'illumcor': None,
                           'imagetyp': None, 'mkfringe': None, 'mkillum': None,
                           'ncombine': None, 'nscanrow': None, 'overscan': None,
                           'rdnoise': None,
                           'readcor': None, 'snoise': None, 'subset': None,
                           'trim': None, 'trimsec': None, 'zerocor': None,
                           'origin': None, 'date': None, 'iraf-tlm': None}
        # Each of these parameters can have a default value, but we
        # start off with all None.
        self.defaults = copy.deepcopy(self.parameters)

        imtyps = "object|zero|dark|flat|illum|fringe|other|comp"
        self._default_imagetyps = imtyps.split('|')
        # any custom image type conversions go here
        self.image_types = {}

        # store the reference to the file
        self.translation_file = translation_file
        # create the instrument translations, if given
        if translation_file is not None:
            with open(translation_file, 'r') as ff:
                # allow for multiple spaces between fields for pretty alignment
                # and use single quotes to mark a field with spaces.
                reader = csv.reader(ff, delimiter=' ', skipinitialspace=True,
                                    quotechar="'")
                for row in reader:
                    # ignore blank or commented lines
                    if len(row) == 0 or row[0][0] == '#':
                        continue
                    # translations can only be 2 columns or 3 with a default
                    # value
                    if len(row) < 2 or len(row) > 3:
                        raise Exception(f"Row '{row}' in instrument "
                                        f"translation file {translation_file} "
                                        f"does not have 2 or 3 values.")
                    # change the value of this parameter and its defaults
                    # if present
                    if row[0].lower() in self.parameters:
                        self.parameters[row[0].lower()] = row[1]
                        if len(row) == 3:
                            self.defaults[row[0].lower()] = row[2]
                    # see if this is a custom image type that needs to map
                    # to one of our recognized types
                    elif row[1].lower() in self._default_imagetyps:
                        self.image_types[row[0]] = row[1].lower()
                    # not sure what to do with this row
                    else:
                        raise Exception(f"Unrecognized row in "
                                        f"{translation_file}. Not sure what to "
                                        f"do with '{row}'.")

    def translate(self, key):
        # XXX: this bit is just for testing Instrument.
        if key not in self.parameters:
            raise Exception(f"Remove when done testing. Needed parameter "
                            f"{key}, but it's not in the instrument parameters"
                            f" dict.")
        if key in self.parameters and self.parameters[key] is not None:
            return self.parameters[key]
        else:
            return key

    def get_default(self, key):
        if key in self.defaults:
            return self.defaults[key]
        else:
            return None

    def get_image_type(self, key):
        if key is None:
            return key
        key = key.strip()
        if key in self.image_types:
            return self.image_types[key]
        else:
            return key


def set_header_value(hdulist, instrument, key, value, comment=None):
    """
    Add a header key/value/comment to an image, or modify the existing value
    if the key already exists.

    If comment is None, the existing comment will remain. To clear
    the comment, set comment to the empty string ''.

    The equivalent of IRAF's hdmput*. (e.g. hdmputi, hdmputr)

    Parameters
    ----------
    hdulist
    instrument : Instrument
    key : str
    value
    comment : str

    """
    # what is the header key in the instrument's language
    key = instrument.translate(key)

    # go through the list of hdus and put it in the first one we find
    found = False
    for hdu in hdulist:
        if key in hdu.header:
            hdu.header.set(key, value, comment)
            found = True
            break
    # if it's not in any headers, put it in the first one
    if not found:
        hdulist[0].header.set(key, value, comment)
    return


def get_header_value(hdulist, instrument, key, default=False):
    """
    Retrieve the value of a header key.

    The equivalent of IRAF's hdmg*. (e.g. hdmgstr) when default is False.
    When default is True, it returns the instrument default value for the
    key, equivalent to IRAF's hdmgdef.

    Parameters
    ----------
    hdulist
    instrument
    key
    default

    Returns
    -------
    value

    """
    # what is the header key in the instrument's language
    key = instrument.translate(key)
    val = None

    if not default:
        # go reverse order because we'll typically want the first instance
        # if something is in multiple headers
        for hdu in hdulist[::-1]:
            if key in hdu.header:
                val = hdu.header[key]
    if val is None:
        val = instrument.get_default(key)

    return val


def delete_header_value(hdulist, instrument, key):
    """
    Delete a header field.

    The equivalent of IRAF's hdmdelf.

    Parameters
    ----------
    hdulist
    instrument
    key

    """
    # what is the header key in the instrument's language
    key = instrument.translate(key)

    # go reverse order because we'll typically want the first instance
    # if something is in multiple headers
    for hdu in hdulist:
        if key in hdu.header:
            hdu.header.remove(key, remove_all=True)
    return


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
    hdulist
    instrument : Instrument

    Returns
    -------

    """

    if not isinstance(instrument, Instrument):
        raise Exception('ccdtypes not given an Instrument object.')

    options = "object|zero|dark|flat|illum|fringe|other|comp".split('|')
    typ = get_header_value(hdulist, instrument, 'imagetyp')
    typ = instrument.get_image_type(typ)

    if typ is None:
        typ = 'none'
    elif typ not in options:
        typ = 'unknown'

    return typ


def ccdsubset(hdulist, instrument):
    """
    Get the subset identifier as the header value of 'subset' (or instrument
    equivalent).

    Parameters
    ----------
    hdulist
    instrument

    Returns
    -------

    """

    if not isinstance(instrument, Instrument):
        raise Exception('ccdsubset not given an Instrument object.')

    subsetstr = get_header_value(hdulist, instrument, 'subset')

    if subsetstr is None:
        subsetstr = ''

    subsetstr = subsetstr.strip()

    # Replace non alphanumeric or '.' characters by '_'
    # since the subset ID is used in forming image names.
    subsetstr = re.sub(r'[^\w.]', '_', subsetstr)

    """
    # This bit was a translation of the original IRAF ccdsubset function
    # that uses the subsets file and turns things into shorter subset strings
    # while also making sure there aren't overlaps. This all seems unnecessary
    # and we should just use the full subset string as the identifier. No
    # point in maintaining some subsets file.

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


def file_new_copy(outstr, in_header, mode='NEW_COPY', overwrite=True,
                  instrument=None):
    """
    Copy a file, adding appropriate header info about IRAF sourcing.

    Meant to be equivalent to immap (output, NEW_COPY, Memi[in]) in IRAF.

    Parameters
    ----------
    outstr
    in_header
    mode
    overwrite
    instrument

    Returns
    -------

    """
    if mode != 'NEW_COPY':
        print('other modes of file_new_copy not supported.')
        return

    ftype = in_header.__filetype__
    if ftype == 'fits':
        # do the actual writing of the copy
        in_header.writeto(outstr, overwrite=overwrite)
        if instrument is None:
            instrument = Instrument()
        # update header parameters in the new file without altering the old
        hdulist = image_open(outstr, mode='update')
        # XXX: update this with a real package name at some point
        orval = 'AIRAF in Python'
        orcomm = 'FITS file originator'
        set_header_value(hdulist, instrument, 'origin', orval, orcomm)
        import datetime
        # don't need microsecond precision in our dates, so leave it out
        dtval = datetime.datetime.now().replace(microsecond=0).isoformat()
        dtcomm = 'Date FITS file was generated'
        set_header_value(hdulist, instrument, 'date', dtval, dtcomm)
        lmcomm = 'Time of last modification'
        set_header_value(hdulist, instrument, 'iraf-tlm', dtval, lmcomm)
        hdulist.close()
    else:
        raise Exception(f'file_new_copy of file type {ftype} not yet '
                        f'implemented.')
    return
