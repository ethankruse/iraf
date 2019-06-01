import copy
import shlex
import os
from ..ccdred import _imagetypes

__all__ = ['Instrument']


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
                           'flatcor': None, 'fringcor': None, 'fringscl': None,
                           'gain': None, 'illumcor': None, 'imagetyp': None,
                           'mkfringe': None, 'mkillum': None, 'ncombine': None,
                           'nscanrow': None, 'overscan': None, 'rdnoise': None,
                           'readcor': None, 'scancor': None, 'snoise': None,
                           'subset': None, 'trim': None, 'trimsec': None,
                           'zerocor': None, 'origin': None, 'date': None,
                           'iraf-tlm': None}
        # Each of these parameters can have a default value, but we
        # start off with all None.
        self.defaults = copy.deepcopy(self.parameters)

        self._default_imagetyps = _imagetypes
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
        if key in _imagetypes:
            return key
        else:
            return 'unknown'
