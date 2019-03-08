from astropy.io import fits

__all__ = ['image_open']


def image_open(image, mode='readonly'):
    """
    Open an image and return an AIRAF image object.

    Supposed to be a function that handles images of all types
    and returns data/headers etc in one consistent format for use in
    AIRAF. Currently only supports FITS and simply returns a standard
    HDUList with a flag indicating it was opened here.

    Always returns a context manager, so can be used in 'with ... as:'
    blocks.

    Will always return an image object. Raises an Exception if the
    image cannot be opened.

    Parameters
    ----------
    image : str
        Image to be opened.
    mode : {'readonly', 'update'}
        Mode the file is to be opened with.

    Returns
    -------
    HDUList object

    Raises
    ------
    Exception
        If the input file cannot be opened or understood.

    """
    hdulist = None
    ftype = ''
    # start with FITS files by default
    try:
        hdulist = fits.open(image, mode=mode)
        ftype = 'fits'
    # catch a FITS problem and move on to other file types
    except IOError:
        pass

    # if we've tried all file types and found one, use it;
    # otherwise raise an error
    if hdulist is not None:
        hdulist.__filetype__ = ftype
    else:
        raise Exception(f'Unable to open image {image}')

    return hdulist
