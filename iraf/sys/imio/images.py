from astropy.io import fits

__all__ = ['image_open']


def image_open(image, mode='readonly'):
    """
    Open an image and return an AIRAF image object.

    Supposed to be a function that handles images of all types
    and returns data/headers etc in one consistent format for use in
    AIRAF. Currently only supports FITS and simply returns a standard
    HDUList with a flag indicating it was opened here.

    Always returns a context managers, so can be used in 'with ... as:'
    blocks.

    Parameters
    ----------
    image : str
        Image to be opened.
    mode : {'readonly', 'update', 'append'}
        Mode the file is to be opened with.

    Returns
    -------
    HDUList object

    """
    hdulist = None
    ftype = ''
    err = None
    # start with FITS files by default
    try:
        hdulist = fits.open(image, mode=mode)
        ftype = 'fits'
    # catch a FITS problem and move on to other file types
    except IOError as ierr:
        err = ierr

    # if we've tried all file types and found one, use it;
    # otherwise raise the most recent error message. (?)
    if hdulist is not None:
        hdulist.__filetype__ = ftype
    else:
        raise err

    return hdulist
