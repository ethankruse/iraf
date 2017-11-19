from astropy.io import fits

__all__ = ['image_open', 'image_close']


def image_open(image):
    """
    Supposed to be a function that handles opening images and returning
    data/headers etc in one consistent format. Currently only supports FITS
    though.

    Parameters
    ----------
    image

    Returns
    -------

    """
    hdulist = None
    ftype = ''
    try:
        hdulist = fits.open(image)
        ftype = 'fits'
    except IOError:
        print("Error reading image {0} ...".format(image))

    if hdulist is not None:
        hdulist.__filetype__ = ftype

    return hdulist


def image_close(image):
    """
    Supposed to be a function that handles closing images of all types.
    Currently just operates as image.close() though.

    Parameters
    ----------
    image

    Returns
    -------

    """
    if image.__filetype__ == 'fits':
        image.close()
