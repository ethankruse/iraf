from astropy.io import fits

__all__ = ['image_open', 'image_close']


def image_open(image, mode='readonly'):
    """
    Open an image and return an AIRAF image object.

    Supposed to be a function that handles opening images and returning
    data/headers etc in one consistent format. Currently only supports FITS
    and simply returns a standard HDUList with a flag indicating it was
    opened here.

    Parameters
    ----------
    image : str
        Image to be opened.
    mode : str
        Mode the file is to be opened with.

    Returns
    -------
    HDUList object
    """
    hdulist = None
    ftype = ''
    try:
        hdulist = fits.open(image, mode=mode)
        ftype = 'fits'
    except IOError:
        print(f"Error reading image {image} ...")

    if hdulist is not None:
        hdulist.__filetype__ = ftype

    return hdulist


def image_close(image):
    """
    Close an AIRAF image.

    Supposed to be a function that handles closing images of all types.
    Currently just operates as image.close() though, assuming it was opened
    with `image_open`.

    Parameters
    ----------
    image : AIRAF image object

    Returns
    -------

    """
    if image.__filetype__ == 'fits':
        image.close()
    else:
        raise Exception(f"{image} of unknown type {image.__filetype__}")
