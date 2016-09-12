from astropy.io import fits

__all__ = ['image_open', 'image_close']


def image_open(image):
    hdulist = None
    try:
        hdulist = fits.open(image)
    except IOError:
        print("Error reading image {0} ...".format(image))

    return hdulist


def image_close(image):
    image.close()
