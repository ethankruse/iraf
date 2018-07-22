import os
import numpy as np
from astropy.io import fits
import iraf
import pytest


def test_image_open_close(tmpdir):
    basedir = str(tmpdir)
    inim = os.path.join(basedir, f'testimg.fits')

    # test the file doesn't exist and this raises an error
    assert not os.path.exists(inim)
    with pytest.raises(IOError):
        iraf.sys.image_open(inim)

    # test opening a fits file
    sz = 5
    arr = np.zeros(sz)
    hdu = fits.PrimaryHDU(arr)
    hdu.writeto(inim, overwrite=True)
    im = iraf.sys.image_open(inim)
    assert im is not None
    assert im.__filetype__ == 'fits'

    # test image closing
    im.__filetype__ = 'bs'
    # unrecognized file types raise errors
    with pytest.raises(Exception):
        iraf.sys.image_close(im)
    im.__filetype__ = 'fits'
    # make sure the file is still open
    assert im[0].data.size == sz
    iraf.sys.image_close(im)
    # file should now be closed and data inaccessible
    with pytest.raises(ValueError):
        print(im[0].data.size)

    # test context manager
