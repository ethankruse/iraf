import os
import numpy as np
from astropy.io import fits
import iraf
import pytest


def test_image_open_close_fits(tmpdir):
    basedir = str(tmpdir)
    inim = os.path.join(basedir, f'testimg.fits')

    # test the file doesn't exist and this raises an error
    assert not os.path.exists(inim)
    with pytest.raises(OSError):
        iraf.sys.image_open(inim)

    # test opening a FITS file
    sz = 5
    arr = np.zeros(sz)
    hdu = fits.PrimaryHDU(arr)
    hdu.writeto(inim, overwrite=True)
    im = iraf.sys.image_open(inim)
    assert im is not None
    assert im.__filetype__ == 'fits'

    # test image closing
    # make sure the file is still open
    assert im[0].data.size == sz
    im.close()
    # file should now be closed and data inaccessible
    with pytest.raises(ValueError):
        print(im[0].data.size)

    # test context manager
    with iraf.sys.image_open(inim) as im:
        # make sure the file is still open
        assert im[0].data.size == sz
    # file should now be closed and data inaccessible
    with pytest.raises(ValueError):
        print(im[0].data.size)

    # make sure you can't change values in 'readonly' mode
    with iraf.sys.image_open(inim, mode='readonly') as im:
        im[0].data += 1
    with iraf.sys.image_open(inim, mode='readonly') as im:
        assert np.allclose(im[0].data, 0)
    # make sure update mode works
    with iraf.sys.image_open(inim, mode='update') as im:
        im[0].data += 1
    with iraf.sys.image_open(inim, mode='readonly') as im:
        assert np.allclose(im[0].data, 1)

    # test bad modes
    with pytest.raises(ValueError):
        iraf.sys.image_open(inim, mode='foo')

    # test a text file that shouldn't be able to be opened
    with open(inim, 'w') as ff:
        ff.write('blah\n')
    with pytest.raises(OSError):
        iraf.sys.image_open(inim)
