import pytest
import numpy as np
from astropy.io import fits
import iraf
import os


def test_combine_basic(tmpdir):
    # create some simple files for testing
    nimg = 5
    nx = 20
    ny = 30
    inputs = []
    basedir = str(tmpdir.mkdir("combine"))
    arr = None
    for ii in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float)
        hdu = fits.PrimaryHDU(arr)
        inim = os.path.join(basedir, f'testimg{ii:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)

    # for the IRAF list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')

    outfile = os.path.join(basedir, 'testout_me.fits')
    iraflist = '@' + inlist

    # explicitly call combine with every parameter set to what we want
    iraf.combine(iraflist, outfile, plfile=None, sigma=None, ccdtype=None,
                 subsets=False, delete=False, method='average',
                 reject='none', project=False, outtype=None, offsets='none',
                 masktype='none', maskvalue=0., blank=0., scale=None, zero=None,
                 weight=None, statsec=None, lthreshold=None, hthreshold=None,
                 nlow=1, nhigh=1, nkeep=1, mclip=True, lsigma=3.0,
                 hsigma=3.0, rdnoise=0., gain=1., snoise=0., sigscale=0.1,
                 pclip=-0.5, grow=0, instrument=None, logfile=None,
                 verbose=False,
                 ssfile=None)

    image_me = fits.open(outfile)
    assert image_me[0].data.shape == (nx, ny)
    assert np.allclose(image_me[0].data, arr)