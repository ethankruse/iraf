import numpy as np
from astropy.io import fits
import iraf
import os
import copy

# explicitly call combine with every parameter set to what we want
defaultargs = {'plfile': None, 'sigmafile': None, 'ccdtype': None,
               'subsets': False, 'delete': False, 'method': 'average',
               'reject': 'none', 'project': False, 'outtype': None,
               'offsets': 'none', 'masktype': 'none', 'maskvalue': 0.,
               'blank': 0., 'scale': None, 'zero': None, 'weight': None,
               'statsec': None, 'lthreshold': None, 'hthreshold': None,
               'nlow': 1, 'nhigh': 1, 'nkeep': 1, 'mclip': True,
               'lsigma': 3.0, 'hsigma': 3.0, 'rdnoise': 0., 'gain': 1.,
               'snoise': 0., 'sigscale': 0.1, 'pclip': -0.5, 'grow': 0,
               'instrument': None, 'logfile': None, 'verbose': False,
               'ssfile': None}

dtypes = []


def simple_inputs(nimg, nx, ny, basedir):
    inputs = []
    for ii in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float)
        hdu = fits.PrimaryHDU(arr)
        inim = os.path.join(basedir, f'testimg{ii:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)
    return inputs


# combine_dir is set up in conftest.py
def test_basic(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    nx = 20
    ny = 30
    inputs = simple_inputs(5, nx, ny, basedir)

    # for the IRAF list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')
    iraflist = '@' + inlist

    outfile = os.path.join(basedir, 'testout_me.fits')

    iraf.combine(iraflist, outfile, **defaultargs)

    outim = fits.open(outfile)
    assert outim[0].data.shape == (nx, ny)
    # need to check name instead of dtype because of endianness problems
    # np.dtype('>f4') != np.dtype('<f4') but both have .name == 'float32'
    # assert arr.dtype.name == outim[0].data.dtype.name

    outim.close()


def test_reject_none(combine_dir):
    pass


def test_reject_minmax(combine_dir):
    pass


def test_plfile(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    nx = 20
    ny = 30
    inputs = simple_inputs(5, nx, ny, basedir)

    # for the IRAF list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')
    iraflist = '@' + inlist

    outfile = os.path.join(basedir, 'testout_me.fits')
    plfile = os.path.join(basedir, 'testout_pl.fits')
    myargs = copy.deepcopy(defaultargs)
    myargs['plfile'] = plfile

    iraf.combine(iraflist, outfile, **myargs)
    # XXX: need to figure out what this is
    assert 1
