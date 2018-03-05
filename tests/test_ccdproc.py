import iraf
import pytest
import os
from pathlib import Path
import numpy as np
from astropy.io import fits
import copy
import time


def test_ccd_section():
    defaults = [5, 6, 7, 8, 9, 10]
    d1, d2, d3, d4, d5, d6 = defaults
    # get the defaults
    x0, x1, xs, y0, y1, ys = iraf.ccd_section('[:,:]', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    x0, x1, xs, y0, y1, ys = iraf.ccd_section(None, defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    x0, x1, xs, y0, y1, ys = iraf.ccd_section('  ', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    # all the options
    x0, x1, xs, y0, y1, ys = iraf.ccd_section('[1:5:3,5:20:-4]',
                                              defaults=defaults)
    assert (x0 == 1 and x1 == 5 and xs == 3 and y0 == 5 and y1 == 20 and
            ys == -4)

    # no step size
    x0, x1, xs, y0, y1, ys = iraf.ccd_section('[1:5,5:20]', defaults=defaults)
    assert (x0 == 1 and x1 == 5 and xs == d3 and y0 == 5 and y1 == 20 and
            ys == d6)

    # single row
    x0, x1, xs, y0, y1, ys = iraf.ccd_section('[9,5:20]', defaults=defaults)
    assert (x0 == 9 and x1 == 9 and xs == d3 and y0 == 5 and y1 == 20 and
            ys == d6)

    with pytest.raises(Exception):
        # test no braces
        iraf.ccd_section('[:,:')
        iraf.ccd_section(':,:]')
        iraf.ccd_section('1:2,5:34')
        # test wrong dimensions
        iraf.ccd_section('[1:5]')
        iraf.ccd_section('[1:5, 1:5, 1:5]')
        # test bad inputs in a dimension
        iraf.ccd_section('[1:3:4:5, :]')
        # test bad size of defaults
        iraf.ccd_section('[:,:]', defaults=[1, 2, 3])


def test_file_handler(combine_dir):
    basedir = str(combine_dir)
    ltxt = 4
    lshorttxt = 3
    ifiles = ['a.txt', 'b.txt', 'c.txt', 'ff.txt', '1.fits', '2.fits', '3.fits']
    samplefiles = []
    # create the files
    for ifile in ifiles:
        ipath = os.path.join(basedir, ifile)
        samplefiles.append(ipath)
        Path(ipath).touch()

    # test None
    assert len(iraf.utils.file_handler(None)) == 0

    # simple test of a file that doesn't exist
    inp = os.path.join(basedir, 'abc.txt')
    out = iraf.utils.file_handler(inp)
    assert len(out) == 1 and out[0] == inp

    # test of @lists
    atlist = os.path.join(basedir, 'atlist.list')
    with open(atlist, 'w') as ff:
        for sample in samplefiles:
            ff.write(sample + '\n')

    assert len(iraf.utils.file_handler('@'+atlist)) == len(samplefiles)

    # test multiple inputs and wildcards
    wild = [os.path.join(basedir, '*.txt'), os.path.join(basedir, '*.fits')]
    assert len(iraf.utils.file_handler(wild[0])) == ltxt
    assert len(iraf.utils.file_handler(wild)) == len(samplefiles)

    wild = os.path.join(basedir, '?.txt')
    assert len(iraf.utils.file_handler(wild)) == lshorttxt

    # make sure home directory expansion works
    assert len(iraf.utils.file_handler('~')[0]) > 1

    # check wildcards in @list names
    atwild = os.path.join(basedir, '*list')
    assert len(iraf.utils.file_handler('@'+atwild)) == len(samplefiles)

    # check wildcards in @list entries
    wild = [os.path.join(basedir, '*.txt'), os.path.join(basedir, '*.fits')]
    with open(atlist, 'w') as ff:
        for iwild in wild:
            ff.write(iwild + '\n')

    assert len(iraf.utils.file_handler('@' + atlist)) == len(samplefiles)

    # check comments and blank lines in @lists
    with open(atlist, 'w') as ff:
        for sample in samplefiles:
            ff.write('#' + sample + '\n\n\n')

    assert len(iraf.utils.file_handler('@' + atlist)) == 0

    # check @lists of @lists
    txtlist = os.path.join(basedir, 'txtlist.list')
    nwrite = 2
    with open(txtlist, 'w') as ff:
        ff.write('hi.txt\n\n')
        ff.write('hi2.txt')

    reallist = os.path.join(basedir, 'reallist.list')
    wild = [os.path.join(basedir, '*.txt'), os.path.join(basedir, '*.fits')]
    with open(reallist, 'w') as ff:
        for iwild in wild:
            ff.write(iwild + '\n')

    with open(atlist, 'w') as ff:
        ff.write('\t@' + txtlist + '\n')
        ff.write('   @' + reallist + '\n')

    tot = nwrite + len(samplefiles)
    assert len(iraf.utils.file_handler('@' + atlist)) == tot


# explicitly call ccdproc with every parameter set to what we want
defaultargs = {'output': None, 'ccdtype': 'object', 'noproc': False,
               'fixpix': False, 'overscan': False, 'trim': False,
               'zerocor': False,
               'darkcor': False, 'flatcor': False, 'illumcor': False,
               'fringecor': False, 'readcor': False, 'scancor': False,
               'readaxis': 'line', 'fixfile': None, 'biassec': None,
               'trimsec': None, 'zero': None, 'dark': None, 'flat': None,
               'illum': None, 'fringe': None, 'minreplace': 1.,
               'scantype': 'shortscan', 'nscan': 1, 'interactive': False,
               'overscan_function': 'legendre', 'order': 1, 'sample': '*',
               'naverage': 1, 'niterate': 1, 'low_reject': 3.,
               'high_reject': 3., 'grow': 0., 'instrument': None,
               'pixeltype': "real", 'logfile': None, 'verbose': False}


def test_zerocor(combine_dir):
    basedir = str(combine_dir)

    # create some simple files for testing
    nx = 50
    ny = 90
    nimg = 3
    baseval = 100

    # create input files
    inputs = []
    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float) * baseval
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)

    # make the input list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')
    inlist = '@' + inlist

    # make the output list
    outlist = os.path.join(basedir, 'outfiles.txt')
    outlists = []
    with open(outlist, 'w') as ff:
        for ifile in inputs:
            of = ifile[:-5] + '.out.fits'
            ff.write(of + '\n')
            outlists.append(of)

    outlist = '@' + outlist

    # make a zero image
    zeroval = 10
    zerofile = os.path.join(basedir, 'testzero.fits')

    arr = np.ones((nx, ny), dtype=float) * zeroval

    hdu = fits.PrimaryHDU(arr)
    hdu.header['imagetyp'] = 'zero'
    hdu.writeto(zerofile, overwrite=True)

    # do my processing
    myargs = copy.deepcopy(defaultargs)
    myargs['output'] = outlist

    # which parts we want to do
    myargs['zerocor'] = True
    myargs['zero'] = zerofile

    iraf.ccdproc(inlist, **myargs)

    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data, baseval - zeroval)
        assert len(hdr[0].header['zerocor']) > 0
        assert len(hdr[0].header['ccdproc']) > 0
        assert hdr[0].header['ccdsec'] == f'[1:{ny:d},1:{nx:d}]'
        hdr.close()

    # test no outputs and replacing the input
    myargs['output'] = None
    iraf.ccdproc(inlist, **myargs)

    oldtimes = []
    for ifile in inputs:
        hdr = fits.open(ifile)
        assert hdr[0].header['ccdsec'] == f'[1:{ny:d},1:{nx:d}]'
        assert len(hdr[0].header['zerocor']) > 0
        oldtimes.append(hdr[0].header['zerocor'])
        hdr.close()

    # wait a second and do it again
    time.sleep(2)
    iraf.ccdproc(inlist, **myargs)
    # make sure nothing has happened since inputs are already processed
    for ii, ifile in enumerate(outlists):
        hdr = fits.open(ifile)
        assert hdr[0].header['zerocor'] == oldtimes[ii]
        hdr.close()

    # test datasec of input image as just a subset of the image

    # add the datasec parameter
    dsec = '[5:30,10:40]'
    inputs = []
    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float) * baseval
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['datasec'] = dsec
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)

    myargs['output'] = outlist
    iraf.ccdproc(inlist, **myargs)

    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data[9:40, 4:30], baseval - zeroval)
        data = hdr[0].data * 1
        data[9:40, 4:30] += zeroval
        assert np.allclose(data, baseval)
        assert hdr[0].header['datasec'] == dsec
        hdr.close()

    # test datasec/ccdsec matches of input and zero image

    # add the datasec/ccdsec parameter
    dsec = '[5:30,10:40]'
    ccdsec = '[50:75,100:130]'
    baseoff = np.arange(ny) + 46  # ccdsec1 - dsec1 + 1
    inputs = []
    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float) * baseval
        arr += baseoff
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['datasec'] = dsec
        hdu.header['ccdsec'] = ccdsec
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)

    # make a zero image
    zeroval = 10
    zerofile = os.path.join(basedir, 'testzero.fits')
    # nx, ny = 50, 90
    zdsec = '[30:90, 5:45]'
    zccdsec = '[20:80, 95:135]'

    zarr = np.ones((nx, ny), dtype=float) * zeroval
    zbaseoff = np.arange(ny) - 9  # ccdsec1 - dsec1 + 1
    zarr += zbaseoff

    hdu = fits.PrimaryHDU(zarr)
    hdu.header['imagetyp'] = 'zero'
    hdu.header['datasec'] = zdsec
    hdu.header['ccdsec'] = zccdsec
    hdu.writeto(zerofile, overwrite=True)

    iraf.ccdproc(inlist, **myargs)

    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data[9:40, 4:30], baseval - zeroval)
        data = hdr[0].data * 1
        data[9:40, 4:30] += zeroval + baseoff[4:30]
        assert np.allclose(data, arr)
        assert hdr[0].header['datasec'] == dsec
        assert hdr[0].header['ccdsec'] == ccdsec
        hdr.close()
