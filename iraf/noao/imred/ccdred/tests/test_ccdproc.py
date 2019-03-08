import iraf
import pytest
import os
import numpy as np
from astropy.io import fits
import copy
import time


def test_ccd_section():
    defaults = [5, 6, 7, 8, 9, 10]
    d1, d2, d3, d4, d5, d6 = defaults
    # get the defaults
    x0, x1, xs, y0, y1, ys = iraf.ccdred.ccd_section('[:,:]', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    x0, x1, xs, y0, y1, ys = iraf.ccdred.ccd_section(None, defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    x0, x1, xs, y0, y1, ys = iraf.ccdred.ccd_section('  ', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    # all the options
    x0, x1, xs, y0, y1, ys = iraf.ccdred.ccd_section('[1:5:3,5:20:-4]',
                                                     defaults=defaults)
    assert (x0 == 1 and x1 == 5 and xs == 3 and y0 == 5 and y1 == 20 and
            ys == -4)

    # no step size
    x0, x1, xs, y0, y1, ys = iraf.ccdred.ccd_section('[1:5,5:20]',
                                                     defaults=defaults)
    assert (x0 == 1 and x1 == 5 and xs == d3 and y0 == 5 and y1 == 20 and
            ys == d6)

    # single row
    x0, x1, xs, y0, y1, ys = iraf.ccdred.ccd_section('[9,5:20]',
                                                     defaults=defaults)
    assert (x0 == 9 and x1 == 9 and xs == d3 and y0 == 5 and y1 == 20 and
            ys == d6)

    # test bad braces
    with pytest.raises(Exception):
        iraf.ccdred.ccd_section('[:,:')
    with pytest.raises(Exception):
        iraf.ccdred.ccd_section(':,:]')
    with pytest.raises(Exception):
        iraf.ccdred.ccd_section('1:2,5:34')

    # test wrong dimensions
    with pytest.raises(Exception):
        iraf.ccdred.ccd_section('[1:5]')
    with pytest.raises(Exception):
        iraf.ccdred.ccd_section('[1:5, 1:5, 1:5]')

    # test bad inputs in a dimension
    with pytest.raises(Exception):
        iraf.ccdred.ccd_section('[1:3:4:5, :]')

    # test bad size of defaults
    with pytest.raises(Exception):
        iraf.ccdred.ccd_section('[:,:]', defaults=[1, 2, 3])


@pytest.mark.parametrize("listtype", iraf.ccdred.imagetypes)
def test_cal_list(tmpdir, listtype):
    basedir = str(tmpdir)

    inst = iraf.ccdred.Instrument()
    # set up the necessary empty lists
    calimages, nscans, caltypes, subsets = [], [], [], []
    scantype = 'shortscan'
    nscan = 1
    scancor = True

    # create some simple files for testing
    nx = 50
    ny = 90
    baseval = 100

    # test files not found
    foo = os.path.join(basedir, 'foo.fits')
    with pytest.raises(Exception):
        iraf.ccdred.cal_list([foo], listtype, inst, calimages, nscans, caltypes,
                             subsets, scantype, nscan, scancor)

    lt = len(iraf.ccdred.imagetypes)
    # test both 'imagetyp' and custom instrument value
    for instval in ['imagetyp', 'ityp']:
        inst.parameters['imagetyp'] = instval

        # create input files
        inputs = []
        # set up the necessary empty lists
        calimages, nscans, caltypes, subsets = [], [], [], []

        # test all possible image types
        for jj in np.arange(len(iraf.ccdred.imagetypes)):
            arr = np.ones((nx, ny), dtype=float) * baseval
            hdu = fits.PrimaryHDU(arr)
            hdu.header[instval] = iraf.ccdred.imagetypes[jj]
            inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
            hdu.writeto(inim, overwrite=True)
            inputs.append(inim)

        iraf.ccdred.cal_list(inputs, listtype, inst, calimages, nscans,
                             caltypes, subsets, scantype, nscan, scancor)

        # XXX: still to do: test the output of the other lists
        if listtype == 'unknown':
            assert (len(calimages) == 5 and len(nscans) == 5 and
                    len(caltypes) == 5 and len(subsets) == 5)
            # make sure there's one of each allowable type
            assert len({'zero', 'dark', 'flat', 'illum', 'fringe'} &
                       set(caltypes)) == 5

        elif listtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
            assert (len(calimages) == lt and len(nscans) == lt and
                    len(caltypes) == lt and len(subsets) == lt)
            # make sure all types are the same as the input type
            for itype in caltypes:
                assert itype == listtype

        else:
            assert (len(calimages) == 0 and len(nscans) == 0 and
                    len(caltypes) == 0 and len(subsets) == 0)

        # XXX: check for duplication


def test_ccdnscan():
    pass


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


def test_basics(tmpdir):
    basedir = str(tmpdir)
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
    outlists = []
    outlist = os.path.join(basedir, 'outfiles.txt')
    with open(outlist, 'w') as ff:
        for ifile in inputs:
            of = ifile[:-5] + '.out.fits'
            outlists.append(of)
            ff.write(of + '\n')
    outlist = '@' + outlist

    # make a zero image
    zeroval = 10
    zerofile = os.path.join(basedir, 'testzero.fits')

    arr = np.ones((nx, ny), dtype=float) * zeroval
    hdu = fits.PrimaryHDU(arr)
    hdu.header['imagetyp'] = 'zero'
    hdu.writeto(zerofile, overwrite=True)

    myargs = copy.deepcopy(defaultargs)
    # which parts we want to do
    myargs['zerocor'] = True
    myargs['zero'] = zerofile

    # test no outputs works
    iraf.ccdred.ccdproc(inputs, **myargs)
    for ifile in inputs:
        with fits.open(ifile) as hdr:
            assert len(hdr[0].header['zerocor']) > 0

    # reset the input files
    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float) * baseval
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)

    # test output gets created
    myargs['output'] = outlists
    iraf.ccdred.ccdproc(inputs, **myargs)
    for ifile in outlists:
        with fits.open(ifile) as hdr:
            assert len(hdr[0].header['zerocor']) > 0
        # remove the file to test the same one again later
        os.remove(ifile)

    # garbage scantype value
    myargs['scantype'] = 'foo'
    with pytest.raises(Exception):
        iraf.ccdred.ccdproc(inputs, **myargs)
    myargs['scantype'] = 'shortscan'

    # test bad output sizes
    myargs['output'] = outlists*2
    with pytest.raises(Exception):
        iraf.ccdred.ccdproc(inputs, **myargs)
    myargs['output'] = outlists[0]
    with pytest.raises(Exception):
        iraf.ccdred.ccdproc(inputs, **myargs)
    # make sure the output files don't exist
    for ifile in outlists:
        assert not os.path.exists(ifile)

    # test input and output lists
    myargs['output'] = outlist
    iraf.ccdred.ccdproc(inlist, **myargs)
    for ifile in outlists:
        with fits.open(ifile) as hdr:
            assert len(hdr[0].header['zerocor']) > 0

    # test if one of the input files doesn't exist
    os.remove(inputs[1])
    with pytest.raises(Exception):
        iraf.ccdred.ccdproc(inlist, **myargs)


def test_zerocor(tmpdir):
    basedir = str(tmpdir)

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
    outlists = []
    for ifile in inputs:
        of = ifile[:-5] + '.out.fits'
        outlists.append(of)

    # make a zero image
    zeroval = 10
    zerofile = os.path.join(basedir, 'testzero.fits')

    arr = np.ones((nx, ny), dtype=float) * zeroval

    hdu = fits.PrimaryHDU(arr)
    hdu.header['imagetyp'] = 'zero'
    hdu.writeto(zerofile, overwrite=True)

    # do my processing
    myargs = copy.deepcopy(defaultargs)
    myargs['output'] = outlists

    # which parts we want to do
    myargs['zerocor'] = True
    myargs['zero'] = zerofile

    iraf.ccdred.ccdproc(inlist, **myargs)

    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data, baseval - zeroval)
        assert len(hdr[0].header['zerocor']) > 0
        assert len(hdr[0].header['ccdproc']) > 0
        assert hdr[0].header['ccdsec'] == f'[1:{ny:d},1:{nx:d}]'
        hdr.close()

    # test no outputs and replacing the input
    myargs['output'] = None
    iraf.ccdred.ccdproc(inlist, **myargs)

    oldtimes = []
    for ifile in inputs:
        hdr = fits.open(ifile)
        assert hdr[0].header['ccdsec'] == f'[1:{ny:d},1:{nx:d}]'
        assert len(hdr[0].header['zerocor']) > 0
        oldtimes.append(hdr[0].header['zerocor'])
        hdr.close()

    # wait a second and do it again
    time.sleep(2)
    iraf.ccdred.ccdproc(inlist, **myargs)
    # make sure nothing has happened since inputs are already processed
    for ii, ifile in enumerate(inputs):
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

    myargs['output'] = outlists
    iraf.ccdred.ccdproc(inlist, **myargs)

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

    iraf.ccdred.ccdproc(inlist, **myargs)

    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data[9:40, 4:30], baseval - zeroval)
        data = hdr[0].data * 1
        data[9:40, 4:30] += zeroval + baseoff[4:30]
        assert np.allclose(data, arr)
        assert hdr[0].header['datasec'] == dsec
        assert hdr[0].header['ccdsec'] == ccdsec
        hdr.close()


# things to test:

# every ccdtype for cal_open. Do this by giving all 5 types of files
# and trying to process all 5 types at once for each ccdtype.

# test having the calibration images in the input list.
