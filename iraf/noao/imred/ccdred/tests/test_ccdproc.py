import copy
import datetime
import os

import numpy as np
import pytest
from astropy.io import fits

import iraf
from .. import ccdproc_routines as ccdr
from ..utils import CCDProcError, CCDProcWarning

# explicitly call ccdproc with every parameter set to what we want
defaultargs = {'output': None, 'ccdtype': 'object', 'noproc': False,
               'fixpix': False, 'overscan': False, 'trim': False,
               'zerocor': False,
               'darkcor': False, 'flatcor': False, 'illumcor': False,
               'fringecor': False, 'readcor': False, 'scancor': False,
               'readaxis': 'line', 'fixfile': None, 'biassec': 'image',
               'trimsec': 'image', 'zero': None, 'dark': None, 'flat': None,
               'illum': None, 'fringe': None, 'minreplace': 1.,
               'scantype': 'shortscan', 'nscan': 1, 'interactive': False,
               'overscan_function': 'legendre', 'order': 1, 'sample': '*',
               'naverage': 1, 'niterate': 1, 'low_reject': 3.,
               'high_reject': 3., 'grow': 0., 'instrument': None,
               'pixeltype': 'real', 'logfile': None, 'verbose': False,
               'overwrite': True}


def test_ccd_section():
    defaults = [5, 6, 7, 8, 9, 10]
    d1, d2, d3, d4, d5, d6 = defaults
    # get the defaults
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[:,:]', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    x0, x1, xs, y0, y1, ys = ccdr.ccd_section(None, defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('  ', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == d6)

    # all the options
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[1:5:3,5:-20:-4]',
                                              defaults=defaults)
    assert (x0 == 1 and x1 == 5 and xs == 3 and y0 == 5 and y1 == -20 and
            ys == -4)

    # no step size
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[1:5,5:20]', defaults=defaults)
    assert (x0 == 1 and x1 == 5 and xs == d3 and y0 == 5 and y1 == 20 and
            ys == d6)

    # single row
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[9,5:20]', defaults=defaults)
    assert (x0 == 9 and x1 == 9 and xs == d3 and y0 == 5 and y1 == 20 and
            ys == d6)

    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[,5:20]', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == 5 and y1 == 20 and
            ys == d6)

    # test missing values
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[:, 1:]', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == 1 and y1 == d5 and
            ys == d6)
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[:, 1::3]', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == 1 and y1 == d5 and
            ys == 3)
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[:, ::3]', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == d5 and
            ys == 3)
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[:, :2:3]', defaults=defaults)
    assert (x0 == d1 and x1 == d2 and xs == d3 and y0 == d4 and y1 == 2 and
            ys == 3)

    # test one dimension
    x0, x1, xs, y0, y1, ys = ccdr.ccd_section('[1:5:2]', defaults=defaults)
    assert (x0 == 1 and x1 == 5 and xs == 2 and y0 == d4 and y1 == d5 and
            ys == d6)

    # test bad braces
    with pytest.raises(ValueError):
        ccdr.ccd_section('[:,:')
    with pytest.raises(ValueError):
        ccdr.ccd_section(':,:]')
    with pytest.raises(ValueError):
        ccdr.ccd_section('1:2,5:34')

    # test non integer input
    with pytest.raises(ValueError):
        ccdr.ccd_section('[a,:]')
    with pytest.raises(ValueError):
        ccdr.ccd_section('[a:b,:]')
    with pytest.raises(ValueError):
        ccdr.ccd_section('[1:3:a,:]')

    with pytest.raises(ValueError):
        ccdr.ccd_section('[2.2,:]')
    with pytest.raises(ValueError):
        ccdr.ccd_section('[2.2:3.3,:]')
    with pytest.raises(ValueError):
        ccdr.ccd_section('[1:3:2.,:]')

    # test wrong dimensions
    with pytest.raises(ValueError):
        ccdr.ccd_section('[1:5, 1:5, 1:5]')

    # test bad inputs in a dimension
    with pytest.raises(ValueError):
        ccdr.ccd_section('[1:3:4:5, :]')

    # test bad size of defaults
    with pytest.raises(IndexError):
        ccdr.ccd_section('[:,:]', defaults=[1, 2, 3])


@pytest.mark.parametrize("listtype", iraf.ccdred._imagetypes)
def test_ccdnscan(tmpdir, listtype):
    basedir = str(tmpdir)

    inst = iraf.ccdred.Instrument()

    # create simple file for testing
    nx = 50
    ny = 90
    baseval = 100
    arr = np.ones((nx, ny), dtype=float) * baseval
    hdu = fits.PrimaryHDU(arr)
    hdu.header['imagetyp'] = listtype
    hdu.header['nscanrow'] = 5
    hdu.header['srowcust'] = 10
    inim = os.path.join(basedir, f'testimg.fits')
    hdu.writeto(inim, overwrite=True)

    # test getting things from the header first
    hdu = iraf.sys.image_open(inim)
    assert ccdr.ccdnscan(hdu, inst, listtype, 'shortscan', 3, True) == 5

    inst.parameters['nscanrow'] = 'srowcust'
    assert ccdr.ccdnscan(hdu, inst, listtype, 'shortscan', 3, True) == 10

    inst.parameters['nscanrow'] = 'foo'
    inst.defaults['nscanrow'] = 7
    assert ccdr.ccdnscan(hdu, inst, listtype, 'shortscan', 3, True) == 7

    # check things not in the image header and no defaults
    inst.defaults['nscanrow'] = None

    # check failures
    if listtype not in ['zero', 'dark', 'flat', 'illum', 'fringe']:
        with pytest.raises(ValueError):
            ccdr.ccdnscan(hdu, inst, listtype, 'noscan', 3, True)

    # test the 4 combos of scantype and scancor
    sscan = 'shortscan'
    scor = True
    res = ccdr.ccdnscan(hdu, inst, listtype, sscan, 3, scor)

    if listtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
        assert res == 1
    else:
        assert res == 3

    sscan = 'shortscan'
    scor = False
    res = ccdr.ccdnscan(hdu, inst, listtype, sscan, 3, scor)

    if listtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
        assert res == 1
    else:
        assert res == 3

    sscan = 'longscan'
    scor = False
    res = ccdr.ccdnscan(hdu, inst, listtype, sscan, 3, scor)

    assert res == 1

    sscan = 'longscan'
    scor = True
    res = ccdr.ccdnscan(hdu, inst, listtype, sscan, 3, scor)

    if listtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
        assert res == 1
    else:
        assert res is None

    hdu.close()


@pytest.mark.parametrize("listtype", iraf.ccdred._imagetypes)
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

    if listtype != 'unknown':
        with pytest.raises(OSError):
            ccdr.cal_list([foo], listtype, inst, calimages, nscans, caltypes,
                          subsets, scantype, nscan, scancor)
    else:
        ccdr.cal_list([foo], listtype, inst, calimages, nscans, caltypes,
                      subsets, scantype, nscan, scancor)
        assert len(calimages) == 0

    lt = len(iraf.ccdred._imagetypes)
    # test both 'imagetyp' and custom instrument value
    for instval in ['imagetyp', 'ityp']:
        inst.parameters['imagetyp'] = instval

        # create input files
        inputs = []
        # set up the necessary empty lists
        calimages, nscans, caltypes, subsets = [], [], [], []

        # test all possible image types
        for jj in np.arange(len(iraf.ccdred._imagetypes)):
            arr = np.ones((nx, ny), dtype=float) * baseval
            hdu = fits.PrimaryHDU(arr)
            hdu.header[instval] = iraf.ccdred._imagetypes[jj]
            hdu.header['subset'] = iraf.ccdred._imagetypes[jj]
            hdu.header['nscanrow'] = jj
            inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
            hdu.writeto(inim, overwrite=True)
            # add each image twice and make sure we only find them once in the
            # list as a check on catching duplicates.
            inputs.append(inim)
            inputs.append(inim)

        ccdr.cal_list(inputs, listtype, inst, calimages, nscans,
                      caltypes, subsets, scantype, nscan, scancor)

        if listtype == 'unknown':
            assert (len(calimages) == 5 and len(nscans) == 5 and
                    len(caltypes) == 5 and len(subsets) == 5)
            # make sure there's one of each allowable type
            assert len({'zero', 'dark', 'flat', 'illum', 'fringe'} &
                       set(caltypes)) == 5
            # make sure the subset and nscan is right
            for ii, ival in enumerate(subsets):
                ind = iraf.ccdred._imagetypes.index(ival)
                assert ival == iraf.ccdred._imagetypes[ind]
                assert nscans[ii] == ind

        elif listtype in ['zero', 'dark', 'flat', 'illum', 'fringe']:
            assert (len(calimages) == lt and len(nscans) == lt and
                    len(caltypes) == lt and len(subsets) == lt)
            # make sure all types are the same as the input type
            for itype in caltypes:
                assert itype == listtype
            # make sure the subset and nscan is right
            for ii, ival in enumerate(subsets):
                ind = iraf.ccdred._imagetypes.index(ival)
                assert ival == iraf.ccdred._imagetypes[ind]
                assert nscans[ii] == ind

        else:
            assert (len(calimages) == 0 and len(nscans) == 0 and
                    len(caltypes) == 0 and len(subsets) == 0)


def test_cal_scan():
    inpath = '/foo/bar/file.txt'
    root, ext = os.path.splitext(inpath)
    assert ccdr.cal_scan(inpath, 1, False) == inpath
    assert ccdr.cal_scan(inpath, 6, False) == inpath
    assert ccdr.cal_scan(inpath, 1, True) == inpath
    assert ccdr.cal_scan(inpath, 6, True) == f"{root}.{6:d}{ext}"
    assert ccdr.cal_scan(inpath, None, True) == f"{root}.{1:d}d{ext}"
    with pytest.raises(ValueError):
        ccdr.cal_scan(inpath, 5.5, True)
    with pytest.raises(ValueError):
        ccdr.cal_scan(inpath, 'str', True)


@pytest.mark.parametrize("imtype", iraf.ccdred._imagetypes)
def test_cal_image(tmpdir, imtype):
    basedir = str(tmpdir)

    inst = iraf.ccdred.Instrument()

    # create simple file for testing
    nx = 50
    ny = 90
    baseval = 100
    arr = np.ones((nx, ny), dtype=float) * baseval
    hdu = fits.PrimaryHDU(arr)
    hdu.header['imagetyp'] = imtype
    hdu.header['nscanrow'] = 5
    hdu.header['srowcust'] = 10
    inim = os.path.join(basedir, f'testimg.fits')
    hdu.writeto(inim, overwrite=True)

    # open the file for feeding into all the tests
    inopen = iraf.sys.image_open(inim)

    # basic test of the function
    calim = os.path.join(basedir, f'testcal.fits')
    cals = ([calim], [1], [imtype], [''])

    # non calibration images fail no matter what
    if imtype not in ['zero', 'dark', 'flat', 'illum', 'fringe']:
        with pytest.raises(ValueError):
            ccdr.cal_image(inopen, inst, imtype, 1, cals, True)
        return
    else:
        assert ccdr.cal_image(inopen, inst, imtype, 1, cals, True) == cals[0][0]
    # don't deal with non-calibration types anymore
    assert imtype in ['zero', 'dark', 'flat', 'illum', 'fringe']

    # for all image types, check error with 0 images
    calimages, nscans, caltypes, subsets = [], [], [], []
    for itype in iraf.ccdred._imagetypes:
        calim = os.path.join(basedir, f'testcal.{itype}.fits')
        if itype != imtype:
            calimages.append(calim)
            nscans.append(1)
            caltypes.append(itype)
            subsets.append('')
    cals = (calimages, nscans, caltypes, subsets)

    with pytest.raises(CCDProcError):
        ccdr.cal_image(inopen, inst, imtype, 1, cals, True)
    # add the image in and make sure we get the expected result
    calim = os.path.join(basedir, f'testcal.{imtype}.fits')
    calimages.append(calim)
    nscans.append(1)
    caltypes.append(imtype)
    subsets.append('')
    cals = (calimages, nscans, caltypes, subsets)
    assert ccdr.cal_image(inopen, inst, imtype, 1, cals, True) == calim

    # check selected image doesn't equal input image
    calimages, nscans, caltypes, subsets = [inim], [1], [imtype], ['']
    cals = (calimages, nscans, caltypes, subsets)
    with pytest.raises(CCDProcError):
        ccdr.cal_image(inopen, inst, imtype, 1, cals, True)

    # fail to reach cal_scan
    calim = os.path.join(basedir, f'testcal.{imtype}.fits')

    cals = ([calim], [6], [imtype], [''])
    with pytest.raises(CCDProcError):
        ccdr.cal_image(inopen, inst, imtype, 1, cals, True)

    cals = ([calim], [5], [imtype], [''])
    with pytest.raises(CCDProcError):
        ccdr.cal_image(inopen, inst, imtype, 6, cals, True)

    # make sure cal_scan is working
    cals = ([calim], [1], [imtype], [''])

    # no scancor, so no cal_scan
    assert ccdr.cal_image(inopen, inst, imtype, 6, cals, False) == calim

    root, ext = os.path.splitext(calim)
    # nscan of None
    retval = f"{root}.1d{ext}"
    assert ccdr.cal_image(inopen, inst, imtype, None, cals, True) == retval

    # nscan not 1
    retval = f"{root}.{6:d}{ext}"
    assert ccdr.cal_image(inopen, inst, imtype, 6, cals, True) == retval

    # for all, check error with 2+ images of the same image type and nscan
    cals = ([calim + '1', calim + '2'], [5, 5], [imtype, imtype], ['', ''])
    with pytest.raises(CCDProcError):
        ccdr.cal_image(inopen, inst, imtype, 1, cals, True)

    # pick the right nscan
    cals = ([calim + '1', calim + '2'], [3, 2], [imtype, imtype], ['', ''])
    assert ccdr.cal_image(inopen, inst, imtype, 2, cals, True) == calim + '2'

    cals = ([calim + '1', calim + '2'], [3, 1], [imtype, imtype], ['', ''])
    assert ccdr.cal_image(inopen, inst, imtype, 2, cals, False) == calim + '2'

    inopen.close()

    # add a subset to our input image
    hdu.header['subset'] = 'red'
    inim = os.path.join(basedir, f'testimg.fits')
    hdu.writeto(inim, overwrite=True)

    # open the file for feeding into all the tests
    inopen = iraf.sys.image_open(inim)

    # only check for correct subset for flat, fringe, illum
    cals = ([calim + '1', calim + '2'], [1, 1], [imtype, imtype],
            ['blue', 'red'])
    if imtype in ['flat', 'illum', 'fringe']:
        cor = calim + '2'
        assert ccdr.cal_image(inopen, inst, imtype, 1, cals, False) == cor
    else:
        with pytest.raises(CCDProcError):
            ccdr.cal_image(inopen, inst, imtype, 1, cals, False)

    inopen.close()


def test_logstring(tmpdir, capsys):
    basedir = str(tmpdir)

    # create simple file for testing
    nx = 50
    ny = 90
    baseval = 100
    arr = np.ones((nx, ny), dtype=float) * baseval
    hdu = fits.PrimaryHDU(arr)
    inim = os.path.join(basedir, f'testimg.fits')
    hdu.writeto(inim, overwrite=True)

    # open the file for feeding into all the tests
    inopen = iraf.sys.image_open(inim)

    instr = 'foo'
    retstr = ccdr.logstring(instr, inopen, False, None)

    # no outputs
    out, err = capsys.readouterr()
    assert out == '' and err == ''

    # make sure the dates are understandable and the same
    now = datetime.datetime.now()
    dfmt = '%Y-%m-%d %H:%M:%S'
    nowstr = now.strftime(dfmt)
    filetime = datetime.datetime.strptime(retstr[:len(nowstr)], dfmt)
    assert np.abs((now - filetime).total_seconds()) <= 2

    # extra character is the space between them
    assert len(retstr) == (len(nowstr) + len(instr) + 1)

    # we got our input string back
    assert retstr[len(nowstr) + 1:] == instr

    # check log files and verbose output
    log = os.path.join(basedir, f'log.txt')
    retstr = ccdr.logstring(instr, inopen, True, log)
    now = datetime.datetime.now()

    # make sure we haven't changed the return string while printing stuff
    assert len(retstr) == (len(nowstr) + len(instr) + 1)

    with open(log, 'r') as ff:
        lines = ff.readlines()
    assert len(lines) == 1
    line = lines[0]

    # make sure output and logfile have the same contents
    out, err = capsys.readouterr()
    assert len(err) == 0
    assert len(out) == len(line)
    assert out == line

    assert out[:len(inim)] == inim
    # account for the newline at the end, but otherwise keep the input string
    assert out[-len(instr) - 1:-1] == instr

    cutstr = out[len(inim) + 2:len(inim) + 2 + len(nowstr)]
    filetime = datetime.datetime.strptime(cutstr, dfmt)
    assert np.abs((now - filetime).total_seconds()) <= 2

    inopen.close()


def test_already_processed(tmpdir):
    basedir = str(tmpdir)

    # create simple file for testing
    nx = 50
    ny = 90
    baseval = 100
    arr = np.ones((nx, ny), dtype=float) * baseval
    hdu = fits.PrimaryHDU(arr)
    hdu.header['trim'] = 'foo'
    hdu.header['custom'] = 'bar'
    inim = os.path.join(basedir, f'testimg.fits')
    hdu.writeto(inim, overwrite=True)

    inst = iraf.ccdred.Instrument()

    # open the file for feeding into all the tests
    inopen = iraf.sys.image_open(inim)

    # default of None is not 'foo'
    assert ccdr.already_processed(inopen, inst, 'trim')

    # default matches what's found in file
    inst.defaults['trim'] = 'foo'
    assert not ccdr.already_processed(inopen, inst, 'trim')

    inst.parameters['trim'] = 'custom'
    # default of 'foo' is not 'bar'
    assert ccdr.already_processed(inopen, inst, 'trim')

    # default matches what's found in file
    inst.defaults['trim'] = 'bar'
    assert not ccdr.already_processed(inopen, inst, 'trim')

    # if it's not in the header
    assert not ccdr.already_processed(inopen, inst, 'zerocor')

    inopen.close()


@pytest.mark.parametrize("imtype", iraf.ccdred._imagetypes)
def test_ccdcheck(tmpdir, imtype):
    basedir = str(tmpdir)

    inst = iraf.ccdred.Instrument()

    flags = copy.deepcopy(defaultargs)

    # create simple file for testing
    nx = 50
    ny = 90
    baseval = 100
    arr = np.ones((nx, ny), dtype=float) * baseval
    basehdu = fits.PrimaryHDU(arr)

    inim = os.path.join(basedir, f'testimg.fits')
    basehdu.writeto(inim, overwrite=True)

    inopen = iraf.sys.image_open(inim)

    # everything should be false
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check trim
    flags['trim'] = True
    assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['trim'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check fixpix
    flags['fixpix'] = True
    assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['fixpix'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check overscan
    flags['overscan'] = True
    assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['overscan'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check readcor
    flags['readcor'] = True
    if imtype in ['zero']:
        assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    else:
        assert not ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['readcor'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check zerocor
    flags['zerocor'] = True
    if imtype not in ['zero']:
        assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    else:
        assert not ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['zerocor'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check darkcor
    flags['darkcor'] = True
    if imtype not in ['zero', 'dark']:
        assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    else:
        assert not ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['darkcor'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check scancor
    flags['scancor'] = True
    if imtype in ['flat']:
        assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    else:
        assert not ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['scancor'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check flatcor
    flags['flatcor'] = True
    if imtype not in ['zero', 'dark', 'flat']:
        assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    else:
        assert not ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['flatcor'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check illumcor
    flags['illumcor'] = True
    if imtype not in ['zero', 'dark', 'flat', 'illum']:
        assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    else:
        assert not ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['illumcor'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    # check fringecor
    flags['fringecor'] = True
    if imtype not in ['zero', 'dark', 'flat', 'illum']:
        assert ccdr.ccdcheck(inopen, inst, imtype, flags)
    else:
        assert not ccdr.ccdcheck(inopen, inst, imtype, flags)
    inopen[0].header['fringcor'] = 'foo'
    assert not ccdr.ccdcheck(inopen, inst, imtype, flags)

    inopen.close()


def test_ccdproc_basics(tmpdir):
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

    # test no outputs works even with overwrite off
    myargs['overwrite'] = False
    iraf.ccdred.ccdproc(inputs, **myargs)
    for ifile in inputs:
        with fits.open(ifile) as hdr:
            assert len(hdr[0].header['zerocor']) > 0
    myargs['overwrite'] = True

    # reset the input files
    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float) * baseval
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)

    # test output files
    myargs['output'] = outlist
    myargs['zerocor'] = False

    # test no output files if nothing is requested to be done
    iraf.ccdred.ccdproc(inputs, **myargs)
    for ifile in outlists:
        assert not os.path.exists(ifile)

    myargs['zerocor'] = True

    iraf.ccdred.ccdproc(inputs, **myargs)
    mtimes = []
    for ifile in outlists:
        with fits.open(ifile) as hdr:
            assert len(hdr[0].header['zerocor']) > 0
        mtimes.append(os.path.getmtime(ifile))

    # test overwrite flag works
    myargs['overwrite'] = False
    with pytest.raises(OSError):
        iraf.ccdred.ccdproc(inputs, **myargs)
    for ii in np.arange(len(outlists)):
        assert os.path.getmtime(outlists[ii]) == mtimes[ii]

    myargs['overwrite'] = True
    iraf.ccdred.ccdproc(inputs, **myargs)
    for ii in np.arange(len(outlists)):
        assert os.path.getmtime(outlists[ii]) != mtimes[ii]

    # remove the files to test the same one again later
    for ifile in outlists:
        os.remove(ifile)

    # garbage scantype value
    myargs['scantype'] = 'foo'
    with pytest.raises(ValueError):
        iraf.ccdred.ccdproc(inputs, **myargs)
    myargs['scantype'] = 'shortscan'

    # garbage readaxis value
    myargs['readaxis'] = 'foo'
    with pytest.raises(ValueError):
        iraf.ccdred.ccdproc(inputs, **myargs)
    myargs['readaxis'] = 'line'

    # no files got created
    for ifile in outlists:
        assert not os.path.exists(ifile)

    # test bad output sizes
    myargs['output'] = outlists * 2
    with pytest.raises(ValueError):
        iraf.ccdred.ccdproc(inputs, **myargs)
    myargs['output'] = outlists[0]
    with pytest.raises(ValueError):
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
        # remove the file to test the same one again later
        os.remove(ifile)

    # test if one of the input files doesn't exist
    os.remove(inputs[1])
    with pytest.raises(ValueError):
        iraf.ccdred.ccdproc(inlist, **myargs)

    # what if that input file is corrupted?
    with open(inputs[1], 'w') as ff:
        ff.write('foo')

    with pytest.warns(CCDProcWarning):
        iraf.ccdred.ccdproc(inlist, **myargs)

    for ii, ifile in enumerate(outlists):
        if ii != 1:
            with fits.open(ifile) as hdr:
                assert len(hdr[0].header['zerocor']) > 0
                # remove the file to test the same one again later
                os.remove(ifile)
        else:
            assert not os.path.exists(ifile)

    # reset the input files
    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float) * baseval
        hdu = fits.PrimaryHDU(arr)
        if jj == 1:
            hdu.header['imagetyp'] = 'none'
        else:
            hdu.header['imagetyp'] = 'object'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)

    # what if the image types don't match?
    with pytest.warns(CCDProcWarning):
        iraf.ccdred.ccdproc(inlist, **myargs)

    for ii, ifile in enumerate(outlists):
        if ii != 1:
            with fits.open(ifile) as hdr:
                assert len(hdr[0].header['zerocor']) > 0
                # remove the file to test the same one again later
                os.remove(ifile)
        else:
            assert not os.path.exists(ifile)


def test_ccdproc_cal_open():
    # XXX: finish
    # go through all ccdtypes.

    # make sure if an object is passed in as a list of a specific type it is
    # processed that way regardless of what's in header

    # make sure that passing in all files as one input list gets things sorted
    # into the right lists as well

    # check subsets are handled appropriately

    # try feeding the same image as multiple types of cals

    # every ccdtype for cal_open. Do this by giving all 5 types of files
    # and trying to process all 5 types at once for each ccdtype.
    pass


def test_ccdproc_noproc():
    # XXX: finish
    # make sure everything is printed

    # make sure no output file is created and the input file isn't changed

    pass


pixeltypes = ['short', 'ushort', 'integer', 'real', 'double', None, 'foo']
realpixtypes = ['short', 'ushort', 'integer', 'real', 'double']


@pytest.mark.parametrize("intype", realpixtypes)
@pytest.mark.parametrize("outtype", pixeltypes)
def test_ccdproc_pixeltype(tmpdir, intype, outtype):
    # make sure rounding works right
    # what does IRAF really do if input is int, zero is float?
    basedir = str(tmpdir)

    # create some simple files for testing
    nx = 50
    ny = 90
    nimg = 3
    baseval = 100

    dd = {'short': np.short, 'ushort': np.ushort, 'integer': np.int32,
          'real': np.single, 'double': np.double}
    ind = dd[intype]

    # create input files
    inputs = []
    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=ind) * baseval
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

    # make a zero image with a float so that the output should be
    # 89.7, so we can test int rounding
    zeroval = 10.3
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
    myargs['pixeltype'] = outtype

    # test all the warnings and errors of trying to convert to requested
    # output pixel types
    if outtype is None or outtype == 'foo':
        with pytest.raises(ValueError):
            iraf.ccdred.ccdproc(inlist, **myargs)
        return
    elif ((outtype == 'short' and intype == 'ushort') or
          (outtype == 'ushort' and intype == 'short') or
          (outtype == 'real' and intype == 'integer') or
          (outtype == 'integer' and intype == 'real')):
        with pytest.warns(CCDProcWarning):
            iraf.ccdred.ccdproc(inlist, **myargs)
    else:
        iraf.ccdred.ccdproc(inlist, **myargs)

    # make sure all the outputs are giving the correct answer in the
    # correct type
    for ifile in outlists:
        with fits.open(ifile) as hdr:
            if outtype == 'short' or outtype == 'ushort':
                # nothing can get downscaled to these, so intype == outtype
                assert np.can_cast(hdr[0].data.dtype, dd[intype],
                                   casting='equiv')
                assert np.allclose(hdr[0].data, dd[intype](baseval - zeroval))
            elif outtype == 'integer':
                # all input integers varieties should be upscaled to "integer"
                if intype in ['short', 'ushort', 'integer']:
                    assert np.can_cast(hdr[0].data.dtype, dd['integer'],
                                       casting='equiv')
                    assert np.allclose(hdr[0].data,
                                       dd['integer'](baseval - zeroval))
                else:
                    # floats can't get downscaled to these, so intype == outtype
                    assert np.can_cast(hdr[0].data.dtype, dd[intype],
                                       casting='equiv')
                    assert np.allclose(hdr[0].data,
                                       dd[intype](baseval - zeroval))
            elif outtype == 'real':
                # all short integers varieties should be upscaled to "real"
                if intype in ['short', 'ushort', 'real']:
                    assert np.can_cast(hdr[0].data.dtype, dd['real'],
                                       casting='equiv')
                    assert np.allclose(hdr[0].data,
                                       dd['real'](baseval - zeroval))
                else:
                    # int and double can't get safely converted to real,
                    # so intype == outtype
                    assert np.can_cast(hdr[0].data.dtype, dd[intype],
                                       casting='equiv')
                    assert np.allclose(hdr[0].data,
                                       dd[intype](baseval - zeroval))
            else:
                # everything should be upscaled to "double"
                assert np.can_cast(hdr[0].data.dtype, dd['double'],
                                   casting='equiv')
                assert np.allclose(hdr[0].data,
                                   dd['double'](baseval - zeroval))


def test_ccdproc_set_sections(tmpdir):
    # XXX: do this
    # some of this might get handled/tested when doing all the other stuff,
    # (e.g. trim is probably tested in test_trim?). so maybe don't need to do
    # everything here.
    pass


def test_ccdproc_wcs():
    # XXX: do this at the bottom of setsections.x
    # set it to fail if the WCS stuff is already in the header?
    # test what it does?
    # otherwise, looks like I just need to put some stuff into
    # the output file header, can play around with what that looks like.
    # need to test it with different dimensions, etc

    # try making a full WCS header in astropy and feeding it through
    # while doing all processing and make sure none of the WCS things
    # change except for the lines added to the header.
    pass


def test_ccdproc_set_trim(tmpdir):
    basedir = str(tmpdir)

    # create some simple files for testing
    nx = 50
    ny = 90
    nimg = 3
    baseval = 100

    # create input files
    inputs = []
    arr = np.arange(nx * ny, dtype=float).reshape((nx, ny)) + baseval
    for jj in np.arange(nimg):
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

    # do my processing
    myargs = copy.deepcopy(defaultargs)
    myargs['output'] = outlists

    # which parts we want to do
    myargs['trim'] = True
    myargs['trimsec'] = None

    iraf.ccdred.ccdproc(inlist, **myargs)

    # trim still happens even if it doesn't actually do anything
    for ifile in outlists:
        with fits.open(ifile) as ff:
            assert 'trim' in ff[0].header
            assert ff[0].data.shape == arr.shape
            assert np.allclose(ff[0].data, arr)

    # add in the data section
    for jj in np.arange(nimg):
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['datasec'] = '[1:60,1:30]'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)

    iraf.ccdred.ccdproc(inlist, **myargs)

    # trim defaults to the data section
    for ifile in outlists:
        with fits.open(ifile) as ff:
            assert 'trim' in ff[0].header
            assert ff[0].data.shape == arr[:30, :60].shape
            assert np.allclose(ff[0].data, arr[:30, :60])

    # test passing in trimsec as an argument
    for jj in np.arange(nimg):
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['datasec'] = '[10:60,10:30]'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)

    # get a failure
    myargs['trimsec'] = '[5:65:2, 5:35:2]'
    with pytest.raises(CCDProcError):
        iraf.ccdred.ccdproc(inlist, **myargs)

    myargs['trimsec'] = '[5:650, 5:35]'
    with pytest.raises(CCDProcError):
        iraf.ccdred.ccdproc(inlist, **myargs)

    myargs['trimsec'] = '[5:65, 5:35]'
    iraf.ccdred.ccdproc(inlist, **myargs)

    # make sure we're using the input trim limits
    for ifile in outlists:
        with fits.open(ifile) as ff:
            assert 'trim' in ff[0].header
            assert ff[0].data.shape == arr[4:35, 4:65].shape
            assert np.allclose(ff[0].data, arr[4:35, 4:65])

    # return to the default trim, which should be "image"
    myargs = copy.deepcopy(defaultargs)
    myargs['output'] = outlists
    myargs['trim'] = True

    # test trimsec in header
    for jj in np.arange(nimg):
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['datasec'] = '[10:60,10:30]'
        hdu.header['trimsec'] = '[4:64, 4:34]'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)

    iraf.ccdred.ccdproc(inlist, **myargs)

    # make sure we're using the input trim limits
    for ifile in outlists:
        with fits.open(ifile) as ff:
            assert 'trim' in ff[0].header
            assert ff[0].data.shape == arr[3:34, 3:64].shape
            assert np.allclose(ff[0].data, arr[3:34, 3:64])

    # repeat just in case we change the trimsec default argument
    myargs['trimsec'] = 'image'
    iraf.ccdred.ccdproc(inlist, **myargs)

    # make sure we're using the input trim limits
    for ifile in outlists:
        with fits.open(ifile) as ff:
            assert 'trim' in ff[0].header
            assert ff[0].data.shape == arr[3:34, 3:64].shape
            assert np.allclose(ff[0].data, arr[3:34, 3:64])

    # test custom instruments
    for jj in np.arange(nimg):
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['mydata'] = '[10:60,10:30]'
        hdu.header['mytrim'] = '[6:64, 6:34]'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)

    inst = iraf.ccdred.Instrument()
    inst.parameters['trimsec'] = 'mytrim'
    inst.parameters['datasec'] = 'mydata'
    myargs['instrument'] = inst

    iraf.ccdred.ccdproc(inlist, **myargs)

    # make sure we're using the custom values
    for ifile in outlists:
        with fits.open(ifile) as ff:
            assert 'trim' in ff[0].header
            assert ff[0].data.shape == arr[5:34, 5:64].shape
            assert np.allclose(ff[0].data, arr[5:34, 5:64])

    # test trimsec inside of datasec
    for jj in np.arange(nimg):
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['datasec'] = '[4:64, 4:34]'
        hdu.header['trimsec'] = '[10:60,10:30]'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)

    myargs = copy.deepcopy(defaultargs)
    myargs['output'] = outlists
    myargs['trim'] = True
    iraf.ccdred.ccdproc(inlist, **myargs)

    # make sure we're using the input trim limits
    for ifile in outlists:
        with fits.open(ifile) as ff:
            assert 'trim' in ff[0].header
            assert ff[0].data.shape == arr[9:30, 9:60].shape
            assert np.allclose(ff[0].data, arr[9:30, 9:60])


def test_ccdproc_verbose_logfile():
    # XXX: do this
    pass


def test_ccdproc_set_fixpix():
    # XXX: do this
    pass


def test_ccdproc_set_overscan():
    # XXX: do this
    pass


def test_ccdproc_set_zero(tmpdir):
    basedir = str(tmpdir)

    # create some simple files for testing
    nx = 50
    ny = 90
    nimg = 3
    baseval = 100

    # create input files
    inputs = []
    for jj in np.arange(nimg):
        arr = np.arange(nx * ny, dtype=float).reshape((nx, ny)) + baseval
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

    zarr = np.arange(nx * ny, dtype=float).reshape((nx, ny)) + zeroval

    hdu = fits.PrimaryHDU(zarr)
    hdu.header['imagetyp'] = 'zero'
    hdu.writeto(zerofile, overwrite=True)

    # do my processing
    myargs = copy.deepcopy(defaultargs)
    myargs['output'] = outlists

    # which parts we want to do
    myargs['zerocor'] = False
    myargs['zero'] = zerofile

    iraf.ccdred.ccdproc(inlist, **myargs)

    # if no processing is desired, no file should be created.
    for ifile in outlists:
        assert not os.path.exists(ifile)

    # add zerocor to the header
    inputs = []
    for jj in np.arange(nimg):
        arr = np.arange(nx * ny, dtype=float).reshape((nx, ny)) + baseval
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['zerocor'] = 'foo'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)

    # which parts we want to do
    myargs['zerocor'] = True
    myargs['zero'] = zerofile

    iraf.ccdred.ccdproc(inlist, **myargs)

    # if zerocor is in the header, no processing should be done
    for ifile in outlists:
        assert not os.path.exists(ifile)

    inputs = []
    arr = np.arange(nx * ny, dtype=float).reshape((nx, ny)) + baseval
    for jj in np.arange(nimg):
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)

    iraf.ccdred.ccdproc(inlist, **myargs)

    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data, arr - zarr)
        assert len(hdr[0].header['zerocor']) > 0
        assert len(hdr[0].header['ccdproc']) > 0
        assert hdr[0].header['ccdsec'] == f'[1:{ny:d},1:{nx:d}]'
        hdr.close()

    # test datasec of input image as just a subset of the image
    dsec = '[5:80,10:40]'
    inputs = []
    for jj in np.arange(nimg):
        arr = np.arange(nx * ny, dtype=float).reshape((nx, ny)) + baseval
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['datasec'] = dsec
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)

    iraf.ccdred.ccdproc(inlist, **myargs)

    # image data sec [5:80, 10:40], CCD sec defaults to [1:76,1:31].
    # zero image data sec is whole image, [1:90,1:50]. zero ccd sec is same.
    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data[9:40, 4:80],
                           arr[9:40, 4:80] - zarr[:31, :76])
        data = hdr[0].data * 1
        data[9:40, 4:80] += zarr[:31, :76]
        assert np.allclose(data, arr)
        assert hdr[0].header['datasec'] == dsec
        assert len(hdr[0].header['zerocor']) > 0
        assert len(hdr[0].header['ccdproc']) > 0
        assert hdr[0].header['ccdsec'] == '[1:76,1:31]'
        hdr.close()

    # add the datasec/ccdsec parameter
    dsec = '[5:80,10:40]'
    ccdsec = '[4:79,11:41]'
    inputs = []
    for jj in np.arange(nimg):
        arr = np.arange(nx * ny, dtype=float).reshape((nx, ny)) + baseval
        hdu = fits.PrimaryHDU(arr)
        hdu.header['imagetyp'] = 'object'
        hdu.header['datasec'] = dsec
        hdu.header['ccdsec'] = ccdsec
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)

    iraf.ccdred.ccdproc(inlist, **myargs)

    # image data sec [5:80, 10:40], CCD sec now [4:79,11:41].
    # zero image data sec is whole image, [1:90,1:50]. zero ccd sec is same.
    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data[9:40, 4:80],
                           arr[9:40, 4:80] - zarr[10:41, 3:79])
        data = hdr[0].data * 1
        data[9:40, 4:80] += zarr[10:41, 3:79]
        assert np.allclose(data, arr)
        assert hdr[0].header['datasec'] == dsec
        assert len(hdr[0].header['zerocor']) > 0
        assert len(hdr[0].header['ccdproc']) > 0
        assert hdr[0].header['ccdsec'] == ccdsec
        hdr.close()

    # test zero image section failures
    zerofile = os.path.join(basedir, 'testzero.fits')
    zdsec = '[30:900, 5:45]'
    zccdsec = '[4:79,11:41]'

    hdu = fits.PrimaryHDU(zarr)
    hdu.header['imagetyp'] = 'zero'
    hdu.header['datasec'] = zdsec
    hdu.header['ccdsec'] = zccdsec
    hdu.writeto(zerofile, overwrite=True)
    # bad data section
    with pytest.raises(CCDProcError):
        iraf.ccdred.ccdproc(inlist, **myargs)

    zerofile = os.path.join(basedir, 'testzero.fits')
    zdsec = '[5:80,10:40]'
    zccdsec = '[4:79,11:42]'

    hdu = fits.PrimaryHDU(zarr)
    hdu.header['imagetyp'] = 'zero'
    hdu.header['datasec'] = zdsec
    hdu.header['ccdsec'] = zccdsec
    hdu.writeto(zerofile, overwrite=True)
    # ccd sec and data sec sizes don't agree
    with pytest.raises(CCDProcError):
        iraf.ccdred.ccdproc(inlist, **myargs)

    zerofile = os.path.join(basedir, 'testzero.fits')
    zdsec = '[5:80,10:40]'
    zccdsec = '[4:79,12:42]'

    hdu = fits.PrimaryHDU(zarr)
    hdu.header['imagetyp'] = 'zero'
    hdu.header['datasec'] = zdsec
    hdu.header['ccdsec'] = zccdsec
    hdu.writeto(zerofile, overwrite=True)
    # ccd sec of zero image doesn't cover that of input images
    with pytest.raises(CCDProcError):
        iraf.ccdred.ccdproc(inlist, **myargs)

    # test datasec/ccdsec matches of input and zero image
    zerofile = os.path.join(basedir, 'testzero.fits')
    zdsec = '[2:82,3:43]'
    zccdsec = '[3:83,10:50]'

    hdu = fits.PrimaryHDU(zarr)
    hdu.header['imagetyp'] = 'zero'
    hdu.header['datasec'] = zdsec
    hdu.header['ccdsec'] = zccdsec
    hdu.writeto(zerofile, overwrite=True)

    iraf.ccdred.ccdproc(inlist, **myargs)

    # image data sec [5:80, 10:40], CCD sec now [4:79,11:41].
    # zero image data sec is [2:82,3:43]. zero ccd sec is [3:80,10:70].
    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data[9:40, 4:80],
                           arr[9:40, 4:80] - zarr[3:34, 2:78])
        data = hdr[0].data * 1
        data[9:40, 4:80] += zarr[3:34, 2:78]
        assert np.allclose(data, arr)
        assert hdr[0].header['datasec'] == dsec
        assert len(hdr[0].header['zerocor']) > 0
        assert len(hdr[0].header['ccdproc']) > 0
        assert hdr[0].header['ccdsec'] == ccdsec
        hdr.close()

    # XXX: make the images closer to square to test the 1-D directions
    # check a 1-D zero image
    zarr = np.arange(ny, dtype=float) + zeroval
    zdsec = '[10:90]'
    zccdsec = '[1:81]'

    hdu = fits.PrimaryHDU(zarr)
    hdu.header['imagetyp'] = 'zero'
    hdu.header['datasec'] = zdsec
    hdu.header['ccdsec'] = zccdsec
    hdu.writeto(zerofile, overwrite=True)

    iraf.ccdred.ccdproc(inlist, **myargs)

    # image data sec [5:80, 10:40], CCD sec now [4:79,11:41].
    # zero image data sec is [10:90]. zero ccd sec is [1:81].
    for ifile in outlists:
        hdr = fits.open(ifile)
        assert np.allclose(hdr[0].data[9:40, 4:80],
                           arr[9:40, 4:80] - zarr[12:88])
        data = hdr[0].data * 1
        data[9:40, 4:80] += zarr[12:88]
        assert np.allclose(data, arr)
        assert hdr[0].header['datasec'] == dsec
        assert len(hdr[0].header['zerocor']) > 0
        assert len(hdr[0].header['ccdproc']) > 0
        assert hdr[0].header['ccdsec'] == ccdsec
        hdr.close()




    # XXX: can ccd sections be negative? what about data sections? in the zero
    #  image only?
    # XXX: need to test trim limits and making sure all ccd/data/trim sections
    # match up

# things to test:

# test having the calibration images in the input list.

# some tests to think about:
# 1D and 3D images
# 1 vs 0 indexing
# what about multiple FITS extensions, do we always require data to be in the
# first one?
# print notes/logs when doing recursive ccdproc?
# go through all arguments and make sure we have tests for all of them
# (e.g. verbose and logfile)
