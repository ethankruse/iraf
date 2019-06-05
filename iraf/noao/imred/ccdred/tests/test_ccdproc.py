import pytest
import os
import numpy as np
from astropy.io import fits
import copy
import time
import datetime
import iraf
from .. import ccdproc_routines as ccdr


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

    # test bad braces
    with pytest.raises(Exception):
        ccdr.ccd_section('[:,:')
    with pytest.raises(Exception):
        ccdr.ccd_section(':,:]')
    with pytest.raises(Exception):
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
    with pytest.raises(Exception):
        ccdr.ccd_section('[1:5]')
    with pytest.raises(Exception):
        ccdr.ccd_section('[1:5, 1:5, 1:5]')

    # test bad inputs in a dimension
    with pytest.raises(Exception):
        ccdr.ccd_section('[1:3:4:5, :]')

    # test bad size of defaults
    with pytest.raises(Exception):
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
        with pytest.raises(Exception):
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
    with pytest.raises(Exception):
        ccdr.cal_list([foo], listtype, inst, calimages, nscans, caltypes,
                      subsets, scantype, nscan, scancor)

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
        with pytest.raises(Exception):
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

    with pytest.raises(Exception):
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
    with pytest.raises(Exception):
        ccdr.cal_image(inopen, inst, imtype, 1, cals, True)

    # fail to reach cal_scan
    calim = os.path.join(basedir, f'testcal.{imtype}.fits')

    cals = ([calim], [6], [imtype], [''])
    with pytest.raises(Exception):
        ccdr.cal_image(inopen, inst, imtype, 1, cals, True)

    cals = ([calim], [5], [imtype], [''])
    with pytest.raises(Exception):
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
    cals = ([calim+'1', calim+'2'], [5, 5], [imtype, imtype], ['', ''])
    with pytest.raises(Exception):
        ccdr.cal_image(inopen, inst, imtype, 1, cals, True)

    # pick the right nscan
    cals = ([calim+'1', calim+'2'], [3, 2], [imtype, imtype], ['', ''])
    assert ccdr.cal_image(inopen, inst, imtype, 2, cals, True) == calim+'2'

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
    cals = ([calim+'1', calim+'2'], [1, 1], [imtype, imtype], ['blue', 'red'])
    if imtype in ['flat', 'illum', 'fringe']:
        cor = calim + '2'
        assert ccdr.cal_image(inopen, inst, imtype, 1, cals, False) == cor
    else:
        with pytest.raises(Exception):
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
    assert retstr[len(nowstr)+1:] == instr

    # check log files and verbose output
    log = os.path.join(basedir, f'log.txt')
    with open(log, 'w') as ff:
        retstr = ccdr.logstring(instr, inopen, True, ff)
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
    assert out[-len(instr)-1:-1] == instr

    cutstr = out[len(inim)+2:len(inim)+2+len(nowstr)]
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

    inopen.close()


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


def test_cal_open():
    # go through all ccdtypes.

    # make sure if an object is passed in as a list of a specific type it is
    # processed that way regardless of what's in header

    # make sure that passing in all files as one input list gets things sorted
    # into the right lists as well

    # check subsets are handled appropriately
    pass


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


# some tests to think about:
# 1D and 3D images
# 1 vs 0 indexing
# print notes/logs when doing recursive ccdproc?
