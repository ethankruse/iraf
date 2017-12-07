import numpy as np
from astropy.io import fits
import iraf
import os
import copy
import pytest

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
               'instrument': None, 'logfile': None, 'verbose': False}

# left to test: sigmafile, ccdtype, subsets, project, offsets,
# masktype, maskvalue, scale, zero, weight, statsec, rdnoise,
# gain, snoise, sigscale, grow, logfile, verbose


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
def test_combine_basic(combine_dir):
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
    outim.close()


def test_reject_none_outtype(combine_dir):
    basedir = str(combine_dir)
    # iraf outtype values and the numpy equivalent (excluding long == int)
    dtypestr = ['short', 'ushort', 'integer', 'real', 'double']
    dtypes = [np.short, np.ushort, np.int_, np.single, np.double]

    # create some simple files for testing
    # use an even number to test the trickier median case
    nimg = 6
    nx = 20
    ny = 30
    methods = ['median', 'average']

    for ii, itype in enumerate(dtypes):
        inputs = []
        innums = []
        arr = None
        for jj in np.arange(nimg)+1:
            innum = jj**2
            # this makes the mean 15.66, which rounds up to make sure we're
            # handling int data types right
            if jj == 6:
                innum += 3
            innums.append(innum)
            arr = np.ones((nx, ny), dtype=itype) * innum
            hdu = fits.PrimaryHDU(arr)
            inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
            hdu.writeto(inim, overwrite=True)
            inputs.append(inim)
        innums = np.array(innums)
        # for the IRAF list
        inlist = os.path.join(basedir, 'infiles.txt')
        with open(inlist, 'w') as ff:
            for ifile in inputs:
                ff.write(ifile + '\n')
        iraflist = '@' + inlist
        outfile = os.path.join(basedir, 'testout_me.fits')

        for imethod in methods:
            myargs = copy.deepcopy(defaultargs)
            myargs['method'] = imethod
            myargs['outtype'] = dtypestr[ii]
            iraf.combine(iraflist, outfile, **myargs)

            outim = fits.open(outfile)
            # this is the simple average and median calculations we
            # should be doing
            if imethod == 'average':
                assert np.allclose(outim[0].data, itype(innums.mean()))
            else:
                assert np.allclose(outim[0].data, itype(np.median(innums)))
            # need to check name instead of dtype because of endianness problems
            # np.dtype('>f4') != np.dtype('<f4') but both have name == 'float32'
            assert arr.dtype.name == outim[0].data.dtype.name
            outim.close()


def test_reject_minmax(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    # use an even number to test nlow + nhigh == nimages
    nimg = 6
    nx = 20
    ny = 30
    inputs = []
    innums = []

    for jj in np.arange(nimg) + 1:
        innum = jj ** 2
        innums.append(innum)
        arr = np.ones((nx, ny), dtype=np.double) * innum
        hdu = fits.PrimaryHDU(arr)
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)
    innums = np.array(innums)
    # for the IRAF list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')
    iraflist = '@' + inlist
    outfile = os.path.join(basedir, 'testout_me.fits')

    myargs = copy.deepcopy(defaultargs)
    myargs['method'] = 'average'
    myargs['reject'] = 'minmax'

    for ilow in np.arange(nimg//2+1):
        for ihigh in np.arange(nimg//2+1):
            myargs['nlow'] = ilow
            myargs['nhigh'] = ihigh

            if ihigh + ilow < nimg:
                iraf.combine(iraflist, outfile, **myargs)

                outim = fits.open(outfile)
                imean = innums[ilow:nimg-ihigh].mean()
                assert np.allclose(outim[0].data, imean)
                outim.close()
            else:
                with pytest.raises(Exception):
                    iraf.combine(iraflist, outfile, **myargs)


def test_reject_pclip(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    nx = 20
    ny = 30
    innums = np.array([5, 5, 18, 19, 20, 21, 22, 23, 40, 40])
    pclips = [-2, 2]
    nkeeps = [2, -2]

    # test with an even and odd number of images for the median/pclip values
    for nimg in [9, 10]:
        inputs = []
        for jj in np.arange(nimg):
            arr = np.ones((nx, ny), dtype=np.double) * innums[jj]
            hdu = fits.PrimaryHDU(arr)
            inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
            hdu.writeto(inim, overwrite=True)
            inputs.append(inim)
        # for the IRAF list
        inlist = os.path.join(basedir, 'infiles.txt')
        with open(inlist, 'w') as ff:
            for ifile in inputs:
                ff.write(ifile + '\n')
        iraflist = '@' + inlist
        outfile = os.path.join(basedir, 'testout_me.fits')

        for iclip in pclips:
            for ikeep in nkeeps:
                myargs = copy.deepcopy(defaultargs)
                myargs['method'] = 'average'
                myargs['reject'] = 'pclip'
                myargs['nkeep'] = ikeep
                myargs['pclip'] = iclip
                myargs['lsigma'] = 3.
                myargs['hsigma'] = 3.

                iraf.combine(iraflist, outfile, **myargs)

                if ikeep > 0:
                    imean = innums[2:-2].mean()
                else:
                    imean = innums[:-2].mean()

                outim = fits.open(outfile)
                assert np.allclose(outim[0].data, imean)
                outim.close()


def test_reject_ccdclip_crreject(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    nx = 20
    ny = 30
    innums = np.array([3., 3., 30., 30.1, 30.2, 30.3, 30.4, 30.7])
    nimg = len(innums)
    inputs = []
    mclip = [True, False]
    rejects = ['ccdclip', 'crreject']
    nkeeps = [1, -1]
    # nm: the column order is readnoise, gain, snoise (defaults 0,1,0)
    # things to test: mclip, lthreshold, hthreshold
    # default sigma is just sqrt(median)

    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=np.double) * innums[jj]
        hdu = fits.PrimaryHDU(arr)
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)
    # for the IRAF list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')
    iraflist = '@' + inlist
    outfile = os.path.join(basedir, 'testout_me.fits')

    for iclip in mclip:
        for ireject in rejects:
            for ikeep in nkeeps:
                myargs = copy.deepcopy(defaultargs)
                myargs['method'] = 'average'
                myargs['reject'] = ireject
                myargs['nkeep'] = ikeep
                myargs['mclip'] = iclip
                # this tests mclip = False and how we calcluate the averages.
                # points will only be removed if the min/max values are
                # removed before calculating the average
                # if mclip = True, this has no effect, things are removed
                # anyway
                myargs['lsigma'] = 4.3

                iraf.combine(iraflist, outfile, **myargs)

                if ikeep < 0 or ireject == 'crreject':
                    imean = innums.mean()
                else:
                    imean = innums[2:].mean()

                outim = fits.open(outfile)
                assert np.allclose(outim[0].data, imean)
                outim.close()


def test_reject_sigclip(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    nx = 20
    ny = 30
    innums = np.array([3., 3., 30., 30.1, 30.2, 30.3, 30.4, 30.7])
    nimg = len(innums)
    inputs = []
    mclip = [True, False]
    nkeeps = [1, -1]

    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=np.double) * innums[jj]
        hdu = fits.PrimaryHDU(arr)
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)
    # for the IRAF list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')
    iraflist = '@' + inlist
    outfile = os.path.join(basedir, 'testout_me.fits')

    for iclip in mclip:
        for ikeep in nkeeps:
            myargs = copy.deepcopy(defaultargs)
            myargs['method'] = 'average'
            myargs['reject'] = 'sigclip'
            myargs['nkeep'] = ikeep
            myargs['mclip'] = iclip
            # need it this low since the sigma is so huge at first
            myargs['lsigma'] = 1.7

            iraf.combine(iraflist, outfile, **myargs)

            if ikeep < 0:
                imean = innums.mean()
            else:
                imean = innums[2:].mean()

            outim = fits.open(outfile)
            assert np.allclose(outim[0].data, imean)
            outim.close()


def test_reject_avsigclip(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    nx = 20
    ny = 30
    innums = np.array([3., 3., 30., 30.01, 30.02, 30.03, 30.04, 30.07])
    nimg = len(innums)
    inputs = []
    mclip = [True, False]
    nkeeps = [1, -1]

    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=np.double) * innums[jj]
        # make each row in the first dimension different
        arr += np.arange(nx)[:, np.newaxis]
        hdu = fits.PrimaryHDU(arr)
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)
    # for the IRAF list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')
    iraflist = '@' + inlist
    outfile = os.path.join(basedir, 'testout_me.fits')

    for iclip in mclip:
        for ikeep in nkeeps:
            myargs = copy.deepcopy(defaultargs)
            myargs['method'] = 'average'
            myargs['reject'] = 'avsigclip'
            myargs['nkeep'] = ikeep
            myargs['mclip'] = iclip
            # need it this low since the sigma is so huge at first
            myargs['lsigma'] = 1.85

            iraf.combine(iraflist, outfile, **myargs)

            if ikeep < 0:
                imean = innums.mean()
            else:
                imean = innums[2:].mean()

            outim = fits.open(outfile)
            # make sure we're doing the summing in the right dimension
            outarr = np.ones(ny) * imean + np.arange(nx)[:, np.newaxis]
            assert np.allclose(outim[0].data, outarr)
            outim.close()


def test_delete(combine_dir):
    basedir = str(combine_dir)
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
    myargs = copy.deepcopy(defaultargs)

    myargs['delete'] = False
    iraf.combine(iraflist, outfile, **myargs)
    assert os.path.exists(outfile)
    for ifile in inputs:
        assert os.path.exists(ifile)

    myargs['delete'] = True
    iraf.combine(iraflist, outfile, **myargs)
    assert os.path.exists(outfile)
    for ifile in inputs:
        assert not os.path.exists(ifile)


def test_threshold_blank(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    nx = 20
    ny = 30
    nimg = 6
    inputs = []
    blank = -3.
    lthresh = [2, 4]
    hthresh = [6, 3]

    for jj in np.arange(nimg):
        arr = np.ones((nx, ny), dtype=float) * (jj + 1)
        hdu = fits.PrimaryHDU(arr)
        inim = os.path.join(basedir, f'testimg{jj:02d}.fits')
        hdu.writeto(inim, overwrite=True)
        inputs.append(inim)
    # for the IRAF list
    inlist = os.path.join(basedir, 'infiles.txt')
    with open(inlist, 'w') as ff:
        for ifile in inputs:
            ff.write(ifile + '\n')
    iraflist = '@' + inlist
    outfile = os.path.join(basedir, 'testout_me.fits')

    myargs = copy.deepcopy(defaultargs)
    myargs['blank'] = blank

    nums = np.arange(nimg) + 1

    for lt in lthresh:
        for ht in hthresh:
            myargs['lthreshold'] = lt
            myargs['hthreshold'] = ht

            iraf.combine(iraflist, outfile, **myargs)

            outim = fits.open(outfile)

            valid = np.where((nums <= ht) & (nums >= lt))[0]
            if valid.size > 0:
                assert np.allclose(outim[0].data, nums[valid].mean())
            else:
                assert np.allclose(outim[0].data, blank)

            outim.close()


def test_plfile(combine_dir):
    basedir = str(combine_dir)
    # create some simple files for testing
    nx = 20
    ny = 30
    nimg = 5
    inputs = simple_inputs(nimg, nx, ny, basedir)

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
    myargs['reject'] = 'minmax'
    myargs['nhigh'] = 0

    nlows = np.arange(nimg)
    for ilow in nlows:
        myargs['nlow'] = ilow
        iraf.combine(iraflist, outfile, **myargs)
        outim = fits.open(plfile)
        assert np.allclose(outim[0].data, ilow)
        # need to check name instead of dtype because of endianness problems
        # np.dtype('>f4') != np.dtype('<f4') but both have name == 'float32'
        assert outim[0].data.dtype.name[:3] == 'int'
        outim.close()
