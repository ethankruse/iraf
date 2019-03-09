import os
import iraf
import pytest
from astropy.io import fits
from datetime import datetime
import numpy as np


def test_set_get_delete_header_value(tmpdir):
    basedir = str(tmpdir)
    fname = os.path.join(basedir, 'header_test.fits')

    inst = iraf.ccdred.Instrument()
    inst.parameters['subset'] = 'filter'
    inst.defaults['subset'] = 'default_filter'
    inst.defaults['exptime'] = -1
    # one we aren't going to put in the header
    inst.defaults['imagetyp'] = 'object'

    # set up a simple file with two headers/hdus
    hdu1 = fits.PrimaryHDU()
    hdu2 = fits.ImageHDU()
    new_hdul = fits.HDUList([hdu1, hdu2])
    new_hdul.writeto(fname, overwrite=True)

    # test adding with a normal key and one that needs to be translated.
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'red')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 1)
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset') == 'red'
        assert ff[0].header.comments['filter'] == ''
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime') == 1
        assert ff[0].header.comments['exptime'] == ''
        assert 'filter' not in ff[1].header and 'exptime' not in ff[1].header
        # test defaults
        assert iraf.ccdred.utils.get_header_value(ff, inst,
                                                  'imagetyp') == 'object'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'imagetyp',
                                                  default=True) == 'object'
        assert (iraf.ccdred.utils.get_header_value(ff, inst, 'subset',
                                                   default=True) ==
                'default_filter')
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime',
                                                  default=True) == -1
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'rdnoise') is None

    # same but adding a comment
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'blue',
                                           comment='c1')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 2, comment='c2')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset') == 'blue'
        assert ff[0].header.comments['filter'] == 'c1'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime') == 2
        assert ff[0].header.comments['exptime'] == 'c2'
        assert 'filter' not in ff[1].header and 'exptime' not in ff[1].header

    # test clearing the comment
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'blue',
                                           comment='')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 2, comment='')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset') == 'blue'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime') == 2
        assert ff[0].header.comments['filter'] == ''
        assert ff[0].header.comments['exptime'] == ''
        assert 'filter' not in ff[1].header and 'exptime' not in ff[1].header

    # test where the keys already exist in both headers.

    # add these to the second header
    with iraf.sys.image_open(fname, mode='update') as ff:
        ff[1].header['filter'] = 'uv'
        ff[1].header['exptime'] = 5
    # make sure we only change the first header
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'red',
                                           comment='com1')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 1,
                                           comment='com2')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert ff[0].header['filter'] == 'red'
        assert ff[0].header.comments['filter'] == 'com1'
        assert ff[0].header['exptime'] == 1
        assert ff[0].header.comments['exptime'] == 'com2'

        assert ff[1].header['filter'] == 'uv'
        assert ff[1].header.comments['filter'] == ''
        assert ff[1].header['exptime'] == 5
        assert ff[1].header.comments['exptime'] == ''

        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset') == 'red'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime') == 1

    # test where the keys only exist in the second header

    # remove them from the first header
    with iraf.sys.image_open(fname, mode='update') as ff:
        del ff[0].header['filter']
        del ff[0].header['exptime']

    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'blue',
                                           comment='c1')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 2,
                                           comment='c2')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset') == 'blue'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime') == 2
        assert ff[1].header.comments['filter'] == 'c1'
        assert ff[1].header.comments['exptime'] == 'c2'
        assert 'filter' not in ff[0].header and 'exptime' not in ff[0].header

    # test updating the value but not the comment
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'green')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 3)
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset') == 'green'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime') == 3
        assert ff[1].header.comments['filter'] == 'c1'
        assert ff[1].header.comments['exptime'] == 'c2'
        assert 'filter' not in ff[0].header and 'exptime' not in ff[0].header

    # test updating the comment but not the value
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', None,
                                           comment='c1mod')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', None,
                                           comment='c2mod')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset') == 'green'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime') == 3
        assert ff[1].header.comments['filter'] == 'c1mod'
        assert ff[1].header.comments['exptime'] == 'c2mod'
        assert 'filter' not in ff[0].header and 'exptime' not in ff[0].header

    # add these back to the first header
    with iraf.sys.image_open(fname, mode='update') as ff:
        ff[0].header['filter'] = 'uv'
        ff[0].header['exptime'] = 5

    # test deleting from headers
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.delete_header_value(ff, inst, 'subset')
        iraf.ccdred.utils.delete_header_value(ff, inst, 'exptime')
        # delete something that's not there
        iraf.ccdred.utils.delete_header_value(ff, inst, 'rdnoise')

    with iraf.sys.image_open(fname, mode='update') as ff:
        assert 'filter' not in ff[0].header and 'exptime' not in ff[0].header
        assert 'filter' not in ff[1].header and 'exptime' not in ff[1].header


def test_ccdsubset(tmpdir):
    basedir = str(tmpdir)

    inst = iraf.ccdred.Instrument()

    # create fake file
    fname = os.path.join(basedir, 'ccdsubsets.fits')
    hdu = fits.PrimaryHDU()
    hdu.header['subset'] = 'red filter'
    hdu.header['filter'] = 'blue filter! @#%_.asdf'
    hdu.writeto(fname, overwrite=True)

    with pytest.raises(Exception):
        with iraf.sys.image_open(fname, mode='update') as ff:
            iraf.ccdred.utils.ccdsubset(ff, 'notinstrument')

    # 'subset' as header value
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.ccdsubset(ff, inst) == 'red_filter'

    # translate 'subset'
    inst.parameters['subset'] = 'filter'
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.ccdsubset(ff, inst) == 'blue_filter______.asdf'

    # no default and not in header
    inst.parameters['subset'] = 'ifilt'
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.ccdsubset(ff, inst) == ''

    # check for default
    inst.defaults['subset'] = '  UB5.3  '
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.ccdsubset(ff, inst) == 'UB5.3'


def test_file_new_copy(tmpdir):
    basedir = str(tmpdir)

    inst = iraf.ccdred.Instrument()

    # create fake file
    fname = os.path.join(basedir, 'ccdnewcopy.fits')
    hdu = fits.PrimaryHDU()
    hdu.header['subset'] = 'red filter'
    hdu.header['filter'] = 'blue filter'
    hdu.writeto(fname, overwrite=True)

    outf = os.path.join(basedir, 'ccdnewcopy_out.fits')

    # test failures
    with iraf.sys.image_open(fname, mode='update') as ff:
        with pytest.raises(Exception):
            iraf.ccdred.utils.file_new_copy(outf, ff, mode='readonly',
                                            instrument=inst, overwrite=False)
        ff.__filetype__ = 'txt'
        with pytest.raises(Exception):
            iraf.ccdred.utils.file_new_copy(outf, ff, mode='NEW_COPY',
                                            instrument=inst, overwrite=False)

    # test simple example
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert not os.path.exists(outf)
        iraf.ccdred.utils.file_new_copy(outf, ff, mode='NEW_COPY',
                                        instrument=None, overwrite=False)
        assert os.path.exists(outf)

    with iraf.sys.image_open(outf, mode='update') as ff:
        assert 'origin' in ff[0].header
        assert 'date' in ff[0].header
        assert 'iraf-tlm' in ff[0].header

        # make sure these are both in UTC.
        dtime = datetime.strptime(ff[0].header['date'], '%Y-%m-%dT%H:%M:%S')
        tlmtime = datetime.strptime(ff[0].header['iraf-tlm'],
                                    '%Y-%m-%dT%H:%M:%S')
        assert abs((dtime - datetime.utcnow()).total_seconds()) < 2
        assert abs((tlmtime - datetime.utcnow()).total_seconds()) < 2

        # replace the date to test overwriting
        ff[0].header['date'] = 'none'

    # test overwrite
    with iraf.sys.image_open(fname, mode='update') as ff:
        with pytest.raises(OSError):
            iraf.ccdred.utils.file_new_copy(outf, ff, mode='NEW_COPY',
                                            instrument=inst, overwrite=False)
    with iraf.sys.image_open(outf, mode='update') as ff:
        assert ff[0].header['date'] == 'none'

    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.file_new_copy(outf, ff, mode='NEW_COPY',
                                        instrument=inst, overwrite=True)
    with iraf.sys.image_open(outf, mode='update') as ff:
        assert ff[0].header['date'] != 'none'

    # test using custom instrument
    inst.parameters['origin'] = 'orig'
    inst.parameters['date'] = 'dt'
    inst.parameters['iraf-tlm'] = 'tlm'
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.file_new_copy(outf, ff, mode='NEW_COPY',
                                        instrument=inst, overwrite=True)

    with iraf.sys.image_open(outf, mode='update') as ff:
        assert 'origin' not in ff[0].header
        assert 'date' not in ff[0].header
        assert 'iraf-tlm' not in ff[0].header

        assert 'orig' in ff[0].header
        assert 'dt' in ff[0].header
        assert 'tlm' in ff[0].header


def test_type_max():
    assert iraf.ccdred.utils.type_max(np.float64, np.float32) == np.float64
    assert iraf.ccdred.utils.type_max(np.float32, np.float64) == np.float64

    assert iraf.ccdred.utils.type_max(np.int32, np.float64) == np.float64
    assert iraf.ccdred.utils.type_max(np.int16, np.int32) == np.int32
    assert iraf.ccdred.utils.type_max(np.uint16, np.int32) == np.int32
    assert iraf.ccdred.utils.type_max(np.float32, np.complex64) == np.complex64

    with pytest.raises(Exception):
        assert iraf.ccdred.utils.type_max(np.uint32, np.int32) == np.int32
