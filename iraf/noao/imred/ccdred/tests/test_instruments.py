import os
import iraf
import pytest
from astropy.io import fits
from datetime import datetime

parameters = {'BPM': None, 'biassec': None, 'ccdmean': None,
              'ccdmeant': None, 'ccdproc': None, 'ccdsec': None,
              'darkcor': None, 'darktime': None, 'datasec': None,
              'exptime': None, 'fixfile': None, 'fixpix': None,
              'flatcor': None, 'fringcor': None, 'gain': None,
              'illumcor': None, 'imagetyp': None, 'mkfringe': None,
              'mkillum': None, 'ncombine': None, 'nscanrow': None,
              'overscan': None, 'rdnoise': None, 'readcor': None,
              'snoise': None, 'subset': None, 'trim': None,
              'trimsec': None, 'zerocor': None, 'origin': None,
              'date': None, 'iraf-tlm': None}

image_types = "object|zero|dark|flat|illum|fringe|other|comp".split('|')


def test_instrument_init(tmpdir):
    # default instrument
    inst = iraf.ccdred.Instrument()
    # test initialized values are right
    for key in inst.parameters:
        assert inst.translate(key) == key and inst.get_default(key) is None
    # make sure this testing list is up to date
    for key in parameters:
        assert key in inst.parameters
    for key in inst.parameters:
        assert key in parameters

    assert inst.translation_file is None
    assert len(inst.image_types) == 0


def test_instrument_translation_file(tmpdir):
    basedir = str(tmpdir)
    # create a translation file
    tfile = os.path.join(basedir, 'trans.txt')

    # test the ccd types
    customs = ["DOME FLAT", "ARC LAMP (7)", "Comp star"]
    with open(tfile, 'w') as ff:
        ff.write('# this is a comment line here\n\n     \t\t \t      \n')
        # test each parameter with a custom parameter name, no default values
        for ipar in parameters:
            cs = ipar + 'custom'
            ff.write(f'{ipar}   \t \t  \t  {cs}\n')
        for ii in range(len(customs)):
            for jj in range(len(image_types)):
                bstr = customs[ii]+image_types[jj]
                # test each image type with different strings
                ff.write(f"'{bstr}'\t{image_types[jj]}\n")
    inst = iraf.ccdred.Instrument(tfile)
    # make sure they got placed properly
    for key in inst.parameters:
        assert (inst.parameters[key] == key + 'custom' and
                inst.defaults[key] is None)
    for ii in range(len(customs)):
        for jj in range(len(image_types)):
            bstr = customs[ii]+image_types[jj]
            assert inst.image_types[bstr] == image_types[jj]
    assert inst.translation_file == tfile
    assert len(inst.image_types) == len(customs) * len(image_types)

    # repeat but with default values this time and spaces in custom parameter
    # names and default values
    tfile = os.path.join(basedir, 'trans2.txt')
    with open(tfile, 'w') as ff:
        ff.write('# this is a comment line here\n\n\n')
        for ipar in parameters:
            cs = ipar + ' custom'
            # wrap in double quotes to test both quote types in one line
            cs = f'"{cs}"'
            ff.write(f"{ipar}     {cs}\t'{ipar} default'\n")
    inst = iraf.ccdred.Instrument(tfile)
    for key in inst.parameters:
        assert (inst.parameters[key] == key + ' custom' and
                inst.defaults[key] == key + ' default')
    assert inst.translation_file == tfile

    # break it with wrong numbers of values per row
    with open(tfile, 'w') as ff:
        ff.write('subset \n')
    with pytest.raises(Exception):
        iraf.ccdred.Instrument(tfile)
    with open(tfile, 'w') as ff:
        ff.write('subset filter default bogus\n')
    with pytest.raises(Exception):
        iraf.ccdred.Instrument(tfile)

    # invalid parameter name
    with open(tfile, 'w') as ff:
        ff.write('subsets filter \n')
    with pytest.raises(Exception):
        iraf.ccdred.Instrument(tfile)

    # nonsense in both/3 values
    with open(tfile, 'w') as ff:
        ff.write('foo filter default \n')
    with pytest.raises(Exception):
        iraf.ccdred.Instrument(tfile)

    # try to add a custom image type with a 'default value'
    with open(tfile, 'w') as ff:
        ff.write('foo dark default \n')
    with pytest.raises(Exception):
        iraf.ccdred.Instrument(tfile)

    # check comments at end of line
    with open(tfile, 'w') as ff:
        ff.write('subset filter default # this is a comment \n')
    inst = iraf.ccdred.Instrument(tfile)
    assert (inst.parameters['subset'] == 'filter' and
            inst.defaults['subset'] == 'default')

    # test repeated values
    with open(tfile, 'w') as ff:
        ff.write('subset filter default\nsubset new def')
    inst = iraf.ccdred.Instrument(tfile)
    assert (inst.parameters['subset'] == 'new' and
            inst.defaults['subset'] == 'def')


# test the three methods of the class
def test_instrument_translate_get_default(tmpdir):
    basedir = str(tmpdir)
    # create a translation file
    tfile = os.path.join(basedir, 'trans.txt')
    with open(tfile, 'w') as ff:
        ff.write('subset filter default\nrdnoise noise\n')
    inst = iraf.ccdred.Instrument(tfile)
    assert inst.translate('subset') == 'filter'
    assert inst.translate('rdnoise') == 'noise'
    assert inst.translate('imagetyp') == 'imagetyp'
    with pytest.raises(KeyError):
        inst.translate('boo')

    assert inst.get_default('subset') == 'default'
    assert inst.get_default('rdnoise') is None
    assert inst.get_default('imagetyp') is None
    with pytest.raises(KeyError):
        inst.get_default('boo')


def test_instrument_get_image_type_ccdtypes(tmpdir):
    basedir = str(tmpdir)
    # create a translation file
    tfile = os.path.join(basedir, 'trans.txt')
    with open(tfile, 'w') as ff:
        ff.write('"sky flat"  flat\n')
    inst = iraf.ccdred.Instrument(tfile)
    assert inst.get_image_type("sky flat") == 'flat'
    assert inst.get_image_type("object") == 'object'
    assert inst.get_image_type("foo") == 'unknown'
    assert inst.get_image_type(None) == 'none'
    assert inst.get_image_type('') == 'none'
    assert inst.get_image_type('\t\n') == 'none'

    # create fake file
    fname = os.path.join(basedir, 'ccdtypestest.fits')
    hdu = fits.PrimaryHDU()
    hdu.header['imagetyp'] = 'object'
    hdu.header['itype'] = 'sky flat'
    hdu.header['itypo'] = 'foo'
    hdu.writeto(fname, overwrite=True)

    with pytest.raises(Exception):
        with iraf.sys.image_open(fname) as ff:
            iraf.ccdred.utils.ccdtypes(ff, 'notinstrument')

    # check ccdtypes
    with iraf.sys.image_open(fname) as ff:
        assert iraf.ccdred.utils.ccdtypes(ff, inst) == 'object'

    inst.parameters['imagetyp'] = 'itype'
    with iraf.sys.image_open(fname) as ff:
        assert iraf.ccdred.utils.ccdtypes(ff, inst) == 'flat'

    # check for correct default value when not in image header
    inst.parameters['imagetyp'] = 'imtype'
    inst.defaults['imagetyp'] = 'dark'
    with iraf.sys.image_open(fname) as ff:
        assert iraf.ccdred.utils.ccdtypes(ff, inst) == 'dark'

    # check for none when not in image header and no default
    inst.defaults['imagetyp'] = None
    with iraf.sys.image_open(fname) as ff:
        assert iraf.ccdred.utils.ccdtypes(ff, inst) == 'none'

    # check for 'unknown' value when not recognized image type
    inst.parameters['imagetyp'] = 'itypo'
    inst.defaults['imagetyp'] = 'dark'
    with iraf.sys.image_open(fname) as ff:
        assert iraf.ccdred.utils.ccdtypes(ff, inst) == 'unknown'


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
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'imagetyp') == 'object'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'imagetyp',
                                            default=True) == 'object'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset',
                                            default=True) == 'default_filter'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime',
                                            default=True) == -1
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'rdnoise') is None

    # same but adding a comment
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'blue', comment='c1')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 2, comment='c2')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'subset') == 'blue'
        assert ff[0].header.comments['filter'] == 'c1'
        assert iraf.ccdred.utils.get_header_value(ff, inst, 'exptime') == 2
        assert ff[0].header.comments['exptime'] == 'c2'
        assert 'filter' not in ff[1].header and 'exptime' not in ff[1].header

    # test clearing the comment
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'blue', comment='')
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
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'red', comment='com1')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 1, comment='com2')
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
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', 'blue', comment='c1')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', 2, comment='c2')
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
        iraf.ccdred.utils.set_header_value(ff, inst, 'subset', None, comment='c1mod')
        iraf.ccdred.utils.set_header_value(ff, inst, 'exptime', None, comment='c2mod')
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
