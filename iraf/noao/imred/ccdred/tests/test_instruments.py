import os
import iraf
import pytest
from astropy.io import fits

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

# XXX: test for return value types in the 3 methods. what if someone sets
# defaults or parameters manually and uses e.g. an int for ncombine? Will that
# break things in other places? Should they always be strings?


def test_instrument_init(tmpdir):
    # default instrument
    inst = iraf.Instrument()
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
    inst = iraf.Instrument(tfile)
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
    inst = iraf.Instrument(tfile)
    for key in inst.parameters:
        assert (inst.parameters[key] == key + ' custom' and
                inst.defaults[key] == key + ' default')
    assert inst.translation_file == tfile

    # break it with wrong numbers of values per row
    with open(tfile, 'w') as ff:
        ff.write('subset \n')
    with pytest.raises(Exception):
        iraf.Instrument(tfile)
    with open(tfile, 'w') as ff:
        ff.write('subset filter default bogus\n')
    with pytest.raises(Exception):
        iraf.Instrument(tfile)

    # invalid parameter name
    with open(tfile, 'w') as ff:
        ff.write('subsets filter \n')
    with pytest.raises(Exception):
        iraf.Instrument(tfile)

    # nonsense in both/3 values
    with open(tfile, 'w') as ff:
        ff.write('foo filter default \n')
    with pytest.raises(Exception):
        iraf.Instrument(tfile)

    # try to add a custom image type with a 'default value'
    with open(tfile, 'w') as ff:
        ff.write('foo dark default \n')
    with pytest.raises(Exception):
        iraf.Instrument(tfile)

    # check comments at end of line
    with open(tfile, 'w') as ff:
        ff.write('subset filter default # this is a comment \n')
    inst = iraf.Instrument(tfile)
    assert (inst.parameters['subset'] == 'filter' and
            inst.defaults['subset'] == 'default')

    # test repeated values
    with open(tfile, 'w') as ff:
        ff.write('subset filter default\nsubset new def')
    inst = iraf.Instrument(tfile)
    assert (inst.parameters['subset'] == 'new' and
            inst.defaults['subset'] == 'def')


# test the three methods of the class
def test_instrument_translate_get_default(tmpdir):
    basedir = str(tmpdir)
    # create a translation file
    tfile = os.path.join(basedir, 'trans.txt')
    with open(tfile, 'w') as ff:
        ff.write('subset filter default\nrdnoise noise\n')
    inst = iraf.Instrument(tfile)
    assert inst.translate('subset') == 'filter'
    assert inst.translate('rdnoise') == 'noise'
    assert inst.translate('imagetyp') == 'imagetyp'
    assert inst.translate('boo') == 'boo'

    assert inst.get_default('subset') == 'default'
    assert inst.get_default('rdnoise') is None
    assert inst.get_default('imagetyp') is None
    assert inst.get_default('boo') is None


def test_instrument_get_image_type(tmpdir):
    basedir = str(tmpdir)
    # create a translation file
    tfile = os.path.join(basedir, 'trans.txt')
    with open(tfile, 'w') as ff:
        ff.write('"sky flat"  flat\n')
    inst = iraf.Instrument(tfile)
    assert inst.get_image_type("sky flat") == 'flat'
    assert inst.get_image_type("object") == 'object'
    assert inst.get_image_type("foo") == 'foo'


def test_set_header_value(tmpdir):
    basedir = str(tmpdir)
    fname = os.path.join(basedir, 'header_test.fits')

    inst = iraf.Instrument()
    inst.parameters['subset'] = 'filter'

    # set up a simple file with two headers/hdus
    hdu1 = fits.PrimaryHDU()
    hdu2 = fits.ImageHDU()
    new_hdul = fits.HDUList([hdu1, hdu2])
    new_hdul.writeto(fname, overwrite=True)

    # test adding with a normal key and one that needs to be translated.
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.set_header_value(ff, inst, 'subset', 'red')
        iraf.set_header_value(ff, inst, 'exptime', 1)
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert ff[0].header['filter'] == 'red'
        assert ff[0].header.comments['filter'] == ''
        assert ff[0].header['exptime'] == 1
        assert ff[0].header.comments['exptime'] == ''
        assert 'filter' not in ff[1].header and 'exptime' not in ff[1].header

    # same but adding a comment
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.set_header_value(ff, inst, 'subset', 'blue', comment='c1')
        iraf.set_header_value(ff, inst, 'exptime', 2, comment='c2')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert ff[0].header['filter'] == 'blue'
        assert ff[0].header.comments['filter'] == 'c1'
        assert ff[0].header['exptime'] == 2
        assert ff[0].header.comments['exptime'] == 'c2'
        assert 'filter' not in ff[1].header and 'exptime' not in ff[1].header

    # test clearing the comment
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.set_header_value(ff, inst, 'subset', 'blue', comment='')
        iraf.set_header_value(ff, inst, 'exptime', 2, comment='')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert ff[0].header['filter'] == 'blue'
        assert ff[0].header.comments['filter'] == ''
        assert ff[0].header['exptime'] == 2
        assert ff[0].header.comments['exptime'] == ''
        assert 'filter' not in ff[1].header and 'exptime' not in ff[1].header

    # test where the keys already exist both headers.

    # add these to the second header
    with iraf.sys.image_open(fname, mode='update') as ff:
        ff[1].header['filter'] = 'uv'
        ff[1].header['exptime'] = 5
    # make sure we only change the first header
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.set_header_value(ff, inst, 'subset', 'red', comment='com1')
        iraf.set_header_value(ff, inst, 'exptime', 1, comment='com2')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert ff[0].header['filter'] == 'red'
        assert ff[0].header.comments['filter'] == 'com1'
        assert ff[0].header['exptime'] == 1
        assert ff[0].header.comments['exptime'] == 'com2'

        assert ff[1].header['filter'] == 'uv'
        assert ff[1].header.comments['filter'] == ''
        assert ff[1].header['exptime'] == 5
        assert ff[1].header.comments['exptime'] == ''

    # test where the keys only exist in the second header

    # remove them from the first header
    with iraf.sys.image_open(fname, mode='update') as ff:
        del ff[0].header['filter']
        del ff[0].header['exptime']

    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.set_header_value(ff, inst, 'subset', 'blue', comment='c1')
        iraf.set_header_value(ff, inst, 'exptime', 2, comment='c2')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert ff[1].header['filter'] == 'blue'
        assert ff[1].header.comments['filter'] == 'c1'
        assert ff[1].header['exptime'] == 2
        assert ff[1].header.comments['exptime'] == 'c2'
        assert 'filter' not in ff[0].header and 'exptime' not in ff[0].header

    # test updating the value but not the comment
    with iraf.sys.image_open(fname, mode='update') as ff:
        iraf.set_header_value(ff, inst, 'subset', 'green')
        iraf.set_header_value(ff, inst, 'exptime', 3)
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert ff[1].header['filter'] == 'green'
        assert ff[1].header.comments['filter'] == 'c1'
        assert ff[1].header['exptime'] == 3
        assert ff[1].header.comments['exptime'] == 'c2'
        assert 'filter' not in ff[0].header and 'exptime' not in ff[0].header

    # test updating the comment but not the value
        iraf.set_header_value(ff, inst, 'subset', None, comment='c1mod')
        iraf.set_header_value(ff, inst, 'exptime', None, comment='c2mod')
    with iraf.sys.image_open(fname, mode='update') as ff:
        assert ff[1].header['filter'] == 'green'
        assert ff[1].header.comments['filter'] == 'c1mod'
        assert ff[1].header['exptime'] == 3
        assert ff[1].header.comments['exptime'] == 'c2mod'
        assert 'filter' not in ff[0].header and 'exptime' not in ff[0].header
