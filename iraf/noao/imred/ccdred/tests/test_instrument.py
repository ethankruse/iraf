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
