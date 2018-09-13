import iraf
import os, copy

parameters = {'BPM': None, 'biassec': None, 'ccdmean': None,
              'ccdmeant': None, 'ccdproc': None, 'ccdsec': None,
              'darkcor': None,
              'darktime': None, 'datasec': None, 'exptime': None,
              'fixfile': None, 'fixpix': None, 'flatcor': None,
              'fringcor': None, 'gain': None, 'illumcor': None,
              'imagetyp': None, 'mkfringe': None, 'mkillum': None,
              'ncombine': None, 'nscanrow': None, 'overscan': None,
              'rdnoise': None,
              'readcor': None, 'snoise': None, 'subset': None,
              'trim': None, 'trimsec': None, 'zerocor': None,
              'origin': None, 'date': None, 'iraf-tlm': None}
image_types = "object|zero|dark|flat|illum|fringe|other|comp".split('|')


def test_instrument_init(tmpdir):
    basedir = str(tmpdir)

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

    imtyps = copy.deepcopy(image_types)
    imtyps.append('blah')
    for it in imtyps:
        assert inst.get_image_type(it) == it
    assert inst.translation_file is None
    assert len(inst.image_types) == 0

    # create a translation file
    tfile = os.path.join(basedir, 'trans.txt')
    customs = ["DOME FLAT", "ARC LAMP", "Comp star"]
    bases = ["zero", "comp", "object"]
    with open(tfile, 'w') as ff:
        ff.write('# this is a comment line here\n\n\n')
        for ipar in parameters:
            cs = ipar + 'custom'
            ff.write(f'{ipar}     {cs}\n')
        for ii in range(len(customs)):
            ff.write(f"'{customs[ii]}'   {bases[ii]}\n")

    inst = iraf.Instrument(tfile)
    for key in inst.parameters:
        assert (inst.translate(key) == key+'custom' and
                inst.get_default(key) is None)
    for ii in range(len(customs)):
        assert inst.get_image_type(customs[ii]) == bases[ii]
    assert inst.translation_file == tfile
    assert len(inst.image_types) == len(customs)

    # repeat but with default values this time and spaces in custom parameter
    # names
    tfile = os.path.join(basedir, 'trans2.txt')
    customs = ["DOME FLAT", "ARC LAMP", "Comp star"]
    bases = ["zero", "comp", "object"]
    with open(tfile, 'w') as ff:
        ff.write('# this is a comment line here\n\n\n')
        for ipar in parameters:
            cs = ipar + ' custom'
            ff.write(f"{ipar}     '{cs}'  '{ipar} default'\n")
        for ii in range(len(customs)):
            ff.write(f"'{customs[ii]}'   {bases[ii]}\n")

    inst = iraf.Instrument(tfile)
    for key in inst.parameters:
        assert (inst.translate(key) == key + ' custom' and
                inst.get_default(key) == key + ' default')
    for ii in range(len(customs)):
        assert inst.get_image_type(customs[ii]) == bases[ii]
    assert inst.translation_file == tfile
    assert len(inst.image_types) == len(customs)

    # tests for bad instrument file construction (tabs, double quotes, wrong
    # parameters, wrong numbers of items on a line)


# test the three methods of the class
def test_new():
    pass
