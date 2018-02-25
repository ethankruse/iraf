import iraf
import pytest
import os
from pathlib import Path


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
