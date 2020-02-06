import os
from pathlib import Path

import iraf


def test_is_iterable():
    # lists, tuples, and dicts are iterable
    assert iraf.utils.is_iterable(['foo', 'bar'])
    assert iraf.utils.is_iterable([])
    assert iraf.utils.is_iterable(('foo', 'bar'))
    assert iraf.utils.is_iterable({'foo': 0, 'bar': 1})
    # strings (and ints, bools, etc.) aren't iterable
    assert not iraf.utils.is_iterable('foo')
    assert not iraf.utils.is_iterable(0)
    assert not iraf.utils.is_iterable(True)


def test_file_handler(tmpdir):
    basedir = str(tmpdir)
    ltxt = 4
    lshorttxt = 3
    bftxt = 3
    ztwotxt = 2
    ifiles = ['a.txt', 'b.txt', 'c.txt', 'ff.txt', '1.fits', '2.fits', '3.fits']
    samplefiles = []
    subdir = os.path.join(basedir, 'sub')
    os.mkdir(subdir)
    # create the files and the same files in a subdirectory
    for ifile in ifiles:
        ipath = os.path.join(basedir, ifile)
        subpath = os.path.join(subdir, ifile)
        samplefiles.append(ipath)
        Path(ipath).touch()
        samplefiles.append(subpath)
        Path(subpath).touch()

    # test None and empty string
    assert len(iraf.utils.file_handler(None)) == 0
    assert len(iraf.utils.file_handler('')) == 0
    assert len(iraf.utils.file_handler('    ')) == 0

    # test of a file that doesn't exist
    inp = os.path.join(basedir, 'abc.txt')
    out = iraf.utils.file_handler(inp)
    assert len(out) == 0
    # test the exists flag
    out = iraf.utils.file_handler(inp, exists=False)
    assert len(out) == 1

    # find single files
    # input list of single files
    out = iraf.utils.file_handler(samplefiles[:4])
    assert len(out) == 4
    # string input of single file
    out = iraf.utils.file_handler(samplefiles[0])
    assert len(out) == 1

    # test of @lists
    atlist = os.path.join(basedir, 'atlist.list')
    with open(atlist, 'w') as ff:
        for sample in samplefiles:
            ff.write(sample + '\n')

    assert len(iraf.utils.file_handler('@'+atlist)) == len(samplefiles)

    # test multiple inputs and wildcards
    wild = [os.path.join(basedir, '*.txt'), os.path.join(basedir, '*.fits')]
    assert len(iraf.utils.file_handler(wild[0])) == ltxt
    assert len(iraf.utils.file_handler(wild)) == len(samplefiles)//2

    wild = os.path.join(basedir, '?.txt')
    assert len(iraf.utils.file_handler(wild)) == lshorttxt

    # test character ranges
    wild = os.path.join(basedir, '[b-f]*.txt')
    assert len(iraf.utils.file_handler(wild)) == bftxt
    wild = os.path.join(basedir, '[0-2].fits')
    assert len(iraf.utils.file_handler(wild)) == ztwotxt

    # test the recursive parameter
    wild = os.path.join(basedir, '**')
    # length of sample files in main dir + sub/ + atlist.list
    totmain = len(samplefiles)//2 + 2
    totsub = len(samplefiles)//2
    assert len(iraf.utils.file_handler(wild)) == totmain
    srch = iraf.utils.file_handler(wild, recursive=True)
    # plus 1 because with recursive on, glob also returns the main dir/
    # as an item for some reason
    assert len(srch) == totmain + totsub + 1

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

    assert len(iraf.utils.file_handler('@' + atlist)) == len(samplefiles)//2

    # check comments and blank lines in @lists
    with open(atlist, 'w') as ff:
        for sample in samplefiles:
            ff.write('#' + sample + '\n\n\n')

    assert len(iraf.utils.file_handler('@' + atlist)) == 0

    # check @lists of @lists
    txtlist = os.path.join(basedir, 'txtlist.list')
    inoexist = 3
    with open(txtlist, 'w') as ff:
        ff.write('hi.txt\ndoes_not_exist.txt\n\n')
        ff.write(os.path.join(basedir, 'sub'+os.sep+'*.txt'))

    reallist = os.path.join(basedir, 'reallist.list')
    wild = [os.path.join(basedir, '*.txt'), os.path.join(basedir, '*.fits')]
    jnoexist = 2
    with open(reallist, 'w') as ff:
        for iwild in wild:
            ff.write(iwild + '\n')

    with open(atlist, 'w') as ff:
        ff.write('\t@' + txtlist + '\n')
        ff.write('   @' + reallist + '\n')

    tot = ltxt + len(samplefiles)//2
    assert len(iraf.utils.file_handler('@' + atlist)) == tot
    # test the exists flag
    tot = inoexist + jnoexist
    assert len(iraf.utils.file_handler('@' + atlist, exists=False)) == tot
