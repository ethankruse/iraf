import iraf
import pytest


def test_ccd_section():
    # get the defaults
    x0, x1, xs, y0, y1, ys = iraf.ccd_section('[:,:]')
    assert (x0 == 0 and x1 is None and xs == 1 and y0 == 0 and y1 is None and
            ys == 1)

    x0, x1, xs, y0, y1, ys = iraf.ccd_section(None)
    assert (x0 == 0 and x1 is None and xs == 1 and y0 == 0 and y1 is None and
            ys == 1)

    x0, x1, xs, y0, y1, ys = iraf.ccd_section('  ')
    assert (x0 == 0 and x1 is None and xs == 1 and y0 == 0 and y1 is None and
            ys == 1)

    # all the options
    x0, x1, xs, y0, y1, ys = iraf.ccd_section('[1:5:3,5:20:-4]')
    assert (x0 == 1 and x1 == 5 and xs == 3 and y0 == 5 and y1 == 20 and
            ys == -4)

    # no step size
    x0, x1, xs, y0, y1, ys = iraf.ccd_section('[1:5,5:20]')
    assert (x0 == 1 and x1 == 5 and xs == 1 and y0 == 5 and y1 == 20 and
            ys == 1)

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
