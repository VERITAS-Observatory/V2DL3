import logging
import pytest

from pyV2DL3.eventdisplay.fillRESPONSE import check_fuzzy_boundary


def test_check_fuzzy_boundary(caplog):

    caplog.clear()
    par = 7.9
    boundary = 7.86
    tolerance_1 = 0.05
    tolerance_2 = 0.001

    check_fuzzy_boundary(par, boundary, tolerance_1)

    assert caplog.record_tuples == [("root", logging.WARNING, "Coordinate tolerance is 0.005 and is within 0.050")]

    with pytest.raises(ValueError):
        check_fuzzy_boundary(par, boundary, tolerance_2)


if __name__ == "__main__":
    test_check_fuzzy_boundary()
