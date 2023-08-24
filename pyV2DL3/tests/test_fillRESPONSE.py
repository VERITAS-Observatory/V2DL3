import pytest

from pyV2DL3.eventdisplay.fillRESPONSE import check_fuzzy_boundary


def test_check_fuzzy_boundary(caplog):

    caplog.clear()
    par_name = "pedvar"
    par = 7.9
    boundary = 7.86
    tolerance_1 = 0.05
    tolerance_2 = 0.001

    check_fuzzy_boundary(par, boundary, tolerance_1, par_name)

    assert check_fuzzy_boundary(par, boundary, tolerance_1, par_name) == 1

    with pytest.raises(ValueError):
        check_fuzzy_boundary(par, boundary, tolerance_2, par_name)

    boundary = 0.0
    assert check_fuzzy_boundary(par, boundary, tolerance_1, par_name) == 0

    boundary = -1.0
    assert check_fuzzy_boundary(par, boundary, tolerance_1, par_name) == 0


if __name__ == "__main__":
    test_check_fuzzy_boundary()
