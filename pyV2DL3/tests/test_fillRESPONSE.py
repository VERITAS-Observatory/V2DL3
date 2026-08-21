import numpy as np
import pytest

from pyV2DL3.eventdisplay.fillRESPONSE import (
    check_fuzzy_boundary,
    check_parameter_range,
)


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


def test_low_zenith_uses_lower_irf_boundary(caplog):
    caplog.clear()

    result = check_parameter_range(
        12.0,
        np.array([20.0, 40.0, 60.0]),
        "zenith",
        use_click=False,
        fuzzy_boundary=0.05,
    )

    assert result == 20.0
    assert "using the lower boundary" in caplog.text


def test_cli_style_fuzzy_boundary_is_selected_per_axis():
    result = check_parameter_range(
        12.0,
        np.array([20.0, 40.0, 60.0]),
        "zenith",
        use_click=False,
        fuzzy_boundary=(("zenith", 0.05),),
    )

    assert result == 20.0


def test_high_zenith_still_requires_fuzzy_boundary():
    with pytest.raises(ValueError):
        check_parameter_range(
            65.0,
            np.array([20.0, 40.0, 60.0]),
            "zenith",
            use_click=False,
            fuzzy_boundary=0.05,
        )


def test_pedvar_lower_boundary_remains_strict():
    with pytest.raises(ValueError):
        check_parameter_range(
            2.0,
            np.array([5.0, 7.0]),
            "pedvar",
            use_click=False,
            fuzzy_boundary=0.05,
        )


if __name__ == "__main__":
    test_check_fuzzy_boundary()
