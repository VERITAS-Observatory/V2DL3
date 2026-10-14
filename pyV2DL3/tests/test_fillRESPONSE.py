import numpy as np
import pytest
from astropy.io import fits

from pyV2DL3.eventdisplay.fillRESPONSE import (
    check_fuzzy_boundary,
    check_parameter_range,
    fill_effective_area,
    fill_energy_migration,
    find_camera_offsets,
)


@pytest.mark.parametrize("offsets", [[0.5], [0.0, 0.5, 1.0]])
def test_response_offset_dimensions(offsets, tmp_path):
    """Serialized responses have one curve per declared offset bin."""
    offsets = np.asarray(offsets)
    theta_low, theta_high = find_camera_offsets(offsets)

    class Interpolator:
        def set_irf(self, name, **kwargs):
            self.name = name

        def interpolate(self, coordinate):
            energy = np.linspace(-1, 1, 60)
            if self.name == "eff":
                return np.full(60, 100 + coordinate[2]), [energy]
            return np.ones((3, 60)), [energy, np.array([0.5, 1.0, 1.5])]

    interpolator = Interpolator()
    ea, _, _ = fill_effective_area(
        "eff", interpolator, offsets, 1, 20, theta_low, theta_high
    )
    migration = fill_energy_migration(
        "hEsysMCRelative2D", interpolator, offsets, 1, 20, theta_low, theta_high
    )
    path = tmp_path / "responses.fits"
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(ea, name="EFFECTIVE AREA"),
        fits.BinTableHDU(migration, name="ENERGY DISPERSION"),
    ]).writeto(path)
    with fits.open(path) as hdus:
        area = hdus[1].data[0]
        dispersion = hdus[2].data[0]
        assert area["EFFAREA"].shape == (len(offsets), 60)
        assert dispersion["MATRIX"].shape == (len(offsets), 3, 60)
        for response in (area, dispersion):
            assert response["THETA_LO"].shape == (len(offsets),)
            assert response["THETA_HI"].shape == (len(offsets),)
        np.testing.assert_allclose(area["EFFAREA"][:, 0], 100 + offsets)
        if len(offsets) == 1:
            np.testing.assert_array_equal(area["THETA_LO"], [0.0])
            np.testing.assert_array_equal(area["THETA_HI"], [1.0])


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
        fuzzy_boundary=0.0,
    )

    assert result == 20.0
    assert "using the lower boundary" in caplog.text


def test_camera_offset_centers_are_converted_to_edges():
    from pyV2DL3.eventdisplay.fillRESPONSE import find_camera_offsets

    theta_low, theta_high = find_camera_offsets(np.array([0.5, 1.0]))
    assert np.array_equal(theta_low, [0.25, 0.75])
    assert np.array_equal(theta_high, [0.75, 1.25])

    theta_low, theta_high = find_camera_offsets(np.array([0.5]))
    assert np.array_equal(theta_low, [0.0])
    assert np.array_equal(theta_high, [1.0])

    theta_low, theta_high = find_camera_offsets(np.array([0.0]))
    assert np.array_equal(theta_low, [0.0])
    assert np.array_equal(theta_high, [0.5])


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
