import numpy as np
import pytest
from astropy.io import fits

from pyV2DL3.eventdisplay.fillRESPONSE import (
    check_fuzzy_boundary,
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
            np.testing.assert_array_equal(area["THETA_HI"], [10.0])


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
