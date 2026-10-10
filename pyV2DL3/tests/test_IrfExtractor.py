import numpy as np

from pyV2DL3.eventdisplay import IrfExtractor
from pyV2DL3.eventdisplay.IrfExtractor import find_closest_az, find_nearest


def test_find_nearest():
    az_centers = np.array(
        [
            187.5,
            206.25,
            232.5,
            255.0,
            277.5,
            300.0,
            322.5,
            345.0,
            7.5,
            30.0,
            52.5,
            75.0,
            97.5,
            120.0,
            142.5,
            161.25,
        ]
    )
    az_50 = find_nearest(az_centers, 50.0)
    az_188 = find_nearest(az_centers, 188.0)

    assert az_50 == 10 and az_188 == 0


def test_find_closest_az():
    # Include an all-azimuth sentinel and bins stored in their original ID order.
    azMins = np.array(
        [
            -1000.0,
            -180.0,
            -157.5,
            -135.0,
            -112.5,
            -90.0,
            -67.5,
            -45.0,
            -22.5,
            0.0,
            22.5,
            45.0,
            67.5,
            90.0,
            112.5,
            135.0,
            150.0,
        ]
    )
    azMaxs = np.array(
        [
            -165.0,
            -150.0,
            -120.0,
            -97.5,
            -75.0,
            -52.5,
            -30.0,
            -7.5,
            15.0,
            37.5,
            60.0,
            82.5,
            105.0,
            127.5,
            150.0,
            172.5,
            1000.0,
        ]
    )
    for azimuth in (359.0, 0.0, 360.0, -1.0):
        assert find_closest_az(azimuth, azMins, azMaxs) == 8
    assert find_closest_az(146.2, azMins, azMaxs) == 15
    assert find_closest_az(320.0, azMins, azMaxs) == 7
    assert find_closest_az(-180.0, azMins, azMaxs) == find_closest_az(180.0, azMins, azMaxs)


def test_extract_irf_accepts_zero_azimuth(monkeypatch):
    calls = {}

    def fake_extract_irf_2d(filename, irf_name, azimuth):
        calls.update(filename=filename, irf_name=irf_name, azimuth=azimuth)
        return "extracted"

    monkeypatch.setattr(IrfExtractor, "extract_irf_2d", fake_extract_irf_2d)

    assert IrfExtractor.extract_irf("effective_area.root", "eff", azimuth=0) == "extracted"
    assert calls == {
        "filename": "effective_area.root",
        "irf_name": "eff",
        "azimuth": 0,
    }


def test_azimuth_mask_preserves_ids_when_rows_are_shuffled():
    from pyV2DL3.eventdisplay.IrfExtractor import _get_az_mask

    class Branch:
        def __init__(self, values): self.values = np.array(values)
        def array(self, library): return self.values
    tree = {"az": Branch([11, 0, 11, 16]),
            "azMin": Branch([22.5, 135, 22.5, -1000]),
            "azMax": Branch([82.5, -165, 82.5, 1000])}
    assert _get_az_mask(50, tree).tolist() == [True, False, True, False]
