import numpy as np

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
    # Paired intervals in the producer's original ID order, including wrapped bins.
    mins = np.array([135, 150, -180] + [-180 + 22.5*i for i in range(1, 14)] + [-1000])
    maxs = np.array([-165, -150, -120] + [-120 + 22.5*i for i in range(1, 14)] + [1000])
    for azimuth, expected in [(1, 9), (50, 11), (165, 0), (180, 1), (359, 9), (0, 9), (360, 9)]:
        assert find_closest_az(azimuth, mins, maxs) == expected
    centres = (mins[:-1] + (maxs[:-1] - mins[:-1]) % 360 / 2) % 360
    for index, centre in enumerate(centres):
        assert find_closest_az(centre, mins, maxs) == index


def test_azimuth_mask_preserves_ids_when_rows_are_shuffled():
    from pyV2DL3.eventdisplay.IrfExtractor import _get_az_mask
    class Branch:
        def __init__(self, values): self.values = np.array(values)
        def array(self, library): return self.values
    tree = {"az": Branch([11, 0, 11, 16]),
            "azMin": Branch([22.5, 135, 22.5, -1000]),
            "azMax": Branch([82.5, -165, 82.5, 1000])}
    assert _get_az_mask(50, tree).tolist() == [True, False, True, False]
