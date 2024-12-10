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
    az_bin_to_store_1 = find_closest_az(146.20, azMins, azMaxs)
    az_bin_to_store_2 = find_closest_az(-180.0, azMins, azMaxs)
    az_bin_to_store_3 = find_closest_az(320.0, azMins, azMaxs)
    assert az_bin_to_store_1 == 14 and az_bin_to_store_2 == 8 and az_bin_to_store_3 == 6


if __name__ == "__main__":
    test_find_nearest()
    test_find_closest_az()
