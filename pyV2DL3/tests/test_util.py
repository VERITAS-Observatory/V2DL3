import numpy as np

from pyV2DL3.eventdisplay.util import getGTI, getRunQuality


def test_getRunQuality():
    vpm_messages = [
        "(VPM) data from database for telescope 1",
        "(VPM) data from database for telescope 2",
        "(VPM) data from database for telescope 3",
        "(VPM) data from database for telescope 4",
        "(VPM) data from database for telescope 5",
    ]

    vpm0 = getRunQuality("")
    vpm = [0] * len(vpm_messages)
    for tel in range(len(vpm_messages)):
        l1 = vpm_messages[tel: tel + 1]
        l1.append("blabla")
        vpm[tel] = getRunQuality(l1)
    l1 = vpm_messages[::2]
    vpm2 = getRunQuality(l1)
    vpm_all = getRunQuality(vpm_messages)

    assert (
        vpm0 == 120
        and vpm[0] == 112
        and vpm[1] == 104
        and vpm[2] == 88
        and vpm[3] == 56
        and vpm2 == 80
        and vpm_all == 0
    )


def test_getGTI():
    run_start_from_reference = 0
    # runno 300 20 0
    BitArray0 = np.array([255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 15, 0,
                          0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, ])

    # runno 0 20 0
    BitArray1 = np.array([0, 0, 240, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, ])

    # runno 580 20 0
    BitArray2 = np.array([255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 15, 0, 0, ])

    # runno 300 20 0
    # runno 0 20 0
    # runno 580 20 0
    BitArray3 = np.array([0, 0, 240, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 15, 0,
                          0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 15, 0, 0, ])

    # runno 300 9 0
    BitArray4 = np.array([255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 15, 224,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, ])

    # runno 3 5 0
    BitArray5 = np.array([7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, ])

    # runno 0 20 0
    # runno 300 9 0
    BitArray6 = np.array([0, 0, 240, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 15, 224,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, ])

    # No time mask
    BitArray7 = np.array([255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                          255, 255, 255, 255, 255, 255, 255, 255, 255, 255, ])

    start0, stop0, ontime_s0 = getGTI(BitArray0, run_start_from_reference)
    start1, stop1, ontime_s1 = getGTI(BitArray1, run_start_from_reference)
    start2, stop2, ontime_s2 = getGTI(BitArray2, run_start_from_reference)
    start3, stop3, ontime_s3 = getGTI(BitArray3, run_start_from_reference)
    start4, stop4, ontime_s4 = getGTI(BitArray4, run_start_from_reference)
    start5, stop5, ontime_s5 = getGTI(BitArray5, run_start_from_reference)
    start6, stop6, ontime_s6 = getGTI(BitArray6, run_start_from_reference)
    start7, stop7, ontime_s7 = getGTI(BitArray7, run_start_from_reference)

    assert all([a == b for a, b in zip(start0, [0, 320])])
    assert all([a == b for a, b in zip(stop0, [300, 600])])
    assert all([a == b for a, b in zip(start1, [20])])
    assert all([a == b for a, b in zip(stop1, [600])])
    assert all([a == b for a, b in zip(start2, [0])])
    assert all([a == b for a, b in zip(stop2, [580])])
    assert all([a == b for a, b in zip(start3, [20, 320])])
    assert all([a == b for a, b in zip(stop3, [300, 580])])
    assert all([a == b for a, b in zip(start4, [0, 309])])
    assert all([a == b for a, b in zip(stop4, [300, 600])])
    assert all([a == b for a, b in zip(start5, [0, 8])])
    assert all([a == b for a, b in zip(stop5, [3, 600])])
    assert all([a == b for a, b in zip(start6, [20, 309])])
    assert all([a == b for a, b in zip(stop6, [300, 600])])
    assert all([a == b for a, b in zip(start7, [0])])
    assert all([a == b for a, b in zip(stop7, [600])])

    assert (
        ontime_s0 == 580
        and ontime_s1 == 580
        and ontime_s2 == 580
        and ontime_s3 == 540
        and ontime_s4 == 591
        and ontime_s5 == 595
        and ontime_s6 == 571
        and ontime_s7 == 600
    )


if __name__ == "__main__":
    test_getRunQuality()
    test_getGTI()
