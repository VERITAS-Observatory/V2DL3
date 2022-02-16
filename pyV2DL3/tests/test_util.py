from pyV2DL3.eventdisplay.util import getRunQuality


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
        l1 = vpm_messages[tel : tel + 1]
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


if __name__ == "__main__":
    test_getRunQuality()
