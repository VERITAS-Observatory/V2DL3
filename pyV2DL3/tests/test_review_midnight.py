import numpy as np
from astropy.time import Time
from pyV2DL3.constant import VTS_REFERENCE_MJD
from pyV2DL3.eventdisplay.fillEVENTS import __get_time_vector

def test_events_use_their_own_utc_date():
    days = np.array([60000, 60001, 60000])
    tod = np.array([86399., 1., 123.])
    result = __get_time_vector(tod, -1, event_mjd=days)
    ref = Time(VTS_REFERENCE_MJD, format="mjd", scale="utc")
    expected = (Time(days, format="mjd", scale="utc") - ref).sec + tod
    np.testing.assert_allclose(result, expected, rtol=0, atol=1e-6)
    assert abs(result[1] - result[0] - 2) < 1e-6
