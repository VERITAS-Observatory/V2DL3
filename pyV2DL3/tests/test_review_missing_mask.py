import pytest
from pyV2DL3.eventdisplay.fillEVENTS import __get_ontime

@pytest.mark.parametrize("stereo", [{}, {"timeMask": {}}])
def test_missing_directory_and_missing_object(stereo):
    start, stop, ontime = __get_ontime({"run_1": {"stereo": stereo}}, 1, 100., 120.5)
    assert start == [100.]
    assert stop == [120.5]
    assert ontime == 20.5
