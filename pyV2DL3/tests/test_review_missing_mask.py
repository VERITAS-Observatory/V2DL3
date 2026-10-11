import logging

import pytest

from pyV2DL3.eventdisplay.fillEVENTS import __get_ontime


@pytest.mark.parametrize("stereo", [{}, {"timeMask": {}}])
def test_missing_directory_and_missing_object_log_fallback(stereo, caplog):
    caplog.set_level(logging.INFO, logger="pyV2DL3.eventdisplay.fillEVENTS")

    start, stop, ontime = __get_ontime({"run_1": {"stereo": stereo}}, 1, 100., 120.5)
    assert start == [100.]
    assert stop == [120.5]
    assert ontime == 20.5

    assert "run 1" in caplog.text
    assert "using the full run interval" in caplog.text
    if "timeMask" not in stereo:
        assert "time mask absent" in caplog.text
    else:
        assert "time mask has no maskBits" in caplog.text
