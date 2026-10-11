import numpy as np
import pytest

from pyV2DL3.eventdisplay.fillEVENTS import __fill_event_list
from pyV2DL3.eventdisplay.util import ZeroLengthEventList


class Tree:
    def arrays(self, library):
        return {"eventNumber": np.array([1]), "RA": np.array([0.]), "IsGamma": np.array([0])}


def test_nonempty_tree_with_no_selected_events(caplog):
    with pytest.raises(ZeroLengthEventList):
        __fill_event_list(
            {"run_1/stereo/DL3EventTree": Tree()}, 1, {"IsGamma": 1}, 0
        )

    assert "Empty event list after applying selection filter" in caplog.text
