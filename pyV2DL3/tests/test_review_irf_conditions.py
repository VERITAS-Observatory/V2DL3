import numpy as np
from pyV2DL3.eventdisplay.fillEVENTS import __fill_event_list, __get_average_event_direction

class Tree:
    def arrays(self, library):
        columns = ["eventNumber", "timeOfDay", "RA", "DEC", "El", "Az", "Energy", "NImages", "Xoff", "Yoff", "ImgSel"]
        tree = {c: np.array([1., 2.]) for c in columns}
        tree.update(MJD=np.array([60000, 60000]), IsGamma=np.array([1, 0]), MeanPedvar=np.array([1., 9.]))
        return tree

def test_pedvar_independent_of_output_filter():
    file = {"run_1/stereo/DL3EventTree": Tree()}
    all_events = __fill_event_list(file, 1, None, 0)
    selected = __fill_event_list(file, 1, {"IsGamma": 1}, 0)
    assert all_events[2] == selected[2] == 5.

def test_pointing_mean_wraps_at_north():
    altitude, azimuth = __get_average_event_direction(np.array([60., 60.]), np.array([359., 1.]))
    assert altitude == 60.
    assert min(abs(azimuth), abs(azimuth-360)) < 1e-10
