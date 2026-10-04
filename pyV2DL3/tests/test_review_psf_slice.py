import numpy as np
from pyV2DL3.eventdisplay.fillRESPONSE import fill_direction_migration

class Interpolator:
    def set_irf(self, *args, **kwargs): pass
    def interpolate(self, coords):
        return np.ones((3, 3)), [np.array([-1., 0., 1.]), np.array([-2., -1., 0.])]

def test_psf_axis_ending_below_upper_cut_is_not_empty():
    result = fill_direction_migration(Interpolator(), [0.5], 1., 20., [0.5], [0.5])
    assert result["ENERG_LO"].shape == (1, 3)
    assert result["RPSF"].shape[-1] == 3
