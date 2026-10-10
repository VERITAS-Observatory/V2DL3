import numpy as np
import pytest
from pyV2DL3.eventdisplay.util import getGTI

@pytest.mark.parametrize("byte", range(256))
def test_all_single_byte_masks(byte):
    start, stop, ontime = getGTI(np.array([byte]), 100)
    decoded = np.zeros(8, dtype=bool)
    assert len(start) == len(stop)
    for a, b in zip(start, stop):
        decoded[int(a-100):int(b-100)] = True
    assert decoded.tolist() == [bool(byte & (1 << i)) for i in range(8)]
    assert ontime == decoded.sum()

def test_partial_second_and_padding():
    start, stop, ontime = getGTI(np.array([255]), 100, nbits=3, run_duration=2.25)
    assert start.tolist() == [100]
    assert stop.tolist() == [102.25]
    assert ontime == 2.25
