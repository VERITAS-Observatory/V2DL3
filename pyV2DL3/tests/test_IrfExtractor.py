import numpy as np
from pytest import approx
from pyV2DL3.eventdisplay.IrfExtractor import find_closest_az
from pyV2DL3.eventdisplay.IrfExtractor import find_nearest

def test_find_nearest():
   az_centers = np.array([
                187.5, 206.25, 232.5, 255., 277.5, 300., 322.5, 345., 
                7.5, 30., 52.5, 75., 97.5, 120., 142.5, 161.25 ])
   az_50 = find_nearest(az_centers, 50.)
   az_188 = find_nearest(az_centers, 188.)

   assert az_50 == 10 and az_188 == 0

def test_find_closest_az():
   azMins = np.array([
                -1000., -180., -157.5, -135.,
                -112.5, -90., -67.5, -45., -22.5,
                0., 22.5, 45., 67.5, 90., 112.5, 135., 150.])
   azMaxs = np.array([
                -165., -150., -120., -97.5, -75., -52.5, -30., 
                -7.5, 15., 37.5, 60., 82.5, 105., 127.5, 150., 
                172.5, 1000. ])
   az_bin_to_store_1 = find_closest_az(146.20, azMins, azMaxs)
   az_bin_to_store_2 = find_closest_az(-180., azMins, azMaxs)
   assert az_bin_to_store_1 == 14 and \
          az_bin_to_store_2 == 8

if __name__ == '__main__':
   test_find_nearest()
   test_find_closest_index()
