import numpy as np
from ROOT import TFile
import sys
from tqdm.auto import tqdm
import uproot


class WrongIrf(Exception):
    def __init__(self, message, errors):
        # Call the base class constructor with the parameters it needs
        super().__init__(message)

        # Now for your custom code...
        self.errors = errors


def produce_tel_list(tel_config):
    """Convert the list of telescopes into a string for FITS header"""
    tel_list = ""
    for tel in tel_config["TelType"]:
        tel_list += "T" + str(tel) + ","
    return tel_list[:-1]


def bin_centers_to_edges(axis):
    """This function assumes bins of equal width"""
    bin_size = axis[1] - axis[0]
    extended_axis = np.insert(axis, 0, axis[0] - bin_size)
    return extended_axis + bin_size / 2.0


def duplicate_dimensions(data):
    """This function duplicates all axes, one by one, when their size is 1."""
    new_data = data
    current_shape = np.shape(new_data)
    for i, dim in enumerate(current_shape):
        if dim == 1:
            new_data = duplicate_dimension(new_data, i)
    return new_data


def getGTI(BitArray, run_start_from_reference):
    """
     Function to decode the time masks which is stored as 'TBits' in anasum ROOT file
     and extract the GTIs for a given run

     Parameters
     ----------
     maskBits :  array of uint8 numbers, read from anasum root file
     run_start_from_reference: Start time of the run in second from reference time

     Retuns
     ------
     gti_start_from_reference : numpy array of start time of GTIs in second from reference time
     gti_end_from_reference: numpy array of stop time of GTIs in second from reference time
     ontime: Total of good times

    """

    n = BitArray.size
    TimeArray_s = []
    for i in range(n):
        TimeArray_s.append(np.binary_repr(BitArray[i]).count('1'))

    duration_s = (n - 1) * 8 + TimeArray_s[-1]
    ontime_s = np.sum(TimeArray_s[0:n])
    print('Duration:', duration_s, '(sec.)', duration_s / 60., '(min.)')
    print('Ontime', ontime_s, '(sec.)', ontime_s / 60., '(min.)')

    gti_start = []
    gti_end = []

    if (TimeArray_s[0] != 0):
        gti_start.append(0)

    for i in range(1, n - 1, 1):
        if ((TimeArray_s[i] == 0) & (TimeArray_s[i - 1] != 0)):
            end = (i * 8 - (8 - TimeArray_s[i - 1]))
            gti_end.append(end)

        if ((TimeArray_s[i] == 0) & (TimeArray_s[i + 1] != 0)):
            start = ((i + 1) * 8 + (8 - TimeArray_s[i + 1]))
            gti_start.append(start)

    if (TimeArray_s[-1] != 0):
        gti_end.append(duration_s)

    print('GTIs start and stop in second since run start:', gti_start, gti_end)

    gti_start_from_reference = np.zeros(np.size(gti_start))
    gti_end_from_reference = np.zeros(np.size(gti_end))

    for i in range(np.size(gti_start)):
        gti_start_from_reference[i] = gti_start[i] + run_start_from_reference
        gti_end_from_reference[i] = gti_end[i] + run_start_from_reference

    return gti_start_from_reference, gti_end_from_reference, ontime_s

def getRunQuality(logdata):
    """
    Function to evaluate the run quality based on VPM data used or not in the evndisp.
    L to R: bit0 not used, bit[1-4] set when VPM data not used for corresponding telescope,
            bit[5-7] reserved for QUALITY flag defined in GADF and not used currently.

    eg: flag: 64  (01000000) means for T1 VPM data are not used
        flag: 120 (01111000) means for T1, T2, T3, T4 VPM data are not used
        flag: 0   (00000000) means for all four telescopes VPM data are used

    Parameters
    ----------
    TMacro read from the anasum ROOT file

    Returns
    -------
    A 8-bit coded integer.

    """

    vpm1 = "(VPM) data from database for telescope 1"
    vpm2 = "(VPM) data from database for telescope 2"
    vpm3 = "(VPM) data from database for telescope 3"
    vpm4 = "(VPM) data from database for telescope 4"

    res1, res2, res3, res4 = "1", "1", "1", "1"

    for line in range(np.size(logdata)):

        if vpm1 in logdata[line]:
            res1 = "0"
        elif vpm2 in logdata[line]:
            res2 = "0"
        elif vpm3 in logdata[line]:
            res3 = "0"
        elif vpm4 in logdata[line]:
            res4 = "0"
        else:
            status = "No VPM data used for any of the telescopes!"

    vpm_used = "0" + res1 + res2 + res3 + res4 + "000"
    flag = int(vpm_used, 2)
    print("Run quality flag: {} (8 Bit code: {})".format(flag, vpm_used))

    return flag
