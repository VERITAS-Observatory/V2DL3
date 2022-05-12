import logging

import numpy as np


class WrongIrf(Exception):
    def __init__(self, message, errors):
        """Call the base class constructor with the parameters it needs"""
        super().__init__(message)

        self.errors = errors


def produce_tel_list(tel_config):
    """Convert the list of telescopes into a string for FITS header"""
    tel_list = "".join("T" + str(tel) + "," for tel in tel_config["TelType"])
    return tel_list[:-1]


def bin_centers_to_edges(axis, logaxis=True):
    """Calculate bin centers of two axes"""
    bin_size = axis[1] - axis[0]
    extended_axis = np.insert(axis, 0, axis[0] - bin_size)
    bin_edges = extended_axis + bin_size / 2.0
    if logaxis:
        bin_low = np.power(10, [bin_edges[:-1]])[0]
        bin_high = np.power(10, [bin_edges[1:]])[0]
    else:
        bin_low = np.array([bin_edges[:-1]])[0]
        bin_high = np.array([bin_edges[1:]])[0]
    return bin_edges, bin_low, bin_high


def duplicate_dimension(data, axis):
    """Duplicate a single axis, assuming it's length is 1."""
    current_shape = np.shape(data)
    corrected_shape = [2 if i == axis else k for i, k in enumerate(current_shape)]
    logging.info(current_shape, corrected_shape)
    tiles = [2 if k == 1 else 1 for k in current_shape]
    return np.tile(data, tiles)


def duplicate_dimensions(data):
    """Duplicate all axes, one by one, when their size is 1."""
    new_data = data
    current_shape = np.shape(new_data)
    for i, dim in enumerate(current_shape):
        if dim == 1:
            new_data = duplicate_dimension(new_data, i)
    return new_data


def getGTI(BitArray, run_start_from_reference):
    """Decode the time masks stored as 'TBits' in anasum file and extract GTIs

    Parameters
    ----------
    maskBits :  array of uint8 numbers, read from anasum root file
    run_start_from_reference: Start time of the run in second
                              from reference time

    Retuns
    ------
    gti_start_from_reference : numpy array of start time of GTIs in second
                               from reference time
    gti_end_from_reference: numpy array of stop time of GTIs in second from
                            reference time
    ontime: Total of good times

    """

    n = BitArray.size
    TimeArray_b = ""
    for i in range(n):
        TimeArray_b = TimeArray_b + np.binary_repr(BitArray[i], width=8)[::-1]

    nbits = len(TimeArray_b)
    for i in range(nbits):
        if (TimeArray_b[-1] == "0"):
            TimeArray_b = TimeArray_b[:-1]
        else:
            break

    duration_s = len(TimeArray_b)
    ontime_s = TimeArray_b.count("1")
    logging.info(
        "Duration: {0:.0f} (sec.) {1:.2f} (min)".format(duration_s, duration_s / 60.0)
    )
    logging.info(
        "Ontime: {0:.0f} (sec.) {1:.2f} (min)".format(ontime_s, ontime_s / 60.0)
    )

    gti_start = []
    gti_end = []

    if TimeArray_b[0] != "0":
        gti_start.append(0)

    for i in range(1, duration_s - 1):

        if (TimeArray_b[i] == "0") and (TimeArray_b[i - 1] == "1"):
            end = i
            gti_end.append(end)

        if (TimeArray_b[i] == "0") and (TimeArray_b[i + 1] == "1"):
            start = i + 1
            gti_start.append(start)

    if (TimeArray_b[-1] != 0):
        gti_end.append(duration_s)

    logging.info(
        "GTIs start and stop in second since run start: {0} {1}".format(
            gti_start, gti_end
        )
    )

    gti_start_from_reference = np.array(gti_start) + run_start_from_reference
    gti_end_from_reference = np.array(gti_end) + run_start_from_reference

    return gti_start_from_reference, gti_end_from_reference, ontime_s


def getRunQuality(logdata, ntel=4):
    """Evaluate the run quality based on VPM data used or not in the evndisp.

    L to R: bit0 not used, bit[1-4] set when VPM data not used for
            corresponding telescope, bit[5-7] reserved for QUALITY
            flag defined in GADF and not used currently.

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

    vpm = 2 ** (ntel + 3) - 2 ** 3
    if np.size(logdata) <= 1:
        return vpm

    for line in range(np.size(logdata)):
        for tel in range(ntel):
            vpm_str = "(VPM) data from database for telescope " + str(tel + 1)
            if vpm_str in logdata[line]:
                vpm &= ~(1 << (tel + 3))

    logging.info(f"Run quality flag: {vpm} (8 Bit code: {vpm:b})")

    return vpm
