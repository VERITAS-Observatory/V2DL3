import logging
import sys

import numpy as np
import uproot


def remove_duplicities(array, atol):
    """remove duplicates in array allowing for a given precision"""
    i = 0
    while i < len(array) - 1:
        i += 1
        if np.isclose(array[i - 1], array[i], atol=atol):
            array = np.delete(array, i - 1)
            i -= 1
    return array


def find_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()


def load_parameter(parameter_name, fast_eff_area, az_mask=None):
    """load effective area parameter

    apply necessary rounding
    """

    if az_mask is None:
        all_par = fast_eff_area[parameter_name].array(library="np")
    else:
        all_par = fast_eff_area[parameter_name].array(library="np")[az_mask]
    par = []
    # round all parameters for correct extraction
    all_par = np.round(all_par, decimals=2)
    par = np.unique(all_par)
    if parameter_name == "pedvar":
        par = remove_duplicities(par, 0.21)
        # replace all_pedvars values with nearest from pedvars list
        # by calculating the difference between each element and taking the min
        all_par = par[abs(all_par[None, :] - par[:, None]).argmin(axis=0)]
    elif parameter_name == "ze":
        par = remove_duplicities(par, 2.0)

    return all_par, par


def find_closest_az(azimuth, azMins, azMaxs):
    """find closest azimuth bin

    note the different conventions for azimuth:
    - anasum file (0..360)
    - EA (-180..180)

    """
    az_centers = np.array((azMaxs[:-1] + azMins[1:]) / 2.0)
    az_centers[az_centers < 0] += 360
    if np.any(az_centers < 0) or np.any(az_centers > 360):
        raise ValueError("IRF azimuth bins not in the range 0-360")
    return find_nearest(az_centers, azimuth)


def get_empty_ndarray(data_dimension):
    """return a zero filled ndarray with the given dimensions"""
    return np.zeros(tuple(data_dimension))


def extract_irf_1d(filename, irf_name, azimuth=None):
    """extract 1D IRF from effective area file

    return a multidimensional array
    """

    fast_eff_area = uproot.open(filename)["fEffAreaH2F"]

    # select az bin and define az mask
    _, azMaxs = load_parameter("azMax", fast_eff_area)
    _, azMins = load_parameter("azMin", fast_eff_area)
    az_bin_to_store = find_closest_az(azimuth, azMins, azMaxs)
    az_mask = fast_eff_area["az"].array(library="np") == az_bin_to_store

    energies = fast_eff_area["e0"].array(library="np")[az_mask]
    irf = fast_eff_area[irf_name].array(library="np")[az_mask]

    all_pedvars, pedvars = load_parameter("pedvar", fast_eff_area, az_mask)
    all_zds, zds = load_parameter("ze", fast_eff_area, az_mask)
    all_Woffs, woffs = load_parameter("Woff", fast_eff_area, az_mask)

    data = get_empty_ndarray([len(irf[0]), len(pedvars), len(zds), len(woffs)])

    for i in range(len(irf)):
        try:
            data[
                :,
                find_nearest(pedvars, all_pedvars[i]),
                find_nearest(zds, all_zds[i]),
                find_nearest(woffs, all_Woffs[i]),
            ] = irf[i]
        except Exception:
            logging.exception("Unexpected error:", sys.exc_info()[0])
            logging.exception("Entry number ", i)
            raise

    axes = {'energies': energies[0], 'pedvars': pedvars, 'zeniths': zds, 'woffs': woffs}

    return data, axes


def read_irf_axis(xy, fast_eff_area, irf_name, az_mask):
    """return irf axis (bin centres)"""

    nbins = fast_eff_area[irf_name + "_bins" + xy].array(library="np")[az_mask]
    c_min = fast_eff_area[irf_name + "_min" + xy].array(library="np")[az_mask]
    c_max = fast_eff_area[irf_name + "_max" + xy].array(library="np")[az_mask]

    if nbins[0] > 0:
        binwidth = (c_max[0] - c_min[0]) / nbins[0] / 2.0
        return np.linspace(c_min[0] + binwidth, c_max[0] - binwidth, nbins[0])

    return None


def extract_irf_2d(filename, irf_name, azimuth=None):
    """extract 2D IRF from effective area file

    return a multidimensional array
    """

    fast_eff_area = uproot.open(filename)["fEffAreaH2F"]

    # select az bin and define az mask
    _, azMaxs = load_parameter("azMax", fast_eff_area)
    _, azMins = load_parameter("azMin", fast_eff_area)
    az_bin_to_store = find_closest_az(azimuth, azMins, azMaxs)
    az_mask = fast_eff_area["az"].array(library="np") == az_bin_to_store

    # IRF axes and values
    irf_dimension_1 = read_irf_axis("x", fast_eff_area, irf_name, az_mask)
    irf_dimension_2 = read_irf_axis("y", fast_eff_area, irf_name, az_mask)
    irf2D = fast_eff_area[irf_name + "_value"].array(library="np")[az_mask]

    # parameter space
    all_pedvars, pedvars = load_parameter("pedvar", fast_eff_area, az_mask)
    all_zds, zds = load_parameter("ze", fast_eff_area, az_mask)
    all_Woffs, woffs = load_parameter("Woff", fast_eff_area, az_mask)

    data = get_empty_ndarray(
        [len(irf_dimension_1), len(irf_dimension_2), len(pedvars), len(zds), len(woffs)]
    )

    for i in range(len(irf2D)):
        irf = np.reshape(irf2D[i], (-1, len(irf_dimension_2)), order="F")
        try:
            data[
                :,
                :,
                find_nearest(pedvars, all_pedvars[i]),
                find_nearest(zds, all_zds[i]),
                find_nearest(woffs, all_Woffs[i]),
            ] = irf
        except Exception:
            logging.exception("Unexpected error:", sys.exc_info()[0])
            logging.exception("Entry number ", i)
            raise

    axes = {
        'irf_dimension_1': np.array(irf_dimension_1),
        'irf_dimension_2': np.array(irf_dimension_2),
        'pedvars': pedvars,
        'zeniths': zds,
        'woffs': woffs
    }

    return data, axes


def extract_irf(filename, irf_name, azimuth=None, irf1d=False):
    """extract IRF from effective area file

    return a multidimensional array
    """

    if not azimuth:
        raise ValueError("Azimuth for IRF extraction not given")

    irf_fn = extract_irf_1d if irf1d else extract_irf_2d
    return irf_fn(filename, irf_name, azimuth)
