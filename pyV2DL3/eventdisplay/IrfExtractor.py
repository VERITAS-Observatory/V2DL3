import logging
import sys

import awkward as ak
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

    Note the different conventions for azimuth:
    - anasum file (0..360)
    - EA (-180..180)

    """
    az_centers = np.array((azMaxs[:-1] + azMins[1:]) / 2.0)
    az_centers[az_centers < 0] += 360
    if np.any(az_centers < 0) or np.any(az_centers > 360):
        logging.error("IRF azimuth bins not in the range 0-360")
        raise ValueError
    return find_nearest(az_centers, azimuth)


def get_empty_ndarray(data_dimension):
    """return a zero filled ndarray with the given dimensions"""
    return np.zeros(tuple(data_dimension))


def _get_az_mask(azimuth, fast_eff_area):
    """Return azimuth mask for given azimuth angle"""
    _, azMaxs = load_parameter("azMax", fast_eff_area)
    _, azMins = load_parameter("azMin", fast_eff_area)
    az_bin_to_store = find_closest_az(azimuth, azMins, azMaxs)
    return fast_eff_area["az"].array(library="np") == az_bin_to_store


def extract_irf_1d(filename, irf_name, azimuth=None):
    """
    Extract 1D IRF from effective area file

    return a multidimensional array
    """

    fast_eff_area = uproot.open(filename)["fEffAreaH2F"]
    az_mask = _get_az_mask(azimuth, fast_eff_area)
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
            logging.error(f"At entry number {i} unexpected error: {sys.exc_info()[0]}")
            raise

    axes = {
        'energies': energies[0],
        'pedvars': pedvars,
        'zeniths': zds,
        'woffs': woffs
    }

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
    """
    Extract 2D IRF from effective area file.

    Returns a multidimensional array with axes:

    - irf_dimension_1
    - irf_dimension_2
    - pedvars
    - zeniths
    - woffs

    For azimuth, select the bin closest to the given azimuth angle.

    Parameters
    ----------
    filename : str
        Path to the effective area file.
    irf_name : str
        Name of the IRF to extract (e.g., 'eff' or 'hEsysMCRelative2D').
    azimuth : float, optional
        Azimuth angle in degrees.

    """

    fast_eff_area = uproot.open(filename)["fEffAreaH2F"]
    az_mask = _get_az_mask(azimuth, fast_eff_area)

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
            logging.error("Unexpected error:", sys.exc_info()[0])
            logging.error("Entry number ", i)
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
        logging.error("Azimuth for IRF extraction not given")
        raise ValueError

    irf_fn = extract_irf_1d if irf1d else extract_irf_2d
    return irf_fn(filename, irf_name, azimuth)


def extract_irf_for_knn(filename, irf_name, irf1d=False, azimuth=None):
    """Extract IRF for KNeighborsRegressor"""
    fast_eff_area = uproot.open(filename)["fEffAreaH2F"]
    az_mask = _get_az_mask(azimuth, fast_eff_area)

    ze = 1. / np.cos(np.radians(fast_eff_area["ze"].array()[az_mask]))
    pedvar = fast_eff_area["pedvar"].array()[az_mask]
    woff = fast_eff_area["Woff"].array()[az_mask]

    if irf1d:
        e0 = fast_eff_area["e0"].array()[az_mask]
        ze_b, pedvar_b, woff_b = ak.broadcast_arrays(e0, ze, pedvar, woff)[1:]
        coords = np.vstack([
            ak.to_numpy(ak.flatten(pedvar_b)),
            ak.to_numpy(ak.flatten(ze_b)),
            ak.to_numpy(ak.flatten(woff_b)),
            ak.to_numpy(ak.flatten(e0)),
        ]).T
        values = ak.to_numpy(ak.flatten(fast_eff_area[irf_name].array()[az_mask]))
    else:
        irf_axis_x = read_irf_axis("x", fast_eff_area, irf_name, az_mask)
        irf_axis_y = read_irf_axis("y", fast_eff_area, irf_name, az_mask)
        values = ak.to_numpy(ak.flatten(fast_eff_area[irf_name + "_value"].array()[az_mask]))

        ze_rep = np.repeat(ak.to_numpy(ze), len(irf_axis_x) * len(irf_axis_y)).astype(np.float32)
        pedvar_rep = np.repeat(ak.to_numpy(pedvar), len(irf_axis_x) * len(irf_axis_y)).astype(np.float32)
        woff_rep = np.repeat(ak.to_numpy(woff), len(irf_axis_x) * len(irf_axis_y)).astype(np.float32)

        xx, yy = np.meshgrid(irf_axis_x, irf_axis_y, indexing='xy')
        irf_dim1 = np.tile(xx.flatten(), len(ze))
        irf_dim2 = np.tile(yy.flatten(), len(ze))

        coords = np.vstack([
            pedvar_rep.flatten(),
            ze_rep.flatten(),
            woff_rep.flatten(),
            irf_dim1.flatten(),
            irf_dim2.flatten(),
        ]).T

    return coords, values
