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
    """Return the closest paired azimuth interval using circular distance.

    Array positions remain associated with their original bin IDs. The
    all-azimuth sentinel is only used if there are no directional bins.
    """
    if azimuth is None or not np.isfinite(azimuth):
        raise ValueError("A finite azimuth is required")
    azMins, azMaxs = np.asarray(azMins), np.asarray(azMaxs)
    if azMins.shape != azMaxs.shape or azMins.size == 0:
        raise ValueError("Azimuth bounds must be paired and nonempty")
    valid = (np.abs(azMins) <= 180) & (np.abs(azMaxs) <= 180)
    indices = np.flatnonzero(valid)
    if not len(indices):
        if len(azMins) == 1 and azMins[0] <= -1000 and azMaxs[0] >= 1000:
            return 0
        raise ValueError("No directional azimuth bins found")
    centres = (azMins[valid] + (azMaxs[valid] - azMins[valid]) % 360 / 2) % 360
    distance = np.abs((centres - azimuth + 180) % 360 - 180)
    return indices[np.argmin(distance)]


def get_empty_ndarray(data_dimension):
    """return a zero filled ndarray with the given dimensions"""
    return np.zeros(tuple(data_dimension))


def _get_az_mask(azimuth, fast_eff_area):
    """Select an azimuth ID while preserving its paired bounds."""
    ids = fast_eff_area["az"].array(library="np")
    records = np.unique(np.column_stack((
        ids,
        fast_eff_area["azMin"].array(library="np"),
        fast_eff_area["azMax"].array(library="np"),
    )), axis=0)
    if len(np.unique(records[:, 0])) != len(records):
        raise ValueError("Inconsistent bounds for an IRF azimuth ID")
    index = find_closest_az(azimuth, records[:, 1], records[:, 2])
    return ids == records[index, 0]


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

    if azimuth is None:
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
