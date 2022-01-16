import numpy as np
from ROOT import TFile
import sys
from tqdm.auto import tqdm
import uproot


def produce_tel_list(tel_config):
    # Convert the list of telescopes into a string for FITS header
    tel_list = ""
    for tel in tel_config["TelType"]:
        tel_list += "T" + str(tel) + ","
    return tel_list[:-1]


def is_close(a, b, rel_tol=1e-04, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class WrongIrf(Exception):
    def __init__(self, message, errors):
        # Call the base class constructor with the parameters it needs
        super().__init__(message)

        # Now for your custom code...
        self.errors = errors


def remove_duplicities(array, atol):
    i = 0
    while i < len(array) - 1:
        i += 1
        if np.isclose(array[i - 1], array[i], atol=atol):
            array = np.delete(array, i - 1)
            i -= 1
    return array


# This seems not to be working on ROOT 5 + python 6.
def graph_to_array_y(graph):
    y = [g for g in graph.GetY()]
    return y


# This seems not to be working on ROOT 5 + python 6.
def graph_to_array_x(graph):
    x = [g for g in graph.GetX()]
    return x


def graph_to_array(graph, nbins):
    x = [graph.GetX()[i] for i in np.arange(0, nbins)]
    y = [graph.GetY()[i] for i in np.arange(0, nbins)]
    return x, y


def bin_edges_to_centers(axis):
    # This function assumes bins of equal width
    bin_size = axis[1] - axis[0]
    return np.delete(axis + bin_size / 2.0, len(axis) - 1)


def bin_centers_to_edges(axis):
    # This function assumes bins of equal width
    bin_size = axis[1] - axis[0]
    extended_axis = np.insert(axis, 0, axis[0] - bin_size)
    return extended_axis + bin_size / 2.0


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def hist2array(h, return_edges=False):
    # extracted TH2D to numpy conversion part from root_numpy
    cName = h.Class_Name()
    if cName == "TH2D" or cName == "TH2F":
        nBinsX = h.GetNbinsX()
        nBinsY = h.GetNbinsY()
        shape = (nBinsY + 2, nBinsX + 2)
        if cName == "TH2F":
            dtype = "f4"
        else:
            dtype = "f8"
        array = np.ndarray(shape=shape, dtype=dtype, buffer=h.GetArray())
    else:
        print(cName)
    array = array[tuple([slice(1, -1) for idim in range(array.ndim)])]
    if return_edges:
        ndims = h.GetDimension()
        axis_getters = ["GetXaxis", "GetYaxis", "GetZaxis"][:ndims]
        edges = []
        for idim, axis_getter in zip(range(ndims), axis_getters):
            ax = getattr(h, axis_getter)(*(()))
            edges.append(np.empty(ax.GetNbins() + 1, dtype=np.double))
            ax.GetLowEdge(edges[-1])
            edges[-1][-1] = ax.GetBinUpEdge(ax.GetNbins())
        return array.T, edges
    else:
        return array.T

def load_parameter(
        parameter_name,
        fast_eff_area ):
    """load effective area parameter

       apply necessary rounding
    """

    all_par = fast_eff_area[parameter_name].array(library="np")
    par = []
    # round all parameters for correct extraction
    all_par = np.round(all_par, decimals=2)
    par = np.unique(all_par)
    if parameter_name is 'pedvar':
        par = remove_duplicities(par, 0.21)
        # replace all_pedvars values with nearest from pedvars list
        # by calculating the difference between each element and taking the min
        all_par = par[
        abs(all_par[None, :] - par[:, None]).argmin(axis=0)
    ]
    elif parameter_name is 'ze':
        par = remove_duplicities(par, 2.0)

    return all_par, par

def read_1d_samplearrays(
    irf_name,
    fast_eff_area,
    all_nbins ):

    entry_with_max_bins = find_nearest(all_nbins, all_nbins.max())
    sample_energies = fast_eff_area['e0'].array(library="np",
                                    entry_start=entry_with_max_bins,
                                    entry_stop=entry_with_max_bins+1)[0]
    sample_irf = fast_eff_area[irf_name].array(library="np",
                                    entry_start=entry_with_max_bins,
                                    entry_stop=entry_with_max_bins+1)[0]
    return sample_energies, sample_irf, entry_with_max_bins

def find_closest_index( indexs,
                        irf_name ):
    """find closest index stored"""

    # true energy IRFs are not index dependent, only lowest index is populated
    if (
        irf_name == "hAngularLogDiffEmc_2D"
        or irf_name == "hEsysMCRelative2DNoDirectionCut"
        or irf_name == "hEsysMCRelative2D"
    ):
        return indexs.min()

    return indexs[find_nearest(indexs, (indexs.min() + indexs.max()) / 2.0)]


def find_closest_az( azimuth,
                     azMins,
                     azMaxs ):
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


def extract_irf(
    filename,
    irf_name,
    azimuth=None):
    """
    extract IRF from effective area file and return a
    multidimensional array 

    """

    print("Extracting IRFs of type: {}".format(irf_name))
    print("    Azimuth", azimuth)
    if not azimuth:
        raise ValueError("Azimuth for IRF extraction not given")

    # List of implemented IRFs
    implemented_irf_names_1d = ["eff", "effNoTh2", "Rec_eff"]
    implemented_irf_names_2d = [
        "hEsysMCRelative2D",
        "hEsysMCRelative2DNoDirectionCut",
        "hAngularLogDiffEmc_2D",
    ]
    # Get both the ROOT effective area TTree and the uproot one (much faster)
    eff_area_file = TFile.Open(filename)
    eff_area_tree = eff_area_file.Get("fEffArea")
    fast_eff_area = uproot.open(filename)["fEffArea"]

    all_zds, zds = load_parameter( "ze", fast_eff_area)
    all_azs, azs = load_parameter( "az", fast_eff_area)
    all_azMins, azMins = load_parameter( "azMin", fast_eff_area)
    all_azMaxs, azMaxs = load_parameter( "azMax", fast_eff_area)
    all_Woffs, woffs = load_parameter( "Woff", fast_eff_area)
    all_pedvars, pedvars = load_parameter( "pedvar", fast_eff_area)
    all_indexs, indexs = load_parameter( "index", fast_eff_area)
    all_nbins, nbins = load_parameter( "nbins", fast_eff_area)

    if len(all_zds) != len(zds) * len(azs) * len(woffs) * len(pedvars) * len(
        indexs
    ):
        raise ValueError(
            "Wrong dimensions extracted from IRF cube."
            + "Probably due to the rounding applied to the IRF coordinates.")

    # For performance, deactivate all branches except the ones needed:
    # Also get the entry with max bins to define the binning in energy
    eff_area_tree.SetBranchStatus("*", 0)
    if irf_name == "eff":
        eff_area_tree.SetBranchStatus("e0", 1)
        eff_area_tree.SetBranchStatus("eff", 1)
    elif irf_name == "effNoTh2":
        eff_area_tree.SetBranchStatus("e0", 1)
        eff_area_tree.SetBranchStatus("effNoTh2", 1)
    elif irf_name == "hEsysMCRelative2D":
        eff_area_tree.SetBranchStatus("hEsysMCRelative2D", 1)
        entry_with_max_bins = 0
    elif irf_name == "hEsysMCRelative2DNoDirectionCut":
        eff_area_tree.SetBranchStatus("hEsysMCRelative2DNoDirectionCut", 1)
        entry_with_max_bins = 0
    elif irf_name == "hAngularLogDiffEmc_2D":
        eff_area_tree.SetBranchStatus("hAngularLogDiffEmc_2D", 1)
        entry_with_max_bins = 0
    else:
        raise Exception("WrongIrfName")

    if irf_name in implemented_irf_names_1d:
        sample_energies, sample_irf, entry_with_max_bins = read_1d_samplearrays(
                                              irf_name,
                                              fast_eff_area,
                                              all_nbins ) 

    # Now we know which entry we need to get in order to have a sample IRF
    #     sample_irf = sample_energies = []
    # Generate sample IRF
    for i, entry in enumerate(eff_area_tree):
        if i == entry_with_max_bins:
            if irf_name == "hEsysMCRelative2D":
                # Migration vs energy bias and true energy
                sample_irf, axes = hist2array(
                    entry.hEsysMCRelative2D, return_edges=True
                )
                # Bin edges (one more entry than migra!) for the true energy and
                # energy bias (Erec/Etrue)
                irf_dimension_1 = bin_edges_to_centers(axes[0])
                irf_dimension_2 = bin_edges_to_centers(axes[1])
            elif irf_name == "hEsysMCRelative2DNoDirectionCut":
                # Migration vs energy bias and true energy, without direction cut
                sample_irf, axes = hist2array(
                    entry.hEsysMCRelative2DNoDirectionCut, return_edges=True
                )
                # Bin edges (one more entry than migra!) for the true energy and
                # energy bias (Erec/Etrue)
                irf_dimension_1 = bin_edges_to_centers(axes[0])
                irf_dimension_2 = bin_edges_to_centers(axes[1])
            elif irf_name == "hAngularLogDiffEmc_2D":
                # PSF vs true energy:
                sample_irf, axes = hist2array(
                    entry.hAngularLogDiffEmc_2D, return_edges=True
                )
                # Bin edges (one more entry than migra!) for the true energy and
                # energy bias (Erec/Etrue)
                irf_dimension_1 = bin_edges_to_centers(axes[0])
                irf_dimension_2 = bin_edges_to_centers(axes[1])
        if i > entry_with_max_bins:
            # Only extract the above defined IRF-name-specific "entry with max bins"
            break
    # Now we should know all the dimensions that need to be stored in the output irf_data
    # * For a given azimuth sore the IRF for the closest value (remove that dimension)
    az_bin_to_store = find_closest_az(azimuth, azMins, azMaxs)
    index_to_store = find_closest_index(indexs, irf_name)
# (GM) why is this not always the lowest value (1.6)?
# (GM) for effNoTh2, a value of 3.4 is used
    print('index to store:', index_to_store, indexs.min())
    # Create data container, filled with zeros, containing the required dimensions to store
    # the IRF for a given coord_tuple. Separated between 1 and 2 dimensions:
    data_shape = []
    if irf_name in implemented_irf_names_2d:
        data_shape.append(len(irf_dimension_1))
        data_shape.append(len(irf_dimension_2))
    else:
        data_shape.append(len(sample_irf))
    # These dimensions will be always included:
    data_shape.append(len(pedvars))
    data_shape.append(len(zds))
    data_shape.append(len(woffs))
    data = np.zeros(data_shape)
    # Iterate over all IRFs within the file. If the entry i is close to the coordinates within
    # the coord_tuple, then store.
    # tqdm for progress bars.
    for i, entry in enumerate(tqdm(eff_area_tree, total=len(all_zds))):
        if all_azs[i] != az_bin_to_store:
            continue
        if not all_indexs[i] == index_to_store:
            continue
        if irf_name == "eff":
            irf = [j for j in entry.eff]
            energies = [j for j in entry.e0]
        elif irf_name == "effNoTh2":
            irf = [j for j in entry.effNoTh2]
            energies = [j for j in entry.e0]
        elif irf_name == "hEsysMCRelative2D":
            irf = hist2array(entry.hEsysMCRelative2D)
        elif irf_name == "hEsysMCRelative2DNoDirectionCut":
            irf = hist2array(entry.hEsysMCRelative2DNoDirectionCut)
        elif irf_name == "hAngularLogDiffEmc_2D":
            irf = hist2array(entry.hAngularLogDiffEmc_2D)
        else:
            raise Exception("WrongIrfName")
        # We separate now 1D and 2D irfs:
        # 2D IRFs:
        if irf_name in implemented_irf_names_2d:
            if np.shape(irf) != np.shape(sample_irf):
                print("2D IRFs of variable size...")
                raise Exception()
            try:
                data[
                    :,
                    :,
                    find_nearest(pedvars, all_pedvars[i]),
                    find_nearest(zds, all_zds[i]),
                    find_nearest(woffs, all_Woffs[i]),
                ] = irf
            except Exception:
                print("Unexpected error:", sys.exc_info()[0])
                print("Entry number ", i)
                raise
        # 1D IRFs:
        else:
            # In case the extracted IRF has less bins than the "sample" one, pad it with zeros:
            if len(irf) < len(sample_irf):
                new_irf = np.zeros(len(sample_energies))
                for j, ener in enumerate(energies):
                    new_irf[find_nearest(sample_energies, ener)] = irf[j]
                irf = new_irf
            try:
                data[
                    :,
                    find_nearest(pedvars, all_pedvars[i]),
                    find_nearest(zds, all_zds[i]),
                    find_nearest(woffs, all_Woffs[i]),
                ] = irf
            except Exception:
                print("Unexpected error:", sys.exc_info()[0])
                print("Entry number ", i)
                print(energies, irf)
                raise

    if irf_name in implemented_irf_names_2d:
        return data, [
            np.array(irf_dimension_1),
            np.array(irf_dimension_2),
            pedvars,
            zds,
            woffs,
        ]

    return data, [
        np.array(sample_energies),
        pedvars,
        zds,
        woffs]


def duplicate_dimension(data, axis):
    # This function duplicates a single axis, assuming it's length is 1.
    current_shape = np.shape(data)
    corrected_shape = [2 if i == axis else k for i, k in enumerate(current_shape)]
    print(current_shape, corrected_shape)
    tiles = [2 if k == 1 else 1 for k in current_shape]
    new_data = np.tile(data, tiles)
    return new_data


def duplicate_dimensions(data):
    # This function duplicates all axes, one by one, when their size is 1.
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
