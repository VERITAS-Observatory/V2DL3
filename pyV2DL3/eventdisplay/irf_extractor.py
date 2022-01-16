import numpy as np
from ROOT import TFile
import sys
from tqdm.auto import tqdm
import uproot


def remove_duplicities(array, atol):
    i = 0
    while i < len(array) - 1:
        i += 1
        if np.isclose(array[i - 1], array[i], atol=atol):
            array = np.delete(array, i - 1)
            i -= 1
    return array


def bin_edges_to_centers(axis):
    # This function assumes bins of equal width
    bin_size = axis[1] - axis[0]
    return np.delete(axis + bin_size / 2.0, len(axis) - 1)


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

