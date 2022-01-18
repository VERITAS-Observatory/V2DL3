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
    return (np.abs(array - value)).argmin()


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
        fast_eff_area,
        az_mask=None):

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
    all_nbins):
    """read energy and irf axes with maximum number of bins"""

    entry_with_max_bins = find_nearest(all_nbins, all_nbins.max())
    sample_energies = fast_eff_area['e0'].array(library="np",
                                    entry_start=entry_with_max_bins,
                                    entry_stop=entry_with_max_bins+1)[0]
    sample_irf = fast_eff_area[irf_name].array(library="np",
                                    entry_start=entry_with_max_bins,
                                    entry_stop=entry_with_max_bins+1)[0]
    return sample_energies, sample_irf, entry_with_max_bins


def read_2d_histograms(
    irf_name,
    eff_area_tree,
    entry_with_max_bins):
    """read irfaxis and sample irfs for given entry"""
    irf_dimension_1 = irf_dimension_2 = None
    axes = None

    for i, entry in enumerate(eff_area_tree):
        if i == entry_with_max_bins:
            if irf_name == "hEsysMCRelative2D":
                # Migration vs energy bias and true energy
                sample_irf, axes = hist2array(
                    entry.hEsysMCRelative2D, return_edges=True
                )
            elif irf_name == "hEsysMCRelative2DNoDirectionCut":
                # Migration vs energy bias and true energy, without direction cut
                sample_irf, axes = hist2array(
                    entry.hEsysMCRelative2DNoDirectionCut, return_edges=True
                )
            elif irf_name == "hAngularLogDiffEmc_2D":
                # PSF vs true energy:
                sample_irf, axes = hist2array(
                    entry.hAngularLogDiffEmc_2D, return_edges=True
                )
        if i > entry_with_max_bins:
            # Only extract the above defined IRF-name-specific "entry with max bins"
            break

    # Bin edges (one more entry than migra!) for the true energy and
    # energy bias (Erec/Etrue)
    if axes:
        irf_dimension_1 = bin_edges_to_centers(axes[0])
        irf_dimension_2 = bin_edges_to_centers(axes[1])

    return irf_dimension_1, irf_dimension_2, sample_irf


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

def extract_irf_1d(
    filename,
    irf_name,
    azimuth=None):
    """
    extract 1D IRF from effective area file and return a
    multidimensional array 
    """

    fast_eff_area = uproot.open(filename)["fEffAreaH2F"]

    # select az bin and define az mask
    _, azMaxs = load_parameter( "azMax", fast_eff_area)
    _, azMins = load_parameter( "azMin", fast_eff_area)
    az_bin_to_store = find_closest_az(azimuth, azMins, azMaxs)
    az_mask = fast_eff_area['az'].array(library="np") == az_bin_to_store

    energies = fast_eff_area["e0"].array(library="np")[az_mask]
    irf = fast_eff_area[irf_name].array(library="np")[az_mask]

    all_zds, zds = load_parameter( "ze", fast_eff_area, az_mask)
    all_Woffs, woffs = load_parameter( "Woff", fast_eff_area, az_mask)
    all_pedvars, pedvars = load_parameter( "pedvar", fast_eff_area, az_mask)

    data_shape = []
    data_shape.append(len(irf[0]))
    data_shape.append(len(pedvars))
    data_shape.append(len(zds))
    data_shape.append(len(woffs))
    data = np.zeros(data_shape)

    for i in range(len(irf)):
        try:
           data[
               :,
               find_nearest(pedvars, all_pedvars[i]),
               find_nearest(zds, all_zds[i]),
               find_nearest(woffs, all_Woffs[i]),
            ] = irf[i]
        except Exception:
            print("Unexpected error:", sys.exc_info()[0])
            print("Entry number ", i)
            raise

    return data, [
        energies[0],
        pedvars,
        zds,
        woffs]

def read_irf_axis(
        coordinate,
        fast_eff_area,
        irf_name,
        az_mask ):
    """return irf axis"""

    nbins = fast_eff_area[irf_name+"_bins"+coordinate].array(library="np")[az_mask]
    c_min = fast_eff_area[irf_name+"_min"+coordinate].array(library="np")[az_mask]
    c_max = fast_eff_area[irf_name+"_max"+coordinate].array(library="np")[az_mask]

    return np.linspace(c_min[0], c_max[0], nbins[0])



def extract_irf_2d(
    filename,
    irf_name,
    azimuth=None):
    """
    extract 2D IRF from effective area file and return a
    multidimensional array 
    """

    fast_eff_area = uproot.open(filename)["fEffAreaH2F"]

    # select az bin and define az mask
    _, azMaxs = load_parameter( "azMax", fast_eff_area)
    _, azMins = load_parameter( "azMin", fast_eff_area)
    az_bin_to_store = find_closest_az(azimuth, azMins, azMaxs)
    az_mask = fast_eff_area['az'].array(library="np") == az_bin_to_store

    # prepare 2D IRF
    irf_dimension_1 = read_irf_axis( "x", fast_eff_area, irf_name, az_mask )
    irf_dimension_2 = read_irf_axis( "y", fast_eff_area, irf_name, az_mask )

    print(irf_dimension_1)
    print(irf_dimension_2)
    irf1D = fast_eff_area[irf_name+"_value"].array(library="np")[az_mask]

    # parameter space
    all_zds, zds = load_parameter( "ze", fast_eff_area)
    all_Woffs, woffs = load_parameter( "Woff", fast_eff_area)
    all_pedvars, pedvars = load_parameter( "pedvar", fast_eff_area)

    return 

    data_shape = []
#    data_shape.append(len(irf_dimension_1)
#    data_shape.append(len(irf_dimension_2)
    data_shape.append(len(pedvars))
    data_shape.append(len(zds))
    data_shape.append(len(woffs))
    data = np.zeros(data_shape)

    for i in range(len(irf1D)):
        irf = np.reshape(irf1D[i], (-1, len(irf_dimension_2)))
        print("DATAAAA", irf.shape)
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

    print("IRF1", irf_dimension_1[0])
    print("IRF2", irf_dimension_2[0])
    print("pedvars", len(pedvars))
    print("zds", len(zds))
    print("woffs", len(woffs))

    return data, [
        np.array(irf_dimension_1[0]),
        np.array(irf_dimension_2[0]),
        pedvars,
        zds,
        woffs]


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

    if irf_name in implemented_irf_names_1d:
        return extract_irf_1d( filename,
                               irf_name,
                               azimuth )

#    return extract_irf_2d( filename,
#   extract_irf_2d( filename,
#                           irf_name,
#                           azimuth )

    # Get both the ROOT effective area TTree and the uproot one (much faster)
    eff_area_file = TFile.Open(filename)
    eff_area_tree = eff_area_file.Get("fEffArea")
    fast_eff_area = uproot.open(filename)["fEffArea"]

    all_zds, zds = load_parameter( "ze", fast_eff_area)
    all_azs, azs = load_parameter( "az", fast_eff_area)
    _, azMins = load_parameter( "azMin", fast_eff_area)
    _, azMaxs = load_parameter( "azMax", fast_eff_area)
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

    # Now we should know all the dimensions that need to be stored in the output irf_data
    # * For a given azimuth sore the IRF for the closest value (remove that dimension)
    az_bin_to_store = find_closest_az(azimuth, azMins, azMaxs)
    index_to_store = find_closest_index(indexs, irf_name)
    print('index to store:', index_to_store, indexs.min())

    data_shape = []

    # 2D histograms
    if irf_name in implemented_irf_names_2d:
        irf_dimension_1, irf_dimension_2, sample_irf = read_2d_histograms(
                                              irf_name,
                                              eff_area_tree,
                                              entry_with_max_bins)
        data_shape.append(len(irf_dimension_1))
        data_shape.append(len(irf_dimension_2))

    # These dimensions will be always included:
    data_shape.append(len(pedvars))
    data_shape.append(len(zds))
    data_shape.append(len(woffs))
    data = np.zeros(data_shape)
    # Iterate over all IRFs within the file. If the entry i is close to the coordinates within
    # the coord_tuple, then store.
    # tqdm for progress bars.
    print('AZS', all_azs, az_bin_to_store)
    print('INDEX', all_indexs, index_to_store)
    # 2D IRFs:
    for i, entry in enumerate(tqdm(eff_area_tree, total=len(all_zds))):
        if all_azs[i] != az_bin_to_store:
            continue
        if not all_indexs[i] == index_to_store:
            continue

        if irf_name == "hEsysMCRelative2D":
            irf = hist2array(entry.hEsysMCRelative2D)
        elif irf_name == "hEsysMCRelative2DNoDirectionCut":
            irf = hist2array(entry.hEsysMCRelative2DNoDirectionCut)
        elif irf_name == "hAngularLogDiffEmc_2D":
            irf = hist2array(entry.hAngularLogDiffEmc_2D)

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

    print("IRF1", len(np.array(irf_dimension_1)))
    print("IRF2", len(np.array(irf_dimension_2)))
    return data, [
        np.array(irf_dimension_1),
        np.array(irf_dimension_2),
        pedvars,
        zds,
        woffs,
    ]
