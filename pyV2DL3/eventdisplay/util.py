import uproot
from root_numpy import hist2array
import numpy as np
import sys
from ROOT import gSystem, TFile, TCanvas, TGraphAsymmErrors, TH1D, TH2D, TGraphAsymmErrors, TProfile


def produce_tel_list(tel_config):
    # Convert the list of telescopes into a string for FITS header
    tel_list = ""
    for tel in tel_config['TelID']:
        tel_list += "T" + str(tel) + ","
    return tel_list[:-1]


def is_close(a, b, rel_tol=1e-04, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


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
    return np.delete(axis + bin_size / 2., len(axis) - 1)


def bin_centers_to_edges(axis):
    # This function assumes bins of equal width
    bin_size = axis[1] - axis[0]
    extended_axis = np.insert(axis, 0, axis[0]-bin_size)
    return extended_axis + bin_size/2.


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def extract_irf(filename, irf_name, azimuth=False, coord_tuple=False,
                return_irf_axes=False, single_index=False):
    print("Extracting IRFs of type: {}".format(irf_name))
    # List of implemented IRFs
    implemented_irf_names_1d = ['eff', 'Rec_eff', 'gEffAreaNoTh2MC', 'gEffAreaNoTh2Rec']
    implemented_irf_names_2d = ['hEsysMCRelative2D', 'hEsysMCRelative2DNoDirectionCut',
                                'hAngularLogDiff_2D', 'hAngularLogDiffEmc_2D']
    # Get both the ROOT effective area TTree and the uproot one (much faster)
    eff_area_file = TFile.Open(filename)
    eff_area_tree = eff_area_file.Get("fEffArea")
    fast_eff_area = uproot.open(filename)['fEffArea']

    # Load parameters from each TTree on arrays with uproot (super fast)
    all_zds = fast_eff_area.array('ze')
    all_azs = fast_eff_area.array('az')
    # all_azMins = fast_eff_area.array('azMin')
    # all_azMaxs = fast_eff_area.array('azMax')
    all_Woffs = fast_eff_area.array('Woff')
    all_pedvars = fast_eff_area.array('pedvar')
    all_indexs = fast_eff_area.array('index')
    all_nbins = fast_eff_area.array('nbins')
    all_rec_nbins = fast_eff_area.array('Rec_nbins')

    # If no coord_tuple is provided, extract the IRF over all dimensions
    azs = indexs = pedvars = zds = woffs = []
    if not coord_tuple:
        zds, zes_counts = np.unique(np.round(fast_eff_area.array('ze'), decimals=2), return_counts=True)
        azs = np.unique(fast_eff_area.array('az'))
        azMins, azMins_counts = np.unique(np.round(fast_eff_area.array('azMin'), decimals=2), return_counts=True)
        azMaxs, azMaxs_counts = np.unique(np.round(fast_eff_area.array('azMax'), decimals=2), return_counts=True)
        woffs, Woffs_counts = np.unique(np.round(fast_eff_area.array('Woff'), decimals=2), return_counts=True)
        pedvars, pedvars_counts = np.unique(np.round(fast_eff_area.array('pedvar'), decimals=2), return_counts=True)
        indexs, indexs_counts = np.unique(np.round(fast_eff_area.array('index'), decimals=2), return_counts=True)
        # IMPORTANT: Remove duplicities (shouldn't be the case, but the values are not stored properly...)
        # We remove the duplicities from the zenith and pedestal values arrays:
        zds = remove_duplicities(zds, 2.0)
        pedvars = remove_duplicities(pedvars, 0.21)
        if len(all_zds) != len(zds) * len(azs) * len(woffs) * len(pedvars) * len(indexs):
            raise ValueError("Wrong dimensions extracted from IRF cube." +
                             "Probably due to the rounding applied to the IRF coordinates.")
    # If a specific coord_tuple is provided, only extract the IRFs wihin that range of dimensions
    else:
        if len(coord_tuple) != 5:
            raise ValueError("coord_tuple needs to contain 5 dimensions, in this order: az, index, pedvar, zd and woff")
        else:
            # Get the coordinates to sample:
            azs = coord_tuple[0]
            indexs = coord_tuple[1]
            pedvars = coord_tuple[2]
            zds = coord_tuple[3]
            woffs = coord_tuple[4]

    print(azs, indexs, pedvars, zds, woffs)

    # For performance, deactivate all branches except the ones needed:
    # Also get the entry with max bins to define the binning in energy
    eff_area_tree.SetBranchStatus("*", 0)
    if irf_name == 'eff':
        eff_area_tree.SetBranchStatus("e0", 1)
        eff_area_tree.SetBranchStatus("eff", 1)
        entry_with_max_bins = find_nearest(all_nbins, all_nbins.max())
    elif irf_name == 'Rec_eff':
        eff_area_tree.SetBranchStatus("Rec_e0", 1)
        eff_area_tree.SetBranchStatus("Rec_eff", 1)
        entry_with_max_bins = find_nearest(all_rec_nbins, all_rec_nbins.max())
    elif irf_name == 'gEffAreaNoTh2MC':
        eff_area_tree.SetBranchStatus("gEffAreaNoTh2MC", 1)
        entry_with_max_bins = find_nearest(all_nbins, all_nbins.max())
    elif irf_name == 'gEffAreaNoTh2Rec':
        eff_area_tree.SetBranchStatus("gEffAreaNoTh2Rec", 1)
        entry_with_max_bins = find_nearest(all_rec_nbins, all_rec_nbins.max())
    elif irf_name == 'hEsysMCRelative2D':
        eff_area_tree.SetBranchStatus("hEsysMCRelative2D", 1)
        entry_with_max_bins = 0
    elif irf_name == 'hEsysMCRelative2DNoDirectionCut':
        eff_area_tree.SetBranchStatus("hEsysMCRelative2DNoDirectionCut", 1)
        entry_with_max_bins = 0
    elif irf_name == 'hAngularLogDiff_2D':
        eff_area_tree.SetBranchStatus("hAngularLogDiff_2D", 1)
        entry_with_max_bins = 0
    elif irf_name == 'hAngularLogDiffEmc_2D':
        eff_area_tree.SetBranchStatus("hAngularLogDiffEmc_2D", 1)
        entry_with_max_bins = 0
    else:
        raise Exception("WrongIrfName")
    # Now we know which entry we need to get in order to have a sample IRF
    #     sample_irf = sample_energies = []
    for i, entry in enumerate(eff_area_tree):
        if i == entry_with_max_bins:
            if irf_name == 'eff':
                sample_irf = [j for j in entry.eff]
                sample_energies = [j for j in entry.e0]
            elif irf_name == 'Rec_eff':
                sample_irf = [j for j in entry.Rec_eff]
                sample_energies = [j for j in entry.Rec_e0]
            elif irf_name == 'effNoTh2':
                sample_irf = [j for j in entry.eff]
                sample_energies = [j for j in entry.e0]
            elif irf_name == 'Rec_effNoTh2':
                sample_irf = [j for j in entry.Rec_eff]
                sample_energies = [j for j in entry.Rec_e0]
            elif irf_name == 'hEsysMCRelative2D':
                # Migration vs energy bias and true energy
                sample_irf, axes = hist2array(entry.hEsysMCRelative2D, return_edges=True)
                # Bin edges (one more entry than migra!) for the true energy and
                # energy bias (Erec/Etrue)
                irf_dimension_1 = bin_edges_to_centers(axes[0])
                irf_dimension_2 = bin_edges_to_centers(axes[1])
            elif irf_name == 'hEsysMCRelative2DNoDirectionCut':
                # Migration vs energy bias and true energy, without direction cut
                sample_irf, axes = hist2array(entry.hEsysMCRelative2DNoDirectionCut, return_edges=True)
                # Bin edges (one more entry than migra!) for the true energy and
                # energy bias (Erec/Etrue)
                irf_dimension_1 = bin_edges_to_centers(axes[0])
                irf_dimension_2 = bin_edges_to_centers(axes[1])
            elif irf_name == 'hAngularLogDiff_2D':
                # PSF vs true energy:
                sample_irf, axes = hist2array(entry.hAngularLogDiff_2D, return_edges=True)
                # Bin edges (one more entry than migra!) for the true energy and
                # energy bias (Erec/Etrue)
                irf_dimension_1 = bin_edges_to_centers(axes[0])
                irf_dimension_2 = bin_edges_to_centers(axes[1])
            elif irf_name == 'hAngularLogDiffEmc_2D':
                # PSF vs true energy:
                sample_irf, axes = hist2array(entry.hAngularLogDiffEmc_2D, return_edges=True)
                # Bin edges (one more entry than migra!) for the true energy and
                # energy bias (Erec/Etrue)
                irf_dimension_1 = bin_edges_to_centers(axes[0])
                irf_dimension_2 = bin_edges_to_centers(axes[1])
            else:
                raise Exception("WrongIrfName")
        if i > entry_with_max_bins:
            break
    # Now we should know all the dimensions that need to be stored in the output irf_data
    # * If an azimuth was given, only store the IRF for the closest value (remove that dimension)
    # * If 'single_index' is True, then only store one index value (average of the simulated ones)
    az_bin_to_store = 0
    if azimuth:
        azMaxs[azMaxs > 180] = 180
        azMins[azMins < -180] = -180
        az_centers = [(azMaxs[i]+azMins[i])/2. for i in np.arange(len(azMaxs))]
        az_bin_to_store = find_nearest(az_centers, azimuth)
    if single_index:
        # Find the closest index to the average value simulated:
        index_to_store = indexs[find_nearest(indexs, (indexs.min()+indexs.max())/2.)]
    # Create data container, filled with zeros, containing the required dimensions to store
    # the IRF for a given coord_tuple. Separated between 1 and 2 dimensions:
    data_shape = []
    if irf_name in implemented_irf_names_2d:
        data_shape.append(len(irf_dimension_1))
        data_shape.append(len(irf_dimension_2))
    else:
        data_shape.append(len(sample_irf))
    # If an azimuth was provided, this dimension will not be included
    if not azimuth:
        data_shape.append(len(azs))
    # This dimension will not be included if we only want to store one single index
    if not single_index:
        data_shape.append(len(indexs))
    # These dimensions will be always included:
    data_shape.append(len(pedvars))
    data_shape.append(len(zds))
    data_shape.append(len(woffs))
    data = np.zeros(data_shape)
    print(np.shape(data))
    # Iterate over all IRFs within the file. If the entry i is close to the coordinates within
    # the coord_tuple, then store.
    for i, entry in enumerate(eff_area_tree):
        if i % 5000 == 0:
            print("{}/{}".format(i, len(all_zds)))
        # Parameters within the effective area files show some fluctuation, therefore
        # we need to use the "isclose".
        if (np.isclose(azs, all_azs[i], atol=0.01).any() and
                np.isclose(indexs, all_indexs[i], atol=0.01).any() and
                np.isclose(pedvars, all_pedvars[i], atol=0.21).any() and
                np.isclose(zds, all_zds[i], atol=2.1).any() and
                np.isclose(woffs, all_Woffs[i], atol=0.05).any()):
            # If an azimuth value is given, only store closest azimuth bin
            if azimuth:
                if all_azs[i] != az_bin_to_store:
                    continue
            # If 'single_index' flag is given, also ignore other index values
            if single_index:
                if not np.isclose(all_indexs[i], index_to_store, atol=0.01):
                    continue
            if irf_name == 'eff':
                irf = [j for j in entry.eff]
                energies = [j for j in entry.e0]
            elif irf_name == 'Rec_eff':
                irf = [j for j in entry.Rec_eff]
                energies = [j for j in entry.Rec_e0]
            elif irf_name == 'gEffAreaNoTh2MC':
                energies, irf = graph_to_array(entry.gEffAreaNoTh2MC, all_nbins[i])
                # irf = graph_to_array_y(entry.gEffAreaNoTh2MC)
                # energies = graph_to_array_x(entry.gEffAreaNoTh2MC)
            elif irf_name == 'gEffAreaNoTh2Rec':
                energies, irf = graph_to_array(entry.gEffAreaNoTh2Rec, all_rec_nbins[i])
                # irf = graph_to_array_y(entry.gEffAreaNoTh2Rec)
                # energies = graph_to_array_x(entry.gEffAreaNoTh2Rec)
            elif irf_name == 'hEsysMCRelative2D':
                irf = hist2array(entry.hEsysMCRelative2D)
            elif irf_name == 'hEsysMCRelative2DNoDirectionCut':
                irf = hist2array(entry.hEsysMCRelative2DNoDirectionCut)
            elif irf_name == 'hAngularLogDiff_2D':
                irf = hist2array(entry.hAngularLogDiff_2D)
            elif irf_name == 'hAngularLogDiffEmc_2D':
                irf = hist2array(entry.hAngularLogDiffEmc_2D)
            else:
                raise Exception("WrongIrfName")
            # We separate now 1D and 2D irfs:
            # 2D IRFs:
            if irf_name in implemented_irf_names_2d:
                if np.shape(irf) == np.shape(sample_irf):
                    # I'm CERTAIN this part can be done more elegantly... For now, good enough
                    if azimuth and single_index:
                        try:
                            data[:, :, find_nearest(pedvars, all_pedvars[i]), find_nearest(zds, all_zds[i]),
                                 find_nearest(woffs, all_Woffs[i])] = irf
                        except Exception:
                            print("Unexpected error:", sys.exc_info()[0])
                            print("Entry number ", i)
                            raise
                    elif azimuth:
                        try:
                            data[:, :, find_nearest(indexs, all_indexs[i]), find_nearest(pedvars, all_pedvars[i]),
                                 find_nearest(zds, all_zds[i]), find_nearest(woffs, all_Woffs[i])] = irf
                        except Exception:
                            print("Unexpected error:", sys.exc_info()[0])
                            print("Entry number ", i)
                            raise
                    elif single_index:
                        try:
                            data[:, :, find_nearest(azs, all_azs[i]), find_nearest(pedvars, all_pedvars[i]),
                                 find_nearest(zds, all_zds[i]), find_nearest(woffs, all_Woffs[i])] = irf
                        except Exception:
                            print("Unexpected error:", sys.exc_info()[0])
                            print("Entry number ", i)
                            raise
                    else:
                        try:
                            data[:, :, find_nearest(azs, all_azs[i]), find_nearest(indexs, all_indexs[i]),
                                 find_nearest(pedvars, all_pedvars[i]),
                                 find_nearest(zds, all_zds[i]), find_nearest(woffs, all_Woffs[i])] = irf
                        except Exception:
                            print("Unexpected error:", sys.exc_info()[0])
                            print("Entry number ", i)
                            raise
                else:
                    print("2D IRFs of variable size...")
                    raise Exception()
            # 1D IRFs:
            else:
                # In case the extracted IRF has less bins than the "sample" one, pad it with zeros:
                if len(irf) < len(sample_irf):
                    new_irf = np.zeros(len(sample_energies))
                    for j, ener in enumerate(energies):
                        new_irf[find_nearest(sample_energies, ener)] = irf[j]
                    irf = new_irf
                if azimuth and single_index:
                    try:
                        data[:, find_nearest(pedvars, all_pedvars[i]), find_nearest(zds, all_zds[i]),
                             find_nearest(woffs, all_Woffs[i])] = irf
                    except Exception:
                        print("Unexpected error:", sys.exc_info()[0])
                        print("Entry number ", i)
                        print(energies, irf)
                        raise
                elif azimuth:
                    try:
                        data[:, find_nearest(indexs, all_indexs[i]), find_nearest(pedvars, all_pedvars[i]),
                         find_nearest(zds, all_zds[i]), find_nearest(woffs, all_Woffs[i])] = irf
                    except Exception:
                        print("Unexpected error:", sys.exc_info()[0])
                        print("Entry number ", i)
                        print(energies, irf)
                        raise
                elif single_index:
                    try:
                        data[:, find_nearest(azs, all_azs[i]), find_nearest(pedvars, all_pedvars[i]),
                             find_nearest(zds, all_zds[i]), find_nearest(woffs, all_Woffs[i])] = irf
                    except Exception:
                        print("Unexpected error:", sys.exc_info()[0])
                        print("Entry number ", i)
                        print(energies, irf)
                        raise
                else:
                    try:
                        data[:, find_nearest(azs, all_azs[i]), find_nearest(indexs, all_indexs[i]),
                         find_nearest(pedvars, all_pedvars[i]),
                         find_nearest(zds, all_zds[i]), find_nearest(woffs, all_Woffs[i])] = irf
                    except Exception:
                        print("Unexpected error:", sys.exc_info()[0])
                        print("Entry number ", i)
                        print(energies, irf)
                        raise
    if return_irf_axes:
        if irf_name in implemented_irf_names_2d:
            if azimuth and single_index:
                return data, [np.array(irf_dimension_1), np.array(irf_dimension_2), pedvars, zds, woffs]
            elif azimuth:
                return data, [np.array(irf_dimension_1), np.array(irf_dimension_2), indexs, pedvars, zds, woffs]
            elif single_index:
                return data, [np.array(irf_dimension_1), np.array(irf_dimension_2), azs, pedvars, zds, woffs]
            else:
                return data, [np.array(irf_dimension_1), np.array(irf_dimension_2), azs, indexs, pedvars, zds, woffs]
        else:
            if azimuth and single_index:
                return data, [np.array(sample_energies), pedvars, zds, woffs]
            elif azimuth:
                return data, [np.array(sample_energies), indexs, pedvars, zds, woffs]
            elif single_index:
                return data, [np.array(sample_energies), azs, pedvars, zds, woffs]
            else:
                return data, [np.array(sample_energies), azs, indexs, pedvars, zds, woffs]
    else:
        return data


def duplicate_dimension(data, axis):
    # This function duplicates a single axis, assuming it's length is 1.
    current_shape = np.shape(data)
    corrected_shape = [2 if i == axis else k for i, k in enumerate(current_shape)]
    print(current_shape, corrected_shape)
    tiles = [ 2 if k==1 else 1 for k in current_shape ]
    new_data = np.tile(data,tiles)
    return new_data


def duplicate_dimensions(data):
    # This function duplicates all axes, one by one, when their size is 1.
    new_data = data
    current_shape = np.shape(new_data)
    for i, dim in enumerate(current_shape):
        if dim == 1:
            new_data = duplicate_dimension(new_data, i)
    return new_data
