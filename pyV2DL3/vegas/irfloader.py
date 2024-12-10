import logging
from bisect import bisect_left

import numpy as np
import ROOT
from scipy.interpolate import RegularGridInterpolator

logger = logging.getLogger(__name__)


def graph_to_array_y(graph):
    return list(graph.GetY())


def graph_to_array_x(graph):
    return list(graph.GetX())


"""
Find closest two values for az, ze and noise axes
"""


def get_axes_edges(az, az_index, ze, ze_index, noise, noise_index):
    # Az
    for low, high in zip(az_index[:-1], az_index[1:]):
        if (az >= low) and (az < high):
            az_low = low
            az_high = high
            break
    if az > az_index[-1]:
        az_low = az_index[-1]
        az_high = az_index[0] + 360
    # Ze
    ze_low = -1
    ze_high = -1
    for low, high in zip(ze_index[:-1], ze_index[1:]):
        if (ze >= low) and (ze < high):
            ze_low = low
            ze_high = high
            break
    if (ze_low < 0) or (ze_high < 0):
        raise Exception("Ze out of range")
    # Noise
    noise_low = -1
    noise_high = -1
    for low, high in zip(noise_index[:-1], noise_index[1:]):
        if (noise >= low) and (noise < high):
            noise_low = low
            noise_high = high
            break
    if (noise_low < 0) or (noise_high < 0):
        raise Exception("Noise out of range")

    return az_low, az_high, ze_low, ze_high, noise_low, noise_high


def get_irf_not_safe(manager, offset_arr, az, ze, noise, pointlike, psf_king=False):
    loaded_offset = []
    ea_data_dict = {}
    ea = []
    ebias_data_dict = {}
    ebias = []
    abias_data_dict = {}
    abias = []
    for offset in offset_arr:
        effectiveAreaParameters = ROOT.VAEASimpleParameterData()
        effectiveAreaParameters.fAzimuth = az
        effectiveAreaParameters.fZenith = ze
        effectiveAreaParameters.fNoise = noise
        effectiveAreaParameters.fOffset = offset
        effectiveAreaParameters = manager.getVectorParamsFromSimpleParameterData(
            effectiveAreaParameters
        )
        ea_dl3 = None
        if pointlike or psf_king:
            ea_dl3 = manager.getEffectiveAreaCurve(effectiveAreaParameters)
            eb_dl3 = manager.getEnergyBias2D(effectiveAreaParameters)
        else:
            ea_dl3 = manager.getEffectiveAreaCurve_DL3_no_theta_cut(
                effectiveAreaParameters
            )
            eb_dl3 = manager.getEnergyBias_DL3(effectiveAreaParameters, False)

        if not ea_dl3:
            continue

        # Get Ebias
        n_bins_x = eb_dl3.GetNbinsX()
        n_bins_y = eb_dl3.GetNbinsY()

        bin_edges_x = [
            eb_dl3.GetXaxis().GetBinLowEdge(i) for i in range(1, n_bins_x + 2)
        ]
        bin_edges_y = [
            eb_dl3.GetYaxis().GetBinLowEdge(i) for i in range(1, n_bins_y + 2)
        ]
        a = np.zeros((n_bins_x, n_bins_y))
        for i in range(1, n_bins_x + 1):
            for j in range(1, n_bins_y + 1):
                bin_content = eb_dl3.GetBinContent(i, j)
                a[i - 1, j - 1] = bin_content
        e = np.vstack((bin_edges_x, bin_edges_y))

        eLow = np.power(10, [e[0][:-1]])[0]
        eHigh = np.power(10, [e[0][1:]])[0]

        bLow = e[1][:-1]
        bHigh = e[1][1:]
        if pointlike or psf_king:
            bLow = np.power(10, [e[1][:-1]])[0]
            bHigh = np.power(10, [e[1][1:]])[0]

        ac = []
        for aa in a:
            ab = aa / (bHigh - bLow) / np.sum(aa) if np.sum(aa) > 0 else aa
            try:
                ac = np.vstack((ac, ab))
            except ValueError:
                ac = ab

        ac = ac.transpose()
        ebias_energy_low = eLow
        ebias_energy_high = eHigh
        ebias_migration_low = bLow
        ebias_migration_high = bHigh
        ebias.append(ac)

        # Get Effective Area
        energy_bin = graph_to_array_x(ea_dl3)
        energy_bin = np.array(energy_bin)
        energy_bin = np.power(10, energy_bin)
        ea_array = graph_to_array_y(ea_dl3)
        ea_array = np.array(ea_array)
        ea_tmp = np.zeros(len(eLow))
        for e, energy in enumerate(energy_bin):
            indx = bisect_left(eLow, energy) - 1
            if (indx >= 0) and (indx < len(eLow)):
                ea_tmp[indx] = ea_array[e]
        ea.append(ea_tmp)
        ea_energy_low = eLow
        ea_energy_high = eHigh

        # Get ABias
        if not pointlike and not psf_king:
            a = np.array(
                [
                    manager.getAngularBias_DL3(effectiveAreaParameters).GetBinContent(i)
                    for i in range(
                        1,
                        manager.getAngularBias_DL3(effectiveAreaParameters).GetNbinsX()
                        + 1,
                    )
                ]
            )
            e = np.array(
                [
                    manager.getAngularBias_DL3(effectiveAreaParameters).GetBinLowEdge(i)
                    for i in range(
                        1,
                        manager.getAngularBias_DL3(effectiveAreaParameters).GetNbinsX()
                        + 2,
                    )
                ]
            )

            eLow = np.power(10, [e[0][:-1]])[0]
            eHigh = np.power(10, [e[0][1:]])[0]

            bLow = np.power(10, [e[1][:-1]])[0]
            bHigh = np.power(10, [e[1][1:]])[0]

            ac = []
            for aa in a:
                if np.sum(aa) > 0:
                    # As the unit is sr^-1 we need to convert y bin size into radian
                    ab = (
                        aa
                        / np.deg2rad(bHigh - bLow)
                        / np.sum(aa)
                        / np.pi
                        / np.deg2rad(bHigh + bLow)
                    )
                else:
                    ab = aa
                try:
                    ac = np.vstack((ac, ab))
                except ValueError:
                    ac = ab

            ac = ac.transpose()
            abias_energy_low = eLow
            abias_energy_high = eHigh
            abias_migration_low = bLow
            abias_migration_high = bHigh
            abias.append(ac)

        loaded_offset.append(offset)

    #  Create EA data
    ea = np.array(ea)
    ea_data_dict["ELow"] = ea_energy_low
    ea_data_dict["EHigh"] = ea_energy_high
    ea_data_dict["ThetaLow"] = np.array(loaded_offset)  # offset_np_arr
    ea_data_dict["ThetaHigh"] = np.array(loaded_offset)  # offset_np_arr
    ea_data_dict["Data"] = ea
    ebias = np.array(ebias)
    ebias_data_dict["ELow"] = ebias_energy_low
    ebias_data_dict["EHigh"] = ebias_energy_high
    ebias_data_dict["ThetaLow"] = np.array(loaded_offset)  # offset_np_arr
    ebias_data_dict["ThetaHigh"] = np.array(loaded_offset)  # offset_np_arr
    ebias_data_dict["MigrationLow"] = ebias_migration_low
    ebias_data_dict["MigrationHigh"] = ebias_migration_high
    ebias_data_dict["Data"] = ebias

    # PSF king fills abias dict via the king parameters (`get_psf_king()`)
    if not pointlike and not psf_king:
        abias = np.array(abias)
        abias_data_dict["ELow"] = abias_energy_low
        abias_data_dict["EHigh"] = abias_energy_high
        abias_data_dict["ThetaLow"] = np.array(loaded_offset)  # offset_np_arr
        abias_data_dict["ThetaHigh"] = np.array(loaded_offset)  # offset_np_arr
        abias_data_dict["MigrationLow"] = abias_migration_low
        abias_data_dict["MigrationHigh"] = abias_migration_high
        abias_data_dict["Data"] = abias

    return ea_data_dict, ebias_data_dict, abias_data_dict


def getIRF(az, ze, noise, event_class, pointlike, psf_king_params=None):
    axis_dict = event_class.axis_dict
    manager = event_class.manager

    # Find closest two values for az, ze and noise axis
    az_low, az_high, ze_low, ze_high, noise_low, noise_high = get_axes_edges(
        az, axis_dict["Azimuth"], ze, axis_dict["Zenith"], noise, axis_dict["Noise"]
    )

    irf_data = []
    offset_index = axis_dict["AbsoluteOffset"]
    for az_i, az_val in [(0, az_low), (1, az_high)]:
        for ze_i, ze_val in [(0, ze_low), (1, ze_high)]:
            for noise_i, noise_val in [(0, noise_low), (1, noise_high)]:
                irf_dict = {"Index": (az_i, ze_i, noise_i)}
                az_shifted = az_val if az_val <= 360 else az_val - 360.0
                ea_dict, ebias_dict, abias_dict = get_irf_not_safe(
                    manager,
                    offset_index,
                    az_shifted,
                    ze_val,
                    noise_val,
                    pointlike,
                    psf_king=psf_king_params is not None,
                )
                irf_dict["EA_Dict"] = ea_dict
                irf_dict["EBias_Dict"] = ebias_dict
                irf_dict["ABias_Dict"] = abias_dict
                irf_data.append(irf_dict)
                # Load values
    # Initialize Data container
    ea_edim = 0
    elow = []
    ehigh = []
    offset_dim = 0
    offset_low = []
    for irf in irf_data:
        ea_data_peek = irf["EA_Dict"]["Data"]
        ebias_data_peek = irf["EBias_Dict"]["Data"]
        if ea_edim < ea_data_peek.shape[1]:
            ea_edim = ea_data_peek.shape[1]
            elow = irf["EA_Dict"]["ELow"]
            ehigh = irf["EA_Dict"]["EHigh"]
        if offset_dim < ea_data_peek.shape[0]:
            offset_dim = ea_data_peek.shape[0]
            offset_low = irf["EA_Dict"]["ThetaLow"]

    ea_array = np.zeros([2, 2, 2, offset_dim, ea_edim])
    ebias_array = np.zeros(
        [
            2,
            2,
            2,
            offset_dim,
            (ebias_data_peek.shape[1]),
            (ebias_data_peek.shape[1]),
        ]
    )

    # Build Interpolator
    for irf in irf_data:
        index = irf["Index"]
        ea_data = irf["EA_Dict"]["Data"]
        ebias_data = irf["EBias_Dict"]["Data"]
        elow_first = irf["EA_Dict"]["ELow"][0]
        offset_low_first = irf["EA_Dict"]["ThetaLow"][0]

        for i, val in enumerate(elow):
            if np.abs(val - elow_first) < 1e-15:
                break
        eindex_low, eindex_high = i, i + len(irf["EA_Dict"]["ELow"])

        for i, val in enumerate(offset_low):
            if val == offset_low_first:
                break
        offset_index_low, offset_index_high = i, i + len(irf["EA_Dict"]["ThetaLow"])

        ea_array[index[0], index[1], index[2]][
            offset_index_low:offset_index_high, eindex_low:eindex_high
        ] = ea_data
        ebias_array[index[0], index[1], index[2]][
            offset_index_low:offset_index_high
        ] = ebias_data
    inter_axis = np.array(
        [[az_low, az_high], [ze_low, ze_high], [noise_low, noise_high]]
    )
    ea_interpolator = RegularGridInterpolator(inter_axis, ea_array)
    ebias_interpolator = RegularGridInterpolator(inter_axis, ebias_array)

    # Now lets actually build the data block to be passed
    # EA
    elow = irf_data[0]["EA_Dict"]["ELow"]
    ehigh = irf_data[0]["EA_Dict"]["EHigh"]
    thetalow = irf_data[0]["EA_Dict"]["ThetaLow"]
    thetahigh = irf_data[0]["EA_Dict"]["ThetaHigh"]
    ea_interpolated = ea_interpolator((az, ze, noise))
    ea_final_data = np.array(
        [(elow, ehigh, thetalow, thetahigh, ea_interpolated)],
        dtype=[
            ("ENERG_LO", ">f4", np.shape(elow)),
            ("ENERG_HI", ">f4", np.shape(ehigh)),
            ("THETA_LO", ">f4", np.shape(thetalow)),
            ("THETA_HI", ">f4", np.shape(thetahigh)),
            ("EFFAREA", ">f4", np.shape(ea_interpolated)),
        ],
    )
    # EBias
    elow = irf_data[0]["EBias_Dict"]["ELow"]
    ehigh = irf_data[0]["EBias_Dict"]["EHigh"]
    thetalow = irf_data[0]["EBias_Dict"]["ThetaLow"]
    thetahigh = irf_data[0]["EBias_Dict"]["ThetaHigh"]
    miglow = irf_data[0]["EBias_Dict"]["MigrationLow"]
    mighigh = irf_data[0]["EBias_Dict"]["MigrationHigh"]
    ebias_interpolated = ebias_interpolator((az, ze, noise))

    ebias_final_data = np.array(
        [(elow, ehigh, miglow, mighigh, thetalow, thetahigh, ebias_interpolated)],
        dtype=[
            ("ENERG_LO", ">f4", np.shape(elow)),
            ("ENERG_HI", ">f4", np.shape(ehigh)),
            ("MIGRA_LO", ">f4", np.shape(miglow)),
            ("MIGRA_HI", ">f4", np.shape(mighigh)),
            ("THETA_LO", ">f4", np.shape(thetalow)),
            ("THETA_HI", ">f4", np.shape(thetahigh)),
            ("MATRIX", ">f4", np.shape(ebias_interpolated)),
        ],
    )
    # ABias
    abias_final_data = None
    if psf_king_params is not None:
        abias_king_dict = get_king_psf_params(
            az, ze, noise, event_class, psf_king_params
        )

        elow = abias_king_dict["ELow"]
        ehigh = abias_king_dict["EHigh"]
        thetalow = abias_king_dict["ThetaLow"]
        thetahigh = abias_king_dict["ThetaHigh"]
        psf_gamma = abias_king_dict["Gamma"]
        psf_sigma = abias_king_dict["Sigma"]

        abias_final_data = np.array(
            [(elow, ehigh, thetalow, thetahigh, psf_gamma, psf_sigma)],
            dtype=[
                ("ENERG_LO", ">f4", np.shape(elow)),
                ("ENERG_HI", ">f4", np.shape(ehigh)),
                ("THETA_LO", ">f4", np.shape(thetalow)),
                ("THETA_HI", ">f4", np.shape(thetahigh)),
                ("GAMMA", ">f4", np.shape(psf_gamma)),
                ("SIGMA", ">f4", np.shape(psf_sigma)),
            ],
        )

    elif not pointlike:
        for irf in irf_data:
            abias_data_peek = irf["ABias_Dict"]["Data"]

        abias_array = np.zeros(
            [
                2,
                2,
                2,
                offset_dim,
                (abias_data_peek.shape[1]),
                (abias_data_peek.shape[2]),
            ]
        )

        # Build Interpolator
        for irf in irf_data:
            abias_data = irf["ABias_Dict"]["Data"]
            abias_array[index[0], index[1], index[2]][
                offset_index_low:offset_index_high
            ] = abias_data

        abias_interpolator = RegularGridInterpolator(inter_axis, abias_array)

        elow = irf_data[0]["ABias_Dict"]["ELow"]
        ehigh = irf_data[0]["ABias_Dict"]["EHigh"]
        thetalow = irf_data[0]["ABias_Dict"]["ThetaLow"]
        thetahigh = irf_data[0]["ABias_Dict"]["ThetaHigh"]
        miglow = irf_data[0]["ABias_Dict"]["MigrationLow"]
        mighigh = irf_data[0]["ABias_Dict"]["MigrationHigh"]
        abias_interpolated = abias_interpolator((az, ze, noise))
        # Flip axis order
        # Axis order:
        # Energy, Theta, Rad
        abias_interpolated = np.transpose(abias_interpolated, axes=(1, 0, 2))
        abias_final_data = np.array(
            [
                (
                    elow,
                    ehigh,
                    thetalow,
                    thetahigh,
                    miglow,
                    mighigh,
                    abias_interpolated,
                )
            ],
            dtype=[
                ("ENERG_LO", ">f4", np.shape(elow)),
                ("ENERG_HI", ">f4", np.shape(ehigh)),
                ("THETA_LO", ">f4", np.shape(thetalow)),
                ("THETA_HI", ">f4", np.shape(thetahigh)),
                ("RAD_LO", ">f4", np.shape(miglow)),
                ("RAD_HI", ">f4", np.shape(mighigh)),
                ("RPSF", ">f4", np.shape(abias_interpolated)),
            ],
        )

    return ea_final_data, ebias_final_data, abias_final_data


"""
Find closest PSF value for each axis from the provided indexes
"""


def get_psf_axes_values(az, az_index, ze, ze_index, noise, noise_index):
    # Az
    az_low = az_index[0]
    az_high = az_index[-1]
    for low, high in zip(az_index[:-1], az_index[1:]):
        if (az >= low) and (az < high):
            az_low = low
            az_high = high
            break
    if abs(az - az_low) > abs(az_high - az):
        az_psf = az_high
    else:
        az_psf = az_low

    # Ze
    ze_low = ze_index[0]
    ze_high = ze_index[-1]
    for low, high in zip(ze_index[:-1], ze_index[1:]):
        if (ze >= low) and (ze < high):
            ze_low = low
            ze_high = high
            break
    if abs(ze - ze_low) > abs(ze_high - ze):
        zen_psf = ze_high
    else:
        zen_psf = ze_low

    # Noise
    noise_low = noise_index[0]
    noise_high = noise_index[-1]
    for low, high in zip(noise_index[:-1], noise_index[1:]):
        if (noise >= low) and (noise < high):
            noise_low = low
            noise_high = high
            break
    if abs(noise - noise_low) > abs(noise_high - noise):
        noise_psf = noise_high
    else:
        noise_psf = noise_low

    return az_psf, zen_psf, noise_psf


"""
Fill PSF table from King PSF parameters
"""


def get_king_psf_params(az, ze, noise, event_class, psf_king_params):
    msw_lower = event_class.msw_lower
    msw_upper = event_class.msw_upper
    if msw_lower == float("-inf") or msw_upper == float("inf"):
        raise Exception(
            "--psf_king currently requires MSW cuts to be defined in your EA file. ",
            "King parameters are defined for particular MSW ranges",
        )

    psf_king_index = psf_king_params["index"]
    az_psf, zen_psf, noise_psf = get_psf_axes_values(
        az,
        psf_king_index["Azimuth"],
        ze,
        psf_king_index["Zenith"],
        noise,
        psf_king_index["Noise"],
    )

    # Init data structs
    abias_king_dict = {}
    full_sigma = []
    full_lambda = []
    offset_arrs = {}
    offset_index = psf_king_index["AbsoluteOffset"]
    for offset in offset_index:
        offset_arrs[str(offset)] = []

    # Search the king params for a line to match our arguments
    for param_values in psf_king_params["values"]:
        # `if` statements are lazily evaluated in Python
        if (
            param_values[0] == zen_psf
            and param_values[2] == noise_psf
            and param_values[3] == az_psf
            and param_values[4] == msw_lower
            and param_values[5] == msw_upper
        ):
            # Add the line to its corresponding offset array in offset_arrs
            for offset in offset_index:
                if param_values[1] == offset:
                    offset_arrs[str(offset)].append(param_values)

    for key in offset_arrs:
        offset_arr = np.array(offset_arrs[key])
        if len(offset_arr) == 0:
            raise Exception(
                "Could not find any matching values in the PSF king parameters file"
            )
        # Build sigma and lambda
        # Offset_arrs is sorted ascending because offset_index was sorted ascending
        full_sigma.append(offset_arr[:, 8])
        full_lambda.append(offset_arr[:, 10])

    # Just need one offset array to sample energy bin values
    first_offset_arr = np.array(offset_arrs[str(offset_index[0])])
    energy_low = first_offset_arr[:, 6]
    energy_high = first_offset_arr[:, 7]
    energy_low_bins = np.power(10, energy_low)
    energy_high_bins = np.power(10, energy_high)

    # Assign loaded values for return
    abias_king_dict["ELow"] = np.array(energy_low_bins)
    abias_king_dict["EHigh"] = np.array(energy_high_bins)
    abias_king_dict["ThetaLow"] = np.array(offset_index)
    abias_king_dict["ThetaHigh"] = np.array(offset_index)
    abias_king_dict["Gamma"] = np.array(full_lambda)
    abias_king_dict["Sigma"] = np.array(full_sigma)

    return abias_king_dict
