from bisect import bisect_left
from ctypes import c_float
import logging

import numpy as np
import ROOT
from root_numpy import hist2array
from scipy.interpolate import RegularGridInterpolator

from pyV2DL3.vegas.load_vegas import VEGASStatus

logger = logging.getLogger(__name__)


def graph_to_array_y(graph):
    return list(graph.GetY())


def graph_to_array_x(graph):
    return list(graph.GetX())


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
            try:
                ea_dl3 = manager.getEffectiveAreaCurve_DL3_no_theta_cut(
                    effectiveAreaParameters
                )
                eb_dl3 = manager.getEnergyBias_DL3(effectiveAreaParameters, False)
            # This will append a clearer message in the event of AttributeError to let 
            # the user know why these methods are not being found.
            except AttributeError as e:
                raise Exception(str(e) 
                    + "\nFull-enclosure is unsupported without --king_function "
                    + "See exception above.")

        if not ea_dl3:
            continue

        # Get Ebias
        a, e = hist2array(eb_dl3, return_edges=True)
        eLow = np.power(10, [e[0][:-1]])[0]
        eHigh = np.power(10, [e[0][1:]])[0]

        if pointlike or psf_king:
            bLow = np.power(10, [e[1][:-1]])[0]
            bHigh = np.power(10, [e[1][1:]])[0]
        else:
            bLow = e[1][:-1]
            bHigh = e[1][1:]

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
            try:
                a, e = hist2array(
                    manager.getAngularBias_DL3(effectiveAreaParameters), return_edges=True
                )
            # This will append a clearer message in the event of AttributeError to let 
            # the user know why these methods are not being found.
            except AttributeError as e:
                raise Exception(str(e) 
                    + "\nFull-enclosure is unsupported without --king_function. "
                    + "See exception above.")

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

    if not pointlike:
        # psf_king fills this through get_psf_king()
        if not psf_king:
            abias = np.array(abias)
            abias_data_dict["ELow"] = abias_energy_low
            abias_data_dict["EHigh"] = abias_energy_high
            abias_data_dict["ThetaLow"] = np.array(loaded_offset)  # offset_np_arr
            abias_data_dict["ThetaHigh"] = np.array(loaded_offset)  # offset_np_arr
            abias_data_dict["MigrationLow"] = abias_migration_low
            abias_data_dict["MigrationHigh"] = abias_migration_high
            abias_data_dict["Data"] = abias

    return ea_data_dict, ebias_data_dict, abias_data_dict


"""
Fill PSF from king function params provided with a king params file (-k <kingfile>)

The `_psf` parameters for this function should be already assigned to a bin value of the file.

See https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#PSF_King_Parameters_File_Format
for the format of the king function parameter lines.
"""
def get_psf_king(psf_king_filename, offset_index, zen_psf, az_psf, noise_psf, msw_range=None):  
    if (msw_range is None 
        or msw_range[0] == float('-inf') 
        or msw_range[1] == float('inf')):
        raise Exception("psf_king currently requires msw ranges to be defined in your EA cuts")

    # Dev exception
    if psf_king_filename is None or psf_king_filename == "":
        raise Exception("get_psf_king was called without a kingfunction file!")

    # Init data structs
    abias_king_dict = {}
    full_sigma = []
    full_lambda = []
    offset_arrs = {}
    for offset in offset_index:
        offset_arrs[str(offset)] = []

    # Search the king params for a line to match our arguments
    with open(psf_king_filename, "r") as king_file:
        for line in king_file:
            line_arr = line.split()
            line_arr = [float(ele) for ele in line_arr]
            # `if` statements are lazily evaluated in Python, so this is efficient.
            if (line_arr[0] == zen_psf
            and line_arr[2] == noise_psf
            and line_arr[3] == az_psf
            and line_arr[4] == msw_range[0]
            and line_arr[5] == msw_range[1]):
                # Add the line to its corresponding offset array in offset_arrs
                for offset in offset_index:
                    if line_arr[1] == offset:
                        offset_arrs[str(offset)].append(line_arr)


    for key in offset_arrs:
        offset_arr = np.array(offset_arrs[key])
        if len(offset_arr) == 0:
            raise Exception("Could not find any matching params in the king params file!")
        # Build sigma and lambda
        # Offset_arrs is sorted ascending because offset_index was sorted ascending
        full_sigma.append(offset_arr[:, 8])
        full_lambda.append(offset_arr[:, 10])

    # Just need one offset array to sample our energy bin values
    first_offset_arr = np.array(offset_arrs[str(offset_index[0])])
    energy_low = first_offset_arr[:, 6]
    energy_high = first_offset_arr[:, 7]
    energy_low_bins = np.power(10, energy_low)
    energy_high_bins = np.power(10, energy_high)
    
    # Assign our loaded values for return
    abias_king_dict['ELow'] = np.array(energy_low_bins)
    abias_king_dict['EHigh'] = np.array(energy_high_bins)
    abias_king_dict['ThetaLow'] = np.array(offset_index)
    abias_king_dict['ThetaHigh'] = np.array(offset_index)
    abias_king_dict['Gamma'] = np.array(full_lambda)
    abias_king_dict['Sigma'] = np.array(full_sigma)

    return abias_king_dict


class IRFLoader:
    def __init__(self, vts_io, irf_to_store, event_class_idx=0, multi_eclass=False):
        # If vegas not loaded. Load vegas
        self.__vegas__ = VEGASStatus()
        self.__vegas__.loadVEGAS()
        self.__manager__ = ROOT.VAEffectiveAreaManager()
        self.__manager__.setUseReconstructedEnergy(False)
        self.__manager__.loadEffectiveAreas(vts_io)
        # This will stay 0 (default) in the case of a single event class
        self.__event_class_idx__ = event_class_idx
        # We will only append a column identifying the event class if there are multiple event classes
        self.__multi_eclass__ = multi_eclass
        # Deal with AbsoluteOffset separately
        self.__axis__ = ["Azimuth", "Zenith", "Noise"]  # , 'AbsoluteOffset']
        # Init default control flows
        self.__pointlike__ = False           
        self.__psf_king_filename__ = None
        self.__psf_king__ = False
        # Load any provided control flows
        if "point-like" in irf_to_store:
            self.__pointlike__ = irf_to_store["point-like"]
        if "psf-king" in irf_to_store:
            self.__psf_king__ = irf_to_store["psf-king"]
            if "psf-king-filename" in irf_to_store:
                self.__psf_king_filename__ = irf_to_store["psf-king-filename"]
            
        self.__buildIndex__()

    def __buildIndex__(self):
        manager = self.__manager__
        if len(manager.fEffectiveAreas) <= 0:
            raise Exception("No effective areas! ")
        index_check = manager.fEffectiveAreas.at(0).fDimensionNames
        for k in self.__axis__:
            if k not in index_check:
                raise Exception("IRF missing axis: {}".format(k))
        index_dict = {"Index": []}

        for i, ea in enumerate(manager.fEffectiveAreas):
            index_dict["Index"].append(i)
            for name, val in zip(ea.fDimensionNames, ea.fDimensionValues):
                if name not in index_dict:
                    index_dict[name] = []
                else:
                    index_dict[name].append(val)

        # Deal with AbsoluteOffset
        if "AbsoluteOffset" not in index_check:
            logger.info("No offset axis available from file. Use 0.5 deg as default.")
            index_dict["AbsoluteOffset"] = []
            for _ in range(len(index_dict["Index"])):
                index_dict["AbsoluteOffset"].append(0.5)

        # Validate Completeness
        axis_dict = {}
        check_num = 1

        for k in self.__axis__ + ["AbsoluteOffset"]:
            check_num *= len(np.unique(index_dict[k]))
            axis_dict[k] = np.sort(np.unique(index_dict[k]))
            if len(axis_dict[k]) < 2 and k != "AbsoluteOffset":
                raise Exception("{} Axis need to have more than two values".format(k))
        self.__axis_dict__ = axis_dict
        self.__index_dict__ = index_dict

        # If using king function, determine bin values from the king params file
        if self.__psf_king__:
            print("Loading PSF bins...")
            psf_king_index = {}
            index_keys = ["Zenith", "AbsoluteOffset", "Noise", "Azimuth"]
            # Initialize dict from keys
            for key in index_keys:
                psf_king_index[key] = []
            # Add unique bin values
            with open(self.__psf_king_filename__, "r") as king_file:
                for line in king_file:
                    line_arr = line.split()
                    if line_arr[0] not in psf_king_index["Zenith"]:
                        psf_king_index["Zenith"].append(line_arr[0])
                    if line_arr[1] not in psf_king_index["AbsoluteOffset"]:
                        psf_king_index["AbsoluteOffset"].append(line_arr[1])
                    if line_arr[2] not in psf_king_index["Noise"]:
                        psf_king_index["Noise"].append(line_arr[2])
                    if line_arr[3] not in psf_king_index["Azimuth"]:
                        psf_king_index["Azimuth"].append(line_arr[3])
            # Cast the strings to floats and sort ascending
            for key in psf_king_index:
                psf_king_index[key] = [float(ele) for ele in psf_king_index[key]]
                psf_king_index[key].sort()
            # Log results to debug mode
            logger.debug("PSF king bins loaded:")
            for key in psf_king_index:
                logger.debug(key + ": " + str(psf_king_index[key]))

            self.__psf_king_index__ = psf_king_index


    def getSafeEnergy(self, az, ze, noise):
        manager = self.__manager__
        effectiveAreaParameters = ROOT.VAEASimpleParameterData()
        effectiveAreaParameters.fAzimuth = az
        effectiveAreaParameters.fZenith = ze
        effectiveAreaParameters.fNoise = noise

        effectiveAreaParameters.fOffset = 0.5
        effectiveAreaParameters = manager.getVectorParamsFromSimpleParameterData(
            effectiveAreaParameters
        )
        minEnergy, maxEnergy = c_float(), c_float()
        # Is it the right way ? what does the offset here provide ?
        manager.getSafeEnergyRange(effectiveAreaParameters, 0.5, minEnergy, maxEnergy)
        return minEnergy.value / 1000.0, maxEnergy.value / 1000.0

    def getIRF(self, az, ze, noise, msw_range=None):
        # Find closest two values for az, ze and noise axis
        # Az
        az_index = self.__axis_dict__["Azimuth"]
        if self.__psf_king__:
            az_index = self.__psf_king_index__["Azimuth"]
        for low, high in zip(az_index[:-1], az_index[1:]):
            if (az >= low) and (az < high):
                az_low = low
                az_high = high
                break
        if az > az_index[-1]:
            az_low = az_index[-1]
            az_high = az_index[0] + 360
        # Ze
        ze_index = self.__axis_dict__["Zenith"]
        ze_low = -1
        ze_high = -1
        for low, high in zip(ze_index[:-1], ze_index[1:]):
            if (ze >= low) and (ze < high):
                ze_low = low
                ze_high = high
                break
        if (ze_low < 0) or (ze_high < 0):
            raise Exception(" Ze out of range")
        # Noise
        noise_index = self.__axis_dict__["Noise"]
        noise_low = -1
        noise_high = -1
        for low, high in zip(noise_index[:-1], noise_index[1:]):
            if (noise >= low) and (noise < high):
                noise_low = low
                noise_high = high
                break
        if (noise_low < 0) or (noise_high < 0):
            raise Exception("Noise out of range")
        # Done finding index values use for interpolation

        irf_data = []
        offset_index = self.__axis_dict__["AbsoluteOffset"]
        for az_i, az_val in [(0, az_low), (1, az_high)]:
            for ze_i, ze_val in [(0, ze_low), (1, ze_high)]:
                for noise_i, noise_val in [(0, noise_low), (1, noise_high)]:
                    irf_dict = {"Index": (az_i, ze_i, noise_i)}
                    az_shifted = az_val if az_val <= 360 else az_val - 360.0
                    ea_dict, ebias_dict, abias_dict = get_irf_not_safe(
                        self.__manager__,
                        offset_index,
                        az_shifted,
                        ze_val,
                        noise_val,
                        pointlike = self.__pointlike__,
                        psf_king = self.__psf_king__
                    )
                    irf_dict["EA_Dict"] = ea_dict
                    irf_dict["EBias_Dict"] = ebias_dict
                    irf_dict["ABias_Dict"] = abias_dict
                    irf_data.append(irf_dict)
                    # Load values
        # Initilaize Data container
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
        dtype=[
                ("ENERG_LO", ">f4", np.shape(elow)),
                ("ENERG_HI", ">f4", np.shape(ehigh)),
                ("THETA_LO", ">f4", np.shape(thetalow)),
                ("THETA_HI", ">f4", np.shape(thetahigh)),
                ("EFFAREA", ">f4", np.shape(ea_interpolated)),
            ]
        # Add event class data if applicable
        if self.__multi_eclass__:
            evtclass = np.full(
                shape=np.shape(ea_interpolated),fill_value=self.__event_class_idx__,dtype=np.uint16)
            dt='i2' #int16
            eclass_tuple = ('EVENT_CLASS', dt, np.shape(ea_interpolated))
            dtype.append(eclass_tuple)
            ea_final_data = np.array(
                [(elow, ehigh, thetalow, thetahigh, ea_interpolated, evtclass)],
                dtype)
        else:
            ea_final_data = np.array(
                [(elow, ehigh, thetalow, thetahigh, ea_interpolated)],
                dtype)

        # EBias
        elow = irf_data[0]["EBias_Dict"]["ELow"]
        ehigh = irf_data[0]["EBias_Dict"]["EHigh"]
        thetalow = irf_data[0]["EBias_Dict"]["ThetaLow"]
        thetahigh = irf_data[0]["EBias_Dict"]["ThetaHigh"]
        miglow = irf_data[0]["EBias_Dict"]["MigrationLow"]
        mighigh = irf_data[0]["EBias_Dict"]["MigrationHigh"]
        ebias_interpolated = ebias_interpolator((az, ze, noise))
        dtype=[
                ("ENERG_LO", ">f4", np.shape(elow)),
                ("ENERG_HI", ">f4", np.shape(ehigh)),
                ("MIGRA_LO", ">f4", np.shape(miglow)),
                ("MIGRA_HI", ">f4", np.shape(mighigh)),
                ("THETA_LO", ">f4", np.shape(thetalow)),
                ("THETA_HI", ">f4", np.shape(thetahigh)),
                ("MATRIX", ">f4", np.shape(ebias_interpolated))
        ]
        # Add event class data if applicable
        if self.__multi_eclass__:
            evtclass = np.full(
                shape=np.shape(ebias_interpolated), fill_value=self.__event_class_idx__, dtype=np.uint16)
            dt='i2' #int16  
            ebias_eclass_tuple = ('EVENT_CLASS', dt, np.shape(ebias_interpolated))
            dtype.append(ebias_eclass_tuple)
            ebias_final_data = np.array(
                [(elow, ehigh, miglow, mighigh, thetalow, thetahigh, ebias_interpolated, evtclass)],
                                          dtype)
        else:
            ebias_final_data = np.array(
                [(elow, ehigh, miglow, mighigh, thetalow, thetahigh, ebias_interpolated)],
                dtype)

        # ABias
        abias_final_data = None
        # If using PSF king. Not interpolated at the moment.
        if self.__psf_king__:
            # Already used our king azimuth values earlier
            if az - az_low > az_high - az:
                az_psf = az_high
            else:
                az_psf = az_low
            # Ze
            ze_index = self.__psf_king_index__["Zenith"]
            ze_low = -1
            ze_high = -1
            for low, high in zip(ze_index[:-1], ze_index[1:]):
                if (ze >= low) and (ze < high):
                    ze_low = low
                    ze_high = high
                    break
            if ze - ze_low > ze_high - ze:
                zen_psf = ze_high
            else:
                zen_psf = ze_low
            if (ze_low < 0) or (ze_high < 0):
                raise Exception("king ze out of range")
            # Noise
            noise_index = self.__psf_king_index__["Noise"]
            noise_low = -1
            noise_high = -1
            for low, high in zip(noise_index[:-1], noise_index[1:]):
                if (noise >= low) and (noise < high):
                    noise_low = low
                    noise_high = high
                    break
            if (noise_low < 0) or (noise_high < 0):
                raise Exception("king noise out of range")
            if noise - noise_low > noise_high - noise:
                noise_psf = noise_high
            else:
                noise_psf = noise_low

            abias_king_dict = get_psf_king(
                self.__psf_king_filename__, self.__psf_king_index__["AbsoluteOffset"], zen_psf, az_psf, noise_psf, msw_range)

            elow = abias_king_dict['ELow']
            ehigh = abias_king_dict['EHigh']
            thetalow = abias_king_dict['ThetaLow']
            thetahigh = abias_king_dict['ThetaHigh']
            psf_gamma = abias_king_dict['Gamma']
            psf_sigma = abias_king_dict['Sigma']

            dtype = [
                ('ENERG_LO', '>f4', np.shape(elow)),
                ('ENERG_HI', '>f4', np.shape(ehigh)),
                ('THETA_LO', '>f4', np.shape(thetalow)),
                ('THETA_HI', '>f4', np.shape(thetahigh)),
                ('GAMMA', '>f4', np.shape(psf_gamma)),
                ('SIGMA', '>f4', np.shape(psf_sigma))
            ]

            if self.__multi_eclass__:
                evtclass = np.full(shape=np.shape(psf_gamma), fill_value=self.__event_class_idx__, dtype=np.uint16)
                dt = 'i2'  # int16
                abias_eclass_tuple = ('EVENT_CLASS', dt, np.shape(psf_gamma))
                dtype.append(abias_eclass_tuple)
                abias_final_data = np.array([(elow, ehigh, thetalow, thetahigh, psf_gamma, psf_sigma, evtclass)],
                                            dtype)
            else:
                abias_final_data = np.array([(elow, ehigh, thetalow, thetahigh, psf_gamma, psf_sigma)],
                                            dtype)

        # else PSF table
        elif not self.__pointlike__:
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
