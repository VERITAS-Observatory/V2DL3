import logging
from ctypes import c_float

import numpy as np
import ROOT

from pyV2DL3.vegas.load_vegas import VEGASStatus
from pyV2DL3.vegas.util import getCuts

logger = logging.getLogger(__name__)

"""
Construct an event class from an effective area file

Event Classes are wrappers for VEGAS effective area files to efficiently
read, store, and validate parameters for event cutting or sorting

To extend event classes to include more parameters, add them
here in order to access them when filling events and IRFs.

https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#Event_Classes
"""


class EffectiveAreaFile(object):
    def __init__(self, effective_area_file):
        self.__vegas__ = VEGASStatus()
        self.__vegas__.loadVEGAS()
        self.effective_area_IO = ROOT.VARootIO(effective_area_file, True)
        self.effective_area_IO.loadTheRootFile()
        self.manager = ROOT.VAEffectiveAreaManager()
        self.manager.setUseReconstructedEnergy(False)
        self.manager.loadEffectiveAreas(self.effective_area_IO)

        # Initialize the cuts parameter names to search for
        cut_searches = [
            "ThetaSquareUpper",
            "MeanScaledWidthLower",
            "MeanScaledWidthUpper",
            "MaxHeightLower",
            "MaxHeightUpper",
            "FoVCutUpper",
            "FoVCutLower",
        ]

        # Initialize corresponding class variables and default values
        self.theta_square_upper = None
        self.msw_lower = float("-inf")
        self.msw_upper = float("inf")
        self.max_height_lower = None
        self.max_height_upper = None
        self.fov_cut_lower = float("-inf")
        self.fov_cut_upper = float("inf")

        # Now load the cuts params
        self.__load_cuts_info__(cut_searches)

        # Build indexes for IRFs
        self.axis_dict, self.index_dict = self.__build_index__()

    """
    Build az, zen, noise, and offset indexes for this file.
    """

    def __build_index__(self):
        manager = self.manager
        axis = ["Azimuth", "Zenith", "Noise"]
        if len(manager.fEffectiveAreas) <= 0:
            raise Exception("No effective areas! ")
        index_check = manager.fEffectiveAreas.at(0).fDimensionNames
        for k in axis:
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
                index_dict["AbsoluteOffset"].append(0)
                index_dict["AbsoluteOffset"].append(10)

        # Validate Completeness
        axis_dict = {}
        check_num = 1

        for k in axis + ["AbsoluteOffset"]:
            check_num *= len(np.unique(index_dict[k]))
            axis_dict[k] = np.sort(np.unique(index_dict[k]))
            if len(axis_dict[k]) < 2 and k != "AbsoluteOffset":
                raise Exception("{} Axis need to have more than two values".format(k))
        return axis_dict, index_dict

    def get_safe_energy(self, az, ze, noise, offset=0.5, st6_configs=None):
        manager = self.manager
        effectiveAreaParameters = ROOT.VAEASimpleParameterData()
        effectiveAreaParameters.fAzimuth = az
        effectiveAreaParameters.fZenith = ze
        effectiveAreaParameters.fNoise = noise
        effectiveAreaParameters.fOffset = offset
        effectiveAreaParameters = manager.getVectorParamsFromSimpleParameterData(
            effectiveAreaParameters
        )
        minEnergy, maxEnergy = c_float(), c_float()
        if st6_configs is not None:
            split_configs = {
                opt.split(" ")[0]: opt.split(" ")[1] for opt in st6_configs if st6_configs is not None
            }
            if "EA_SafeEnergyRangeMethod" in split_configs.keys():
                safe_energy_method = str(split_configs["EA_SafeEnergyRangeMethod"])
                ea_uncertainty = float(split_configs["EA_MaxEffectiveAreaUncertainty"])
                energy_bias = float(split_configs["EA_MaxAllowedEnergyBias"])
                logger.debug(
                    f"Loaded st6 options EA_SafeEnergyRangeMethod: {safe_energy_method}, "
                    f"EA_MaxEffectiveAreaUncertainty: {ea_uncertainty}, "
                    f"EA_MaxAllowedEnergyBias: {energy_bias}"
                )
                self.manager.setOption("EA_SafeEnergyRangeMethod", safe_energy_method)
                self.manager.setOption("EA_MaxEffectiveAreaUncertainty", ea_uncertainty)
                self.manager.setOption("EA_MaxAllowedEnergyBias", energy_bias)

        manager.getSafeEnergyRange(effectiveAreaParameters, 0.5, minEnergy, maxEnergy)
        return minEnergy.value / 1000.0, maxEnergy.value / 1000.0

    """
    Loads and stores the effective area file's cuts parameters values
    """

    def __load_cuts_info__(self, cut_searches):
        # This dict will only contain keys from the found cuts.
        for cuts in self.effective_area_IO.loadTheCutsInfo():
            ea_cut_dict = getCuts(cuts.fCutsFileText, cut_searches)

        # FoVCuts are optional
        if "FoVCutLower" in ea_cut_dict:
            self.fov_cut_lower = float(ea_cut_dict["FoVCutLower"])
        if "FoVCutUpper" in ea_cut_dict:
            self.fov_cut_upper = float(ea_cut_dict["FoVCutUpper"])

        # MSW cuts are optional
        if "MeanScaledWidthLower" in ea_cut_dict:
            self.msw_lower = float(ea_cut_dict["MeanScaledWidthLower"])
        if "MeanScaledWidthUpper" in ea_cut_dict:
            self.msw_upper = float(ea_cut_dict["MeanScaledWidthUpper"])

        if self.msw_lower >= self.msw_upper:
            raise Exception(
                "MeanScaledWidthLower: "
                + str(self.msw_lower)
                + " must be < MeanScaledWidthUpper: "
                + str(self.msw_upper)
            )

        # Theta^2 cut is required
        if "ThetaSquareUpper" in ea_cut_dict:
            self.theta_square_upper = float(ea_cut_dict["ThetaSquareUpper"])
        else:
            raise Exception("ThetaSquareUpper not found in EA cuts parameters")
