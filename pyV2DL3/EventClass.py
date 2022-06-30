import ROOT
import logging

from pyV2DL3.vegas.util import getCuts
logger = logging.getLogger(__name__)

"""
Construct an event class from an effective area cuts info

To extend event classes to include your parameters, start by adding them
here in order to access them when filling events and IRFs.

https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#Event_Classes
"""

class EventClass(object):
    def __init__(self, effective_area, irfs_to_store, override_cuts_validation=False):
        self.effective_area_IO = ROOT.VARootIO(effective_area, True)
        self.__override_cuts_validation__ = override_cuts_validation

        # Initialize the parameter names to search for
        self.cut_searches = ["ThetaSquareUpper",
                             "MeanScaledWidthLower",
                             "MeanScaledWidthUpper",
                             "MaxHeightLower",
                             "MaxHeightUpper",
                             "FoVCutUpper",
                             "FoVCutLower",
                             ]
        # Initialize corresponding class variables
        self.theta_square_upper = None
        self.msw_lower = None
        self.msw_upper = None
        self.max_height_lower = None
        self.max_height_upper = None
        self.fov_cut_lower = None
        self.fov_cut_upper = None

        # This dict will only contain keys from the found cuts.
        self.effective_area_IO.loadTheRootFile()
        for cuts in self.effective_area_IO.loadTheCutsInfo():
            ea_cut_dict = getCuts(cuts.fCutsFileText, self.cut_searches)

    # ------------- Handle params found in this EA's cuts -------------
        exception_messages = []
        # Theta^2 cut is currently required
        if "ThetaSquareUpper" in ea_cut_dict:
            self.theta_square_upper = float(ea_cut_dict["ThetaSquareUpper"])
            if irfs_to_store["point-like"] and self.theta_square_upper > 1:
                exception_messages.append(
                    "ThetaSquareUpper" + str(self.theta_square_upper) + " must be <= 1 for point-like")
            elif irfs_to_store["full-enclosure"] and self.theta_square_upper < 128881:
                exception_messages.append("ThetaSquareUpper: " + str(
                    self.theta_square_upper) + " must be >= 128881 (359 degrees) for full-enclosure")
        else:
            exception_messages.append("ThetaSquareUpper not found")

        # MSW cuts are optional
        if "MeanScaledWidthLower" in ea_cut_dict:
            self.msw_lower = float(ea_cut_dict["MeanScaledWidthLower"])
        else:
            self.msw_lower = float('-inf')
        if "MeanScaledWidthUpper" in ea_cut_dict:
            self.msw_upper = float(ea_cut_dict["MeanScaledWidthUpper"])
        else:
            self.msw_upper = float('inf')
        if (self.msw_lower >= self.msw_upper):
            exception_messages.append("MeanScaledWidthLower: " + str(
                self.msw_lower) + " must be < MeanScaledWidthUpper: " + str(self.msw_upper))

        # MaxHeightLower/Upper are optional, but must match event cuts if present
        if "MaxHeightLower" in ea_cut_dict:
            self.max_height_lower = float(ea_cut_dict["MaxHeightLower"])
        if "MaxHeightUpper" in ea_cut_dict:
            self.max_height_upper = float(ea_cut_dict["MaxHeightUpper"])

        # FoVCuts should match events if present
        if "FoVCutLower" in ea_cut_dict:
            self.fov_cut_lower = float(ea_cut_dict["FoVCutLower"])
        if "FoVCutUpper" in ea_cut_dict:
            self.fov_cut_upper = float(ea_cut_dict["FoVCutUpper"])

        # Raise exception with all exception messages for the user.
        if len(exception_messages) > 0:
            exception_string = ""
            for error in exception_messages:
                exception_string += error + '\n'
            if self.__override_cuts_validation__:
                logger.info(exception_messages)
            else:
                raise Exception(exception_string
                                + "The above problems were found with your EA cuts\n"
                                + "If you are absolutely certain that they are not an issue, "
                                + "you may use -o to override cut validation")


    """
    Called on VDS construction to ensure that events do not cut more deeply than this EA
    """

    def validate_vegas_cuts(self, vegas_file_IO):
        exception_messages = []
        vegas_cuts = vegas_file_IO.loadTheCutsInfo()
        for cuts in vegas_cuts:
            event_cuts = getCuts(cuts.fCutsFileText, self.cut_searches)

        # Events T^2 must match EA T^2
        if "ThetaSquareUpper" in event_cuts:
            if float(event_cuts["ThetaSquareUpper"]) != self.theta_square_upper:
                exception_messages.append(
                    "events ThetaSquareUpper:" + event_cuts["ThetaSquareUpper"] + " != EA ThetaSquareUpper: " + str(self.theta_square_upper))
        else:
            exception_messages.append(
                "ThetaSquareUpper not found in event cuts")

        # Event MSW cuts must not be deeper than EA MSW cuts
        if "MeanScaledWidthLower" in event_cuts:
            if float(event_cuts["MeanScaledWidthLower"]) > self.msw_lower:
                exception_messages.append(
                    "events MeanScaledWidthLower: " + event_cuts["MeanScaledWidthLower"] + " > EA MeanScaledWidthLower: " + str(self.msw_lower))
        if "MeanScaledWidthUpper" in event_cuts:
            if float(event_cuts["MeanScaledWidthUpper"]) < self.msw_upper:
                exception_messages.append(
                    "events MeanScaledWidthUpper: " + event_cuts["MeanScaledWidthUpper"] + " > EA MeanScaledWidthUpper: " + str(self.msw_upper))

        # Events MaxHeightLower/Upper are optional but must match EA
        if "MaxHeightLower" in event_cuts:
            if float(event_cuts["MaxHeightLower"]) != self.max_height_lower:
                exception_messages.append(
                    "events MaxHeightLower: " + event_cuts["MaxHeightLower"] + " != EA MaxHeightLower: " + str(self.max_height_lower))
        elif self.max_height_lower is not None:
            exception_messages.append(
                "MaxHeightLower was found in EA but not event cuts")
        if "MaxHeightUpper" in event_cuts:
            if float(event_cuts["MaxHeightUpper"]) != self.max_height_upper:
                exception_messages.append(
                    "events MaxHeightUpper: " + event_cuts["MaxHeightUpper"] + " != EA MaxHeightUpper: " + str(self.max_height_upper))
        elif self.max_height_upper is not None:
            exception_messages.append(
                "MaxHeightUpper was found in EA but not event cuts")

        # Events FoV cuts are optional but must match EA -- update to fovcut
        # if "FoVCutLower" in event_cuts:
        #     if float(event_cuts["FoVCutLower"]) != self.fov_cut_lower:
        #         exception_messages.append(
        #             "events FoVCutLower: " + event_cuts["FoVCutLower"] + " != EA FoVCutLower: " + str(self.fov_cut_lower))
        # elif self.fov_cut_lower is not None:
        #     exception_messages.append(
        #         "FoVCutLower was found in EA but not event cuts")
        # if "FoVCutUpper" in event_cuts:
        #     if float(event_cuts["FoVCutUpper"]) != self.fov_cut_upper:
        #         exception_messages.append(
        #             "events FoVCutUpper: " + event_cuts["FoVCutUpper"] + " != EA FoVCutUpper: " + str(self.fov_cut_upper))
        # elif self.fov_cut_upper is not None:
        #     exception_messages.append(
        #         "FoVCutUpper was found in EA but not event cuts")

        # Raise exception with all exception messages for the user.
        if len(exception_messages) > 0:
            exception_string = ""
            for error in exception_messages:
                exception_string += error + '\n'
            if self.__override_cuts_validation__:
                logger.info(exception_messages)
            else:
                raise Exception(exception_string
                                + "The above problems were found when validating your event cuts "
                                + "against your EA cuts\n"
                                + "If you are absolutely certain that they are not an issue, "
                                + "you may use -o to override cut validation")

    def __del__(self):
        cpy_nonestring = "<class 'CPyCppyy_NoneType'>"
        if self.effective_area_IO is not None:
            if str(type(self.effective_area_IO)) != cpy_nonestring:
                self.effective_area_IO.closeTheRootFile()
