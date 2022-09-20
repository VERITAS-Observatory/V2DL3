import ROOT
import logging

from pyV2DL3.vegas.util import getCuts
logger = logging.getLogger(__name__)

"""
Construct an event class from an effective area cuts info.

To extend event classes to include your parameters, start by adding them
here in order to access them when filling events and IRFs.

https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#Event_Classes
"""


class EventClass(object):
    def __init__(self, effective_area):
        self.effective_area_IO = ROOT.VARootIO(effective_area, True)

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

        # MSW cuts are optional
        if "MeanScaledWidthLower" in ea_cut_dict:
            self.msw_lower = float(ea_cut_dict["MeanScaledWidthLower"])
        else:
            # Assign +/- inf so that comparisons will still work when filling events
            self.msw_lower = float('-inf')
        if "MeanScaledWidthUpper" in ea_cut_dict:
            self.msw_upper = float(ea_cut_dict["MeanScaledWidthUpper"])
        else:
            self.msw_upper = float('inf')
        if (self.msw_lower >= self.msw_upper):
            raise Exception("MeanScaledWidthLower: " + str(
                self.msw_lower) + " must be < MeanScaledWidthUpper: " + str(self.msw_upper))

    def __del__(self):
        cpy_nonestring = "<class 'CPyCppyy_NoneType'>"
        if self.effective_area_IO is not None:
            if str(type(self.effective_area_IO)) != cpy_nonestring:
                self.effective_area_IO.closeTheRootFile()
