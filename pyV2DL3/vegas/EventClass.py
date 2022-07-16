import ROOT
import logging

from pyV2DL3.vegas.util import getCuts
logger = logging.getLogger(__name__)

"""
Construct an event class from an effective area cuts info

Event Classes are wrappers for VEGAS effective area files to efficiently
read, store, and validate the parameters needed for cuts or sorting

To extend event classes to include more parameters, add them
here in order to access them when filling events and IRFs.

https://veritas.sao.arizona.edu/wiki/V2dl3_dev_notes#Event_Classes
"""


class EventClass(object):
    def __init__(self, effective_area):
        self.effective_area_IO = ROOT.VARootIO(effective_area, True)

        # Initialize the cuts parameter names to search for
        cut_searches = [
                        "ThetaSquareUpper",
                        "MeanScaledWidthLower", "MeanScaledWidthUpper",  # MSW
                        "MaxHeightLower", "MaxHeightUpper",              # Max height
                        "FoVCutUpper", "FoVCutLower",                    # Field of view
                        ]
                        
        # Initialize corresponding class variables
        self.theta_square_upper = None
        self.msw_lower = None
        self.msw_upper = None
        self.max_height_lower = None
        self.max_height_upper = None
        self.fov_cut_lower = None
        self.fov_cut_upper = None

        # Now load the cuts params
        self.__load_cuts_info__(cut_searches)

        
    def __del__(self):
        cpy_nonestring = "<class 'CPyCppyy_NoneType'>"
        if self.effective_area_IO is not None:
            if str(type(self.effective_area_IO)) != cpy_nonestring:
                self.effective_area_IO.closeTheRootFile()


    """
    Loads and stores the effective area's cuts parameters values
    """
    def __load_cuts_info__(self, cut_searches):
        # This dict will only contain keys from the found cuts.
        self.effective_area_IO.loadTheRootFile()
        for cuts in self.effective_area_IO.loadTheCutsInfo():
            ea_cut_dict = getCuts(cuts.fCutsFileText, cut_searches)

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