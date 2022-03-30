import logging
import numpy as np
import os.path
from pyV2DL3.eventdisplay.IrfExtractor import extract_irf
from pyV2DL3.eventdisplay.util import duplicate_dimensions
from pyV2DL3.eventdisplay.util import WrongIrf
from scipy.interpolate import RegularGridInterpolator


class IrfInterpolator:
    def __init__(self, filename, azimuth):
        self.implemented_irf_names_1d = ["eff", "Rec_eff", "effNoTh2", "Rec_effNoTh2"]
        self.implemented_irf_names_2d = [
            "hEsysMCRelative2D",
            "hEsysMCRelative2DNoDirectionCut",
            "hAngularLogDiff_2D",
            "hAngularLogDiffEmc_2D",
        ]
        self.irf_name = ""
        self.azimuth = azimuth

        if os.path.isfile(filename):
            self.filename = filename
        else:
            raise FileNotFoundError

    def set_irf(self, irf_name):
        """Check consistency of IRF name"""
        if (
            irf_name in self.implemented_irf_names_1d
            or irf_name in self.implemented_irf_names_2d
        ):
            self.irf_name = irf_name
            self.__load_irf()
        else:
            logging.exception(
                "The irf you entered: {} is either wrong or not implemented.".format(
                    irf_name
                )
            )
            raise WrongIrf

    def __load_irf(self):
        """Load IRFs from effective area file"""

        logging.info(
            "Extracting IRFs of type: {0} for azimuth {1} deg".format(
                self.irf_name, np.array2string(self.azimuth, precision=2)
            )
        )
        irf_data, irf_axes = extract_irf(
            self.filename,
            self.irf_name,
            irf1d=(self.irf_name in self.implemented_irf_names_1d),
            azimuth=self.azimuth,
        )
        # This is an important technical step:
        # the regular grid interpolator does not accept
        # interpolating on a dimension with size = 1.
        # Make sure that there are no size 1 dimensions.
        # Do the same with the axes:
        irf_data = duplicate_dimensions(irf_data)
        # Also the coordinates of the axes need to be in increasing order.
        for i, axis in enumerate(irf_axes):
            if len(axis) == 1:
                irf_axes[i] = np.concatenate(
                    (axis.flatten(), axis.flatten() + 0.01), axis=None
                )
            print(i, irf_axes[i])
           
            
        irf_axes[-2] = 1. /np.cos(np.radians(irf_axes[-2]))
        print(irf_axes[-2])

        self.irf_data = irf_data
        self.irf_axes = irf_axes 
        print('IRF axes:', irf_axes)
        self.interpolator = RegularGridInterpolator(self.irf_axes, self.irf_data)

    def interpolate(self, coordinate):
        for c in coordinate:
            logging.debug("Interpolating coordinates: {0:.2f}".format(c))
            #print("cordinate: ",c)
        
        print('co:',coordinate)
        coordinate[1] = 1. / np.cos(np.radians(coordinate[1])) 
        print('co:',coordinate) 

        # The interpolation is slightly different for 1D or 2D IRFs.
        if self.azimuth == 0:
            if len(coordinate) != 4:
                raise ValueError(
                    "IRF interpolation: for azimuth 0, require 4 coordinates "
                    "(azimuth,  pedvar, zenith, offset)"
                )
        elif len(coordinate) != 3:
            raise ValueError(
                "IRF Interpolation: Require 3 coordinates (pedvar, zenith, offset)"
            )

        if self.irf_name in self.implemented_irf_names_2d:
            # In this case, the interpolator needs to interpolate over 2 dimensions:
            xx, yy = np.meshgrid(self.irf_axes[0], self.irf_axes[1])
            interpolated_irf = self.interpolator((xx, yy, *coordinate))
            return interpolated_irf, [self.irf_axes[0], self.irf_axes[1]]
        elif self.irf_name in self.implemented_irf_names_1d:
            # In this case, the interpolator needs to interpolate only
            # over 1 dimension (true energy):
            interpolated_irf = self.interpolator((self.irf_axes[0], *coordinate))
            return interpolated_irf, [self.irf_axes[0]]
        else:
            logging.exception(
                "The irf you entered: {}" " is not available.".format(self.irf_name)
            )
            raise WrongIrf
