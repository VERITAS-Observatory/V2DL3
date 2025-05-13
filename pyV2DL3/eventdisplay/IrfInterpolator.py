import logging
import os.path

import click
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from sklearn.neighbors import KNeighborsRegressor

from pyV2DL3.eventdisplay.IrfExtractor import extract_irf, extract_irf_for_knn
from pyV2DL3.eventdisplay.util import WrongIrf, duplicate_dimensions


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
        self.interpolator = None
        self.interpolator_name = "KNeighborsRegressor"

        if os.path.isfile(filename):
            self.filename = filename
        else:
            logging.error(f"IRF file not found {filename}")
            raise FileNotFoundError

    def set_irf(self, irf_name, **kwargs):
        """Check consistency of IRF name"""
        if (
            irf_name in self.implemented_irf_names_1d
            or irf_name in self.implemented_irf_names_2d
        ):
            self.irf_name = irf_name
            self._load_irf(**kwargs)
        else:
            logging.error(
                "The irf you entered: {} is either wrong or not implemented.".format(
                    irf_name
                )
            )
            raise WrongIrf

    def _load_irf(self, **kwargs):
        """Load IRFs from effective area file"""

        logging.info(
            "Extracting IRFs of type: {0} for azimuth {1} deg".format(
                self.irf_name, np.array2string(self.azimuth, precision=2)
            )
        )

        if self.interpolator_name == "KNeighborsRegressor":
            self._load_irf_for_knn(**kwargs)
        else:
            self._load_irf_for_regular_grid_interpolator(**kwargs)

    def _load_irf_for_knn(self, **kwargs):
        """Load IRFs from file for KNeighborsRegressor"""

        coords, values = extract_irf_for_knn(
            self.filename,
            self.irf_name,
            irf1d=self.irf_name in self.implemented_irf_names_1d,
            azimuth=self.azimuth,
        )
        self.interpolator = KNeighborsRegressor(n_neighbors=5, weights="distance")
        self.interpolator.fit(coords, values)

        # Set irf_axes based on the coordinates
        if self.irf_name in self.implemented_irf_names_1d:
            self.irf_axes = [np.unique(coords[:, 3])]  # energy axis
        else:
            self.irf_axes = [
                np.unique(coords[:, 3]),  # x dimension
                np.unique(coords[:, 4])   # y dimension
            ]
        logging.debug(str(("IRF axes:", self.irf_axes)))

    def _load_irf_for_regular_grid_interpolator(self, **kwargs):
        """Load IRFs from file for RegularGridInterpolator"""

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
        zenith_axis = None
        for i, axis in enumerate(irf_axes):
            if len(irf_axes[axis]) == 1:
                irf_axes[axis] = np.concatenate(
                    (irf_axes[axis].flatten(), irf_axes[axis].flatten() + 0.01), axis=None
                )

            if axis == 'zeniths':
                irf_axes['zeniths'] = np.cos(np.radians(irf_axes['zeniths']))[::-1]
                zenith_axis = i
                logging.debug("zenith axis index: {}".format(zenith_axis))

        if zenith_axis is None:
            logging.error("zenith axis not found in irf_axes")
            raise ValueError

        self.irf_data = np.flip(irf_data, axis=zenith_axis)
        self.irf_axes = list(irf_axes.values())
        logging.debug(str(("IRF axes:", irf_axes)))

        if kwargs.get("use_click", True):
            clk = click.get_current_context()
            extrapolation = clk.params["force_extrapolation"]
        else:
            extrapolation = kwargs.get("force_extrapolation", False)

        if extrapolation:
            self.interpolator = RegularGridInterpolator(
                self.irf_axes, self.irf_data, bounds_error=False, fill_value=None)
        else:
            self.interpolator = RegularGridInterpolator(self.irf_axes, self.irf_data)

    def interpolate(self, coordinate):
        coordinate[1] = np.cos(np.radians(coordinate[1]))
        for c in coordinate:
            logging.debug("Interpolating coordinates: {0:.2f}".format(c))

        # The interpolation is slightly different for 1D or 2D IRFs.
        if self.azimuth == 0:
            if len(coordinate) != 4:
                logging.error(
                    "IRF interpolation: for azimuth 0, require 4 coordinates "
                    "(azimuth,  pedvar, zenith, offset)"
                )
                raise ValueError
        elif len(coordinate) != 3:
            logging.error(
                "IRF Interpolation: Require 3 coordinates (pedvar, zenith, offset)"
            )
            raise ValueError

        if self.irf_name in self.implemented_irf_names_2d:
            return self._interpolate_2d(coordinate, self.irf_axes)
        elif self.irf_name in self.implemented_irf_names_1d:
            return self._interpolate_1d(coordinate, self.irf_axes[0])
        else:
            logging.error("The requested %s" " is not available.", self.irf_name)
            raise WrongIrf

    def _interpolate_2d(self, coordinate, irf_axes):
        """Interpolate IRF for 2D IRFs."""
        xx, yy = np.meshgrid(self.irf_axes[0], self.irf_axes[1])
        xx_flat = xx.flatten()
        yy_flat = yy.flatten()

        try:
            if self.interpolator_name == "KNeighborsRegressor":
                predict_coords = np.array([
                    [coordinate[0], coordinate[1], coordinate[2], x, y]
                    for x, y in zip(xx_flat, yy_flat)
                ])
                interpolated_irf = self.interpolator.predict(predict_coords)
                interpolated_irf = interpolated_irf.reshape(xx.shape)
            else:
                interpolated_irf = self.interpolator((xx, yy, *coordinate))
        except ValueError:
            raise ValueError("IRF interpolation failed for axis %s", self.irf_name)
        return interpolated_irf, [irf_axes[0], irf_axes[1]]

    def _interpolate_1d(self, coordinate, irf_axis):
        """Interpolate IRF for 1D IRFs (energy axis only)."""
        try:
            if self.interpolator_name == "KNeighborsRegressor":
                interpolated_irf = self.interpolator.predict(
                    np.array([[coordinate[0], coordinate[1], coordinate[2], e] for e in irf_axis])
                )
            else:
                interpolated_irf = self.interpolator((irf_axis, *coordinate))
        except ValueError:
            raise ValueError("IRF interpolation failed for axis %s", self.irf_name)

        return interpolated_irf, [irf_axis]
