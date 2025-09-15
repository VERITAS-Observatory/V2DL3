import logging
import os.path

import click
import numpy as np
from scipy.interpolate import (
    LinearNDInterpolator,
    RBFInterpolator,
    RegularGridInterpolator,
)
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MinMaxScaler

from pyV2DL3.eventdisplay.IrfExtractor import extract_irf, extract_irf_for_knn
from pyV2DL3.eventdisplay.util import WrongIrf, duplicate_dimensions


class IrfInterpolator:
    def __init__(self, filename, azimuth, interpolator_name):
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
        self.interpolator_name = interpolator_name
        if interpolator_name not in ["KNeighborsRegressor", "RegularGridInterpolator",
                                     "LinearNDInterpolator", "RBFInterpolator"]:
            raise ValueError(
                "The interpolator you entered: {} is either wrong or not implemented.".format(
                    interpolator_name
                )
            )

        if os.path.isfile(filename):
            self.filename = filename
        else:
            logging.error(f"IRF file not found {filename}")
            raise FileNotFoundError

    def set_irf(self, irf_name, **kwargs):
        """Check consistency of IRF name

        Parameters:
        -----------
        irf_name : str
            Name of the IRF to load
        fill_empty_bins : bool, optional
            Whether to fill empty bins in RegularGridInterpolator (default: True)
        """
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
            self._load_irf_for_knn()
        elif self.interpolator_name == "LinearNDInterpolator":
            self._load_irf_for_linear_nd_interpolator()
        elif self.interpolator_name == "RBFInterpolator":
            self._load_irf_for_rbf_interpolator()
        else:
            self._load_irf_for_regular_grid_interpolator(**kwargs)

    def _load_irf_for_knn(self):
        """Load IRFs from file for KNeighborsRegressor"""

        coords, values = extract_irf_for_knn(
            self.filename,
            self.irf_name,
            irf1d=self.irf_name in self.implemented_irf_names_1d,
            azimuth=self.azimuth,
        )
        self.interpolator = make_pipeline(
            MinMaxScaler(), KNeighborsRegressor(n_neighbors=5, weights='distance')
        )
        self.interpolator.fit(coords, values)

        if self.irf_name in self.implemented_irf_names_1d:
            self.irf_axes = [np.unique(coords[:, 3])]  # energy axis
        else:
            self.irf_axes = [
                np.unique(coords[:, 3]),  # x dimension
                np.unique(coords[:, 4])   # y dimension
            ]
        logging.debug(str(("IRF axes:", self.irf_axes)))

    def _load_irf_for_rbf_interpolator(self):
        """Load IRFs from file for RBFInterpolator with performance optimizations"""

        coords, values = extract_irf_for_knn(
            self.filename,
            self.irf_name,
            irf1d=self.irf_name in self.implemented_irf_names_1d,
            azimuth=self.azimuth,
        )

        # Store the original coordinates for later use
        self.coords = coords
        self.values = values

        # Scale coordinates to improve numerical stability and performance
        scaled_coords = coords
        # Store axes for later use
        if self.irf_name in self.implemented_irf_names_1d:
            self.irf_axes = [np.unique(coords[:, 3])]  # energy axis

            # For 1D IRFs, use multi-quadric kernel which provides better smoothness
            self.interpolator = RBFInterpolator(
                scaled_coords,
                values,
                neighbors=5,
                kernel='linear',       # Or 'thin_plate_spline'
                epsilon=3,
                degree=0,             # Linear polynomial helps global behavior
            )

        else:
            x_axis = np.unique(coords[:, 3])
            y_axis = np.unique(coords[:, 4])
            self.irf_axes = [x_axis, y_axis]
            xx, yy = np.meshgrid(x_axis, y_axis, indexing='xy')
            self.mesh_points = np.array([[
                0, 0, 0, x, y] for x, y in zip(xx.flatten(), yy.flatten())])

            # For 2D IRFs, use more robust settings to avoid singular matrix
            self.interpolator = RBFInterpolator(
                scaled_coords,
                values,
                neighbors=8,
                kernel='cubic',    # More stable than thin_plate_spline
                epsilon=3.0,       # Higher epsilon for better conditioning
                degree=-1,         # No polynomial term, more stable
            )

        logging.debug(str(("IRF axes:", self.irf_axes)))

    def _load_irf_for_linear_nd_interpolator(self):
        """Load IRFs from file for LinearNDInterpolator with performance optimizations"""

        coords, values = extract_irf_for_knn(
            self.filename,
            self.irf_name,
            irf1d=self.irf_name in self.implemented_irf_names_1d,
            azimuth=self.azimuth,
        )

        # Store the original coordinates for later use
        self.coords = coords
        self.values = values

        # OPTIMIZATION 1: Use rescaled coordinates
        # Scale coordinates to improve numerical stability and performance
        from sklearn.preprocessing import MinMaxScaler
        self.scaler = MinMaxScaler()
        scaled_coords = self.scaler.fit_transform(coords)

        # OPTIMIZATION 2: Cache meshgrid
        if self.irf_name in self.implemented_irf_names_1d:
            self.irf_axes = [np.unique(coords[:, 3])]  # energy axis
        else:
            x_axis = np.unique(coords[:, 3])
            y_axis = np.unique(coords[:, 4])
            self.irf_axes = [x_axis, y_axis]
            # Pre-compute the meshgrid for faster interpolation
            xx, yy = np.meshgrid(x_axis, y_axis, indexing='xy')
            self.mesh_points = np.array([[
                0, 0, 0, x, y] for x, y in zip(xx.flatten(), yy.flatten())])

        # Create the interpolator with optimized parameters
        self.interpolator = LinearNDInterpolator(scaled_coords, values, fill_value=0)

        logging.debug(str(("IRF axes:", self.irf_axes)))

    def _load_irf_for_regular_grid_interpolator(self, **kwargs):
        """Load IRFs from file for RegularGridInterpolator with optional bin filling"""

        irf_data, irf_axes = extract_irf(
            self.filename,
            self.irf_name,
            irf1d=(self.irf_name in self.implemented_irf_names_1d),
            azimuth=self.azimuth,
        )

        # Check if we should artificially zero out a slice for testing
        zero_out_slice = kwargs.get("zero_out_slice", None)
        if zero_out_slice is not None:
            # Determine offset dimension based on IRF type and dimensionality
            if self.irf_name in ["eff", "Rec_eff", "effNoTh2", "Rec_effNoTh2"]:
                # For effective area, offset is always the last dimension
                offset_dim = irf_data.ndim - 1
                logging.info(f"Effective area IRF: Using dimension {offset_dim} as offset dimension")
            else:
                # For other IRF types, follow standard convention
                offset_dim = 2 if irf_data.ndim > 2 else 1
                logging.info(f"Standard IRF: Using dimension {offset_dim} as offset dimension")

            if 0 <= zero_out_slice < irf_data.shape[offset_dim]:
                # Create slice indices
                idx = [slice(None)] * irf_data.ndim
                idx[offset_dim] = zero_out_slice

                # Store the original values for logging
                before_zero = irf_data[tuple(idx)].copy()
                nonzero_count = np.sum(before_zero > 0)
                nonzero_max = np.max(before_zero) if nonzero_count > 0 else 0

                # Zero out the slice
                irf_data[tuple(idx)] = 0

                # Verify the zero-out worked
                after_zero = irf_data[tuple(idx)]
                if np.sum(after_zero > 0) > 0:
                    logging.error(f"ERROR: Failed to zero out offset slice {zero_out_slice}")
                else:
                    logging.info(f"TEST MODE: Successfully zeroed out offset slice {zero_out_slice}")
                    logging.info(f"  Zeroed {nonzero_count} values (max was {nonzero_max:.2f})")

                # For debugging: show IRF shape
                logging.info(f"IRF data shape: {irf_data.shape}")

        # Flag to determine if bin filling should be applied
        fill_empty_bins = kwargs.get("fill_empty_bins", True)

        # Fill empty bins within the data hull before interpolation if enabled
        if fill_empty_bins:
            logging.info("Filling empty IRF bins before interpolation")
            logging.warning(
                "WARNING: Filling empty bins is not the standard method and introduces additional "
                "systematic uncertainties (approximately 10%) due to linear interpolation"
            )
            irf_data = self._fill_empty_bins(irf_data)
        else:
            logging.info("Using original IRF data without filling empty bins")

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
                irf_axes['zeniths'] = 1. / np.cos(np.radians(irf_axes['zeniths']))[::-1]
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

    def _fill_empty_bins(self, data):
        """Fill empty bins using an average energy profile as a template"""

        # Create a copy of the data for filling
        filled_data = data.copy()

        # =====================================================================
        # PHASE 0: CREATE AVERAGE ENERGY PROFILE TO DETERMINE VALID RANGES
        # =====================================================================
        print("Phase 0: Creating average energy profile to determine valid ranges")

        # Determine offset dimension based on IRF type and dimensionality - FIXED!
        if self.irf_name in ["eff", "Rec_eff", "effNoTh2", "Rec_effNoTh2"]:
            # For effective area, offset is always the last dimension
            offset_dim = data.ndim - 1
            print(f"  Effective area IRF: Using dimension {offset_dim} as offset dimension")
        else:
            # For other IRF types, follow standard convention
            offset_dim = 2 if data.ndim > 2 else 1
            print(f"  Standard IRF: Using dimension {offset_dim} as offset dimension")

        energy_dim = 0  # Energy is always the first dimension

        # Create boundary mask for regions that should stay as zeros
        boundary_mask = np.zeros_like(filled_data, dtype=bool)

        # Calculate threshold for identifying significant values
        data_max = np.max(data)
        significance_threshold = data_max * 0.01  # 1% of max value
        print(f"  Using threshold of {significance_threshold:.2f} for boundary detection")

        # Create a combined energy profile by averaging across all offset slices
        # First, create a mask of where values are significant
        # For each non-offset dimension, find the valid energy range

        # Collect valid energy ranges across all non-zeroed slices
        valid_energy_points = np.zeros(data.shape[energy_dim])
        valid_energy_counts = np.zeros(data.shape[energy_dim])

        # Loop through all offset slices
        for offset_pos in range(data.shape[offset_dim]):
            # Skip completely zeroed slices (like our test slice)
            offset_idx = [slice(None)] * data.ndim
            offset_idx[offset_dim] = offset_pos
            offset_slice = data[tuple(offset_idx)]

            if np.max(offset_slice) <= 0:
                print(f"  Skipping offset {offset_pos} (all zeros)")
                continue

            # Loop through all combinations of other dimensions (excluding energy and offset)
            other_dims = [d for d in range(1, data.ndim) if d != offset_dim]
            other_shapes = [data.shape[d] for d in other_dims]

            for idx in np.ndindex(*other_shapes):
                # Convert flat index to multi-dimensional index for the other dimensions
                full_idx = [slice(None)] * data.ndim  # Start with all slices
                full_idx[offset_dim] = offset_pos     # Set offset dimension

                # Set other dimensions (excluding energy and offset)
                dim_counter = 0
                for d in other_dims:
                    full_idx[d] = idx[dim_counter]
                    dim_counter += 1

                # Get energy profile for this specific combination
                energy_profile = data[tuple(full_idx)]

                # Find where values are significant
                significant_mask = energy_profile > significance_threshold

                # Add to our running total
                valid_energy_points += significant_mask.astype(float)
                valid_energy_counts += 1.0

        # Calculate the average profile (percentage of slices where each energy bin is valid)
        if np.sum(valid_energy_counts) > 0:
            avg_energy_profile = valid_energy_points / valid_energy_counts

            # Print the average profile for debugging
            print("  Average energy profile (% of slices where energy is valid):")
            for e in range(len(avg_energy_profile)):
                if e % 5 == 0:  # Print every 5th bin to avoid too much output
                    print(f"    Energy bin {e}: {avg_energy_profile[e]:.1%}")

            # Determine energy threshold - consider an energy bin valid if it's
            # significant in at least 20% of the valid slices
            valid_energy_threshold = 0.20
            energy_validity_mask = avg_energy_profile >= valid_energy_threshold

            # Find the energy range where bins are consistently valid
            valid_energy_indices = np.where(energy_validity_mask)[0]

            if len(valid_energy_indices) > 0:
                min_valid_energy = valid_energy_indices.min()
                max_valid_energy = valid_energy_indices.max()

                print(f"  Global valid energy range: {min_valid_energy} to {max_valid_energy}")

                # Apply this valid energy range to all slices
                # Mark energies outside this range as boundary (should stay zero)
                for idx in np.ndindex(*[data.shape[d] for d in range(1, data.ndim)]):
                    # Mark lower energy bins as boundary
                    for e in range(0, min_valid_energy):
                        boundary_idx = tuple([e] + list(idx))
                        boundary_mask[boundary_idx] = True

                    # Mark higher energy bins as boundary
                    for e in range(max_valid_energy + 1, data.shape[energy_dim]):
                        boundary_idx = tuple([e] + list(idx))
                        boundary_mask[boundary_idx] = True
            else:
                print("  WARNING: No valid energy range found")
        else:
            print("  WARNING: No valid slices found for energy profile")

        print(f"  Identified {np.sum(boundary_mask)} boundary zeros to preserve")

        # =====================================================================
        # PHASE 1: LINEAR SLICE INTERPOLATION (RESPECTING ENERGY BOUNDARIES)
        # =====================================================================
        print("Phase 1: Filling empty or partially filled slices (respecting energy boundaries)")

        # Skip if offset dimension is too small
        if data.shape[offset_dim] < 3:
            print(f"Offset dimension {offset_dim} has fewer than 3 points, skipping slice interpolation")
            return filled_data

        print(f"Processing offset dimension {offset_dim}")

        # Look at each interior offset slice
        for pos in range(1, data.shape[offset_dim]-1):
            idx = [slice(None)] * data.ndim
            idx[offset_dim] = pos
            slice_data = filled_data[tuple(idx)]

            # Get neighboring slices
            left_idx = [slice(None)] * data.ndim
            left_idx[offset_dim] = pos - 1
            left_slice = filled_data[tuple(left_idx)]

            right_idx = [slice(None)] * data.ndim
            right_idx[offset_dim] = pos + 1
            right_slice = filled_data[tuple(right_idx)]

            # Get energy boundary mask for this slice
            slice_boundary_mask = boundary_mask[tuple(idx)]

            # Only fill zeros that aren't part of the energy boundary
            fillable_zeros = (slice_data == 0) & ~slice_boundary_mask
            zero_ratio = np.sum(fillable_zeros) / slice_data.size

            # Check if this is likely an artificially zeroed slice (high zero ratio)
            is_artificial_zero = zero_ratio > 0.5  # Over 50% zeros suggests artificial zeroing

            if zero_ratio > 0:
                if is_artificial_zero:
                    print(f"  Offset slice {pos} has {zero_ratio:.1%} fillable zeros - likely artificially zeroed")
                else:
                    print(f"  Offset slice {pos} has {zero_ratio:.1%} fillable zeros")

                # Only fill zeros that have valid values on both sides
                min_valid_value = data_max * 0.001  # 0.1% of maximum as minimum valid value
                both_valid = (left_slice > min_valid_value) & (right_slice > min_valid_value) & fillable_zeros

                if np.any(both_valid):
                    # Simple linear interpolation
                    filled_data[tuple(idx)][both_valid] = 0.5 * left_slice[both_valid] + 0.5 * right_slice[both_valid]
                    print(f"    Interpolated {np.sum(both_valid)} zeros with both neighbors")

        return filled_data

    def interpolate(self, coordinate):
        coordinate[1] = 1. / np.cos(np.radians(coordinate[1]))
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
            logging.error(f"The requested {self.irf_name} is not available.")
            raise WrongIrf

    def _interpolate_2d(self, coordinate, irf_axes):
        """Interpolate IRF for 2D IRFs."""
        xx, yy = np.meshgrid(irf_axes[0], irf_axes[1], indexing='xy')
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
            elif self.interpolator_name in ["LinearNDInterpolator", "RBFInterpolator"]:
                # Scale the input coordinates if using RBF or LinearND
                if hasattr(self, 'scaler'):
                    # Build coordinate array for all points to predict
                    predict_coords = np.array([
                        [coordinate[0], coordinate[1], coordinate[2], x, y]
                        for x, y in zip(xx_flat, yy_flat)
                    ])
                    # Scale the coordinates
                    predict_coords = self.scaler.transform(predict_coords)
                    # Apply interpolation
                    interpolated_irf = self.interpolator(predict_coords)
                    interpolated_irf = interpolated_irf.reshape(xx.shape)
                else:
                    # Fallback if scaler is not available
                    predict_coords = np.array([
                        [coordinate[0], coordinate[1], coordinate[2], x, y]
                        for x, y in zip(xx_flat, yy_flat)
                    ])
                    interpolated_irf = self.interpolator(predict_coords)
                    interpolated_irf = interpolated_irf.reshape(xx.shape)
            else:
                interpolated_irf = self.interpolator((xx, yy, *coordinate))
        except ValueError:
            raise ValueError(f"IRF interpolation failed for axis {self.irf_name}")
        return interpolated_irf, [irf_axes[0], irf_axes[1]]

    def _interpolate_1d(self, coordinate, irf_axis):
        """Interpolate IRF for 1D IRFs (energy axis only)."""
        try:
            if self.interpolator_name == "KNeighborsRegressor":
                interpolated_irf = self.interpolator.predict(
                    np.array([[coordinate[0], coordinate[1], coordinate[2], e] for e in irf_axis])
                )
            elif self.interpolator_name in ["LinearNDInterpolator", "RBFInterpolator"]:
                # Scale the input coordinates if using RBF or LinearND
                if hasattr(self, 'scaler'):
                    predict_coords = np.array([[coordinate[0], coordinate[1], coordinate[2], e]
                                              for e in irf_axis])
                    predict_coords = self.scaler.transform(predict_coords)
                    interpolated_irf = self.interpolator(predict_coords)
                else:
                    predict_coords = np.array([[coordinate[0], coordinate[1], coordinate[2], e]
                                              for e in irf_axis])
                    interpolated_irf = self.interpolator(predict_coords)
            else:
                interpolated_irf = self.interpolator((irf_axis, *coordinate))
        except ValueError:
            raise ValueError(f"IRF interpolation failed for axis {self.irf_name}")

        return interpolated_irf, [irf_axis]
