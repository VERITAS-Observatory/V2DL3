import os

import numpy as np
import pandas as pd
import uproot
import xarray as xr


class XIRFExtractor:
    implemented_irf_names_1d = [
        "eff", "Rec_eff", "effNoTh2", "Rec_effNoTh2"
    ]
    implemented_irf_names_2d = [
        "hEsysMCRelative2D",
        "hEsysMCRelative2DNoDirectionCut",
        "hAngularLogDiff_2D",
        "hAngularLogDiffEmc_2D",
    ]
    
    @staticmethod
    def _convert_index_to_dims(irf_ds, **indCoords):
        # Time to turn indexed coordinates into dimensions
        # https://stackoverflow.com/a/70873363
        #
        names, coords = tuple(zip(*indCoords.items()))
        return irf_ds.assign_coords(
            {
                "index": pd.MultiIndex.from_arrays(
                    coords, names=names
                )
            }
        ).unstack("index").fillna(0)
    
    @staticmethod
    def _get_2d_irf_axis(in_min, in_max, in_bins):
        binwidth = abs(in_max - in_min) / in_bins / 2.0
        return np.linspace(in_max + binwidth, in_max - binwidth, in_bins)
    
    @staticmethod
    def _to_uniform_numpy(in_arr):
        return np.array(tuple(in_arr), dtype=in_arr[0].dtype)
    
    def _extract(self, *leafs):
        return self.tree_ptr.arrays(leafs, library="np").values()
    
    def _get_1d_irf(self, irf_name):
        az_min, az_max, pedvar, woff, ze, e0, az, oned_irf = self._extract(
            "azMin", "azMax", "pedvar", "Woff", "ze", "e0", "az", irf_name
        )
        az_centers = (az_max + az_min) / 2.0
        az_centers[az_centers < 0] += 360
        e0 = self._to_uniform_numpy(e0)
        oned_irf = self._to_uniform_numpy(oned_irf)
        irf_ds = xr.Dataset(
            {
                "irf": (["index", "energy"], oned_irf),
            },
            coords={
                "energy": e0[0],
            },
            attrs={
                "irf_filename": self.file_name,
                "irf_name": irf_name,  # Additional data
            },
        )
        
        # Time to turn indexed coordinates into dimensions
        # https://stackoverflow.com/a/70873363
        
        return self._convert_index_to_dims(
            irf_ds, az=az_centers, pedvar=pedvar,
            ze=ze, woff=woff
        )
    
    def _get_2d_irf(self, irf_name):
        az_min, az_max, pedvar, woff, ze, az, \
        twod_irf, x_bins, x_mins, \
        x_maxs, y_bins, y_mins, y_maxs = self._extract(
            "azMin", "azMax", "pedvar", "Woff", "ze", "az",
            f"{irf_name}_value", f"{irf_name}_binsx",
            f"{irf_name}_minx", f"{irf_name}_maxx",
            f"{irf_name}_binsy", f"{irf_name}_miny", f"{irf_name}_maxy"
        )
        az_centers = (az_max + az_min) / 2.0
        az_centers[az_centers < 0] += 360
        twod_irf = self._to_uniform_numpy(twod_irf)
        twod_irf = twod_irf.reshape((-1, x_bins[0], y_bins[0]))
        
        x_axis = self._get_2d_irf_axis(x_mins[0], x_maxs[0], int(x_bins[0]))
        y_axis = self._get_2d_irf_axis(y_mins[0], y_maxs[0], int(y_bins[0]))
        
        irf_ds = xr.Dataset(
            {
                "irf": (["index", "irf_x", "irf_y"], twod_irf),
            },
            coords={
                "irf_x": x_axis,
                "irf_y": y_axis,
            },
            attrs={
                "irf_filename": self.file_name,
                "irf_name": irf_name,  # Additional data
            },
        )
        
        return self._convert_index_to_dims(
            irf_ds, az=az_centers, pedvar=pedvar,
            ze=ze, woff=woff
        )
    
    def __init__(self, in_filename):
        self.file_name = in_filename
        self.tree_ptr = uproot.open(f"{in_filename}:fEffAreaH2F")
    
    def __call__(self, irf_name):
        try:
            if irf_name in self.implemented_irf_names_1d:
                return self._get_1d_irf(irf_name)
            if irf_name in self.implemented_irf_names_2d:
                return self._get_2d_irf(irf_name)
        except uproot.KeyInFileError:
            raise KeyError(
                f"IRF not found! {irf_name}; "
                f"All keys in {self.file_name}:fEffAreaH2F; "
                f"{self.tree_ptr.keys()}"
            )
        raise ValueError(f"Unknown IRF: {irf_name}!")


class XIRFManager:
    def __init__(self, in_filename):
        self.cache = {}
        self.irf_extractor = None
        self.filename = in_filename
    
    @property
    def filename(self):
        return self._filename
    
    @filename.setter
    def filename(self, value):
        self._filename = value
        if value is not None and os.path.isfile(value):
            self.irf_extractor = XIRFExtractor(value)
        else:
            self._filename = None
            self.irf_extractor = None
            self.cache = {}
    
    def clear_cache(self):
        self.cache = {}
    
    def interpolate(self, irf_name, **kwargs):
        """A wrapper function for xarray.interp
         passes all keyword arguments to xarray.interp
         for the generated irf dataset object"""
        if self.irf_extractor is None:
            raise ValueError("Please specify a IRF root file first!")
        
        if irf_name not in self.cache:
            self.cache[irf_name] = self.irf_extractor(irf_name)
        c_irf = self.cache[irf_name]
        
        return c_irf.interp(**kwargs)
        