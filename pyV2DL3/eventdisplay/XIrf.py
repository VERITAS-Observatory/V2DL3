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
        
        return self._convert_index_to_dims(
            irf_ds, azimuth=az_centers, pedvar=pedvar,
            zenith=ze, woff=woff
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
            irf_ds, azimuth=az_centers, pedvar=pedvar,
            zenith=ze, woff=woff
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


class InterpolateInput:
    """This class turns the event dictionary into
    a wrapper object for the pandas Dataframe with
    properties given as keys in `exposed`;
    the values of `exposed` tells the class which
    dictionary keys correspond to which property name.
    For more flexibility, if the value is callable then
    the entire dictionary is passed into it and can be manipulated.
    Once built, some nice functions for helping create
    input values for the locations where the data is in the
    IRF hypercube are implemented.
    """
    
    exposed = {
        "azimuth": "AZ",
        "zenith": lambda x: 90. - x["ALT"],
        "pedvar": "PEDVAR",
        "woff": "WOFF"
    }
    
    def __getattr__(self, item):
        return getattr(self.data, item)
    
    def __getitem__(self, item):
        return self.data[item]
    
    def __init__(self, evt_dict):
        dvars_d = {
            k: v(evt_dict) if callable(v) else evt_dict[v]
            for k, v in self.exposed.items()
        }
        self.data = pd.DataFrame(dvars_d.copy())
    
    def limits(self, in_attr):
        c_attr = self[in_attr]
        return c_attr.min(), c_attr.max()
    
    def range(self, in_attr):
        a_min, a_max = self.limits(in_attr)
        return a_max - a_min
    
    def bin(self,
            in_attr,
            n_bins=None,
            bin_width=None,
            equal_sized=False,
            retbins=False
            ):
        
        if not n_bins and not bin_width:
            raise ValueError("Bin Specification Undefined!")
        
        if equal_sized:
            if not n_bins:
                t_width = self.range(in_attr)
                c_bins = int((t_width // bin_width) + 1)
            else:
                c_bins = n_bins
        else:
            c_min, c_max = self.limits(in_attr)
            
            if not bin_width:
                t_width = self.range(in_attr)
                bin_width = t_width / n_bins
            
            c_bins = []
            c_cur = c_min
            while c_cur < c_max:
                c_bins.append(c_cur)
                c_cur += bin_width
            c_bins.append(c_max)
        
        c_series = self[in_attr].copy(deep=True)
        
        c_df = pd.DataFrame({in_attr: c_series})
        if isinstance(c_bins, int):
            labels = list(range(c_bins))
        else:
            labels = list(range(len(c_bins) - 1))
        cut_rv = pd.cut(
            c_series,
            c_bins,
            labels=labels,
            retbins=retbins
        )
        if retbins:
            c_df["bin"] = cut_rv[0]
            return c_df, cut_rv[1], labels
        else:
            c_df["bin"] = cut_rv
            return c_df, labels
    
    @staticmethod
    def sorted_unique(in_arr, d_places=2):
        base = np.float_power(10, d_places)
        return np.sort(np.floor(in_arr * base).unique() / base)
    
    def make_bins(self, in_attr, **kwargs):
        a_data, a_bins = self.bin(in_attr, **kwargs)
        return [
            self.sorted_unique(a_data[a_data["bin"] == bn][in_attr]).mean()
            for bn in a_bins
        ]
    
    def __call__(self, ped_bin=None, az_bin=None, ze_bin=10.):
        if ped_bin is None:
            ped_rv = [self.sorted_unique(self.pedvar).mean()]
        else:
            ped_rv = self.make_bins("pedvar", bin_width=ped_bin)
        
        if az_bin is None:
            az_rv = [self.sorted_unique(self.azimuth).mean()]
        else:
            az_rv = self.make_bins("azimuth", bin_width=az_bin)
        
        if ze_bin is None:
            ze_rv = [self.sorted_unique(self.zenith).mean()]
        else:
            ze_rv = self.make_bins("zenith", bin_width=ze_bin)
        
        woff_rv = self.sorted_unique(self.woff)
        
        return ped_rv, az_rv, ze_rv, woff_rv


class XIRFManager:
    def __init__(self, in_filename, evt_dict):
        self.irf_extractor = None
        self.filename = in_filename
        self.ii = InterpolateInput(evt_dict)
        self.irf_cache = {}
    
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
            self.irf_cache = {}
    
    def interpolate(self, irf_name, **kwargs):
        """A wrapper function for xarray.interp
         passes all keyword arguments to xarray.interp
         for the generated irf dataset object"""
        if self.irf_extractor is None:
            raise ValueError("Please specify a IRF root file first!")
        
        if irf_name in self.irf_cache:
            c_irf = self.irf_cache[irf_name]
        else:
            c_irf = self.irf_extractor(irf_name)
            self.irf_cache[irf_name] = c_irf
        
        return c_irf.interp(**kwargs)
    
    def fill_direction_migration(self):
        ped, az, ze, woff = self.ii()
        
        return self.interpolate(
            "hAngularLogDiffEmc_2D",
            pedvar=ped,
            azimuth=az,
            zenith=ze,
            woff=woff
        )


def load_test():
    from pyV2DL3.eventdisplay.fillEVENTS import __fillEVENTS__
    irf_filename = "effArea-v486-auxv01-CARE_June2020-Cut-NTel2-PointSource-Moderate-TMVA-BDT-GEO-V6_2012_2013a-ATM62-T1234-testCI.root"
    root_file = "64080.anasum.root"
    rv = __fillEVENTS__(root_file)
    event_dict = rv[2]
    irf_ex = XIRFManager(irf_filename, event_dict)
    return irf_ex.fill_direction_migration()
