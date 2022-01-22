from pyV2DL3.eventdisplay.fillEVENTS import __fillEVENTS__
from pyV2DL3.eventdisplay.fillRESPONSE import __fillRESPONSE__
from pyV2DL3.vegas.fillEVENTS_not_safe import __fillEVENTS_not_safe__
from pyV2DL3.vegas.fillRESPONSE_not_safe import __fillRESPONSE_not_safe__


def empty_function(*args, **kwargs):
    None


class VtsDataSourceMeta(type):
    
    @staticmethod
    def build_series_function(*in_sub_fns):
        def new_fn(self, *args, **kwargs):
            for sub_fn in in_sub_fns:
                getattr(self, sub_fn)(*args, **kwargs)
        
        return new_fn
    
    @staticmethod
    def build_fill_response_fn(in_ft):
        vegas_attrs = ("__ea_file__", "__azimuth__", "__zenith__", "__noise__", "__irf_to_store__")
        ed_attrs = ("__evt_file__", "__ea_file__", "__azimuth__", "__zenith__", "__pedvar__", "__offset__", "__irf_to_store__")
        
        in_fill_fn = __fillRESPONSE_not_safe__ if in_ft == "VEGAS" else __fillRESPONSE__
        in_obj_attrs = vegas_attrs if in_ft == "VEGAS" else ed_attrs
        
        def fill_response(self, *_):
            attr_in = tuple(getattr(self, aa) for aa in in_obj_attrs)
            self.__response__ = in_fill_fn(*attr_in)
        
        return fill_response
    
    @staticmethod
    def build_fill_evt_fn(in_ft):
        in_fill_evt = __fillEVENTS_not_safe__ if in_ft == "VEGAS" else __fillEVENTS__
        
        def fill_evt(self, *args, **kwargs):
            self.gti, self.ea_config, self.evts = in_fill_evt(self.__evt_file__, *args, **kwargs)
        
        return fill_evt
    
    def __new__(mcs, clsname, bases, attrs, file_type="VEGAS"):
        attrs['fill_data'] = mcs.build_series_function("__fill_evt__", "__fill_gti__", "__fill_response__")
        
        res_cls_fn = "_fill_response" if "__fill_response__" in attrs else "__fill_response__"
        evt_cls_fn = "_fill_evt" if "__fill_evt__" in attrs else "__fill_evt__"
        attrs[res_cls_fn] = mcs.build_fill_response_fn(file_type)
        attrs[evt_cls_fn] = mcs.build_fill_evt_fn(file_type)
        
        return type.__new__(mcs, clsname, bases, attrs)


class VtsDataSource(object):
    __slots__ = (
        "__irf_to_store__", "__data_source_name__", "__evt__",
        "__gti__", "__response__", "__azimuth__", "__zenith__",
        "__pedvar__", "__offset__", "__noise__", "__evt_file__",
        "__ea_file__"
    )
    
    def __init__(self, source_name):
        self.__irf_to_store__ = {'point-like': True, 'full-enclosure': False}
        self.__data_source_name__ = source_name
        print("Reconstruction type:", source_name)
        
        self.__evt__ = {}
        self.__gti__ = {}
        self.__response__ = {}
        
        # Auxiliary storage
        self.__azimuth__ = 0
        self.__zenith__ = 0
        self.__pedvar__ = 0
        self.__offset__ = 0
        self.__noise__ = 0
    
    __fill_gti__ = empty_function
    
    @property
    def gti(self):
        return self.__gti__
    
    @gti.setter
    def gti(self, other):
        self.__gti__ = other
    
    @property
    def ea_config(self):
        return {"azimuth": self.__azimuth__, "zenith": self.__zenith__, "noise": self.__noise__}
    
    @ea_config.setter
    def ea_config(self, other):
        self.__azimuth__ = other['azimuth']
        self.__zenith__ = other['zenith']
        self.__noise__ = other['noise']
    
    @property
    def evts(self):
        return self.__evt__
    
    @evts.setter
    def evts(self, other):
        self.__evt__ = other
    
    def get_source_name(self):
        return self.__data_source_name__
    
    def get_evt_data(self):
        return self.__evt__
    
    def get_gti_data(self):
        return self.__gti__
    
    def get_response_data(self):
        return self.__response__
    
    def set_irfs_to_store(self, irf_to_store):
        self.__irf_to_store__ = irf_to_store
