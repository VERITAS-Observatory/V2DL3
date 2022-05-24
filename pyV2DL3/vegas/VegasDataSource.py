from typing import List
import ROOT

from pyV2DL3.vegas.fillEVENTS_not_safe import __fillEVENTS_not_safe__
from pyV2DL3.vegas.fillRESPONSE_not_safe import __fillRESPONSE_not_safe__
from pyV2DL3.vegas.load_vegas import VEGASStatus
from pyV2DL3.VtsDataSource import VtsDataSource
from pyV2DL3.EventClass import EventClass
from pyV2DL3.vegas.util import loadUserCuts


class VegasDataSource(VtsDataSource):
    def __init__(self, etv_file, ea_file, 
                 event_classes: List[EventClass] = None, 
                 user_cut_file=None,
                 reco_type=1, 
                 store_msw_msl=False,
                 ):
        super(VegasDataSource, self).__init__("VEGAS", etv_file, ea_file)

        # Loading VEGAS if not already done so
        self.vegas_status = VEGASStatus()
        self.vegas_status.loadVEGAS()
        self.__evt_file__ = ROOT.VARootIO(etv_file, True)
        self.__event_classes__ = event_classes
        if user_cut_file is not None:
            # Load user cut data. See loadUserCuts() in util.py for possible keys.
            self.__user_cuts__ = loadUserCuts(user_cut_file)
        else:
            self.__user_cuts__ = None
        self.__reco_type__ = reco_type
        self.__store_msw_msl__ = store_msw_msl
        # Developer exceptions; a user should be unable to trigger these
        # Ensures that this VDS was constructed with an EA xor EventClass(s)
        if ea_file is None and event_classes is None:
            ea_exception = "Running V2DL3 without effective area file(s) is currently unsupported."
            ea_exception += "\n Remove this exception if you are implementing this behavior."
            raise Exception(ea_exception)
        elif ea_file is not None and event_classes is not None:
            raise Exception("VegasDataSource was somehow constructed with both an effective area"
                           +"file and event classes")

        if ea_file is not None:
            self.__ea_file__ = ROOT.VARootIO(ea_file, True)
        elif event_classes is not None:
            for ec in event_classes:
                ec.validate_vegas_cuts(self.__evt_file__)
            
            
        # Auxiliary storage
        self.__azimuth__ = 0
        self.__zenith__ = 0
        self.__noise__ = 0


    def __del__(self):
        """Close the root files
        
        These typechecks will prevent the user from having their true exception
        buried by a CPyCppyy exception on program exit.
        
        This seems to happen when the file is already closed, which 
        seems to happen in most cases upon program exception.
        """
        cpy_nonestring = "<class 'CPyCppyy_NoneType'>"

        if str(type(self.__evt_file__)) != cpy_nonestring:
            self.__evt_file__.closeTheRootFile()

        if self.__ea_file__ is not None:
            if str(type(self.__ea_file__)) != cpy_nonestring:
                self.__ea_file__.closeTheRootFile()


    def __fill_evt__(self):
        gti, ea_config, evts = __fillEVENTS_not_safe__(self.__evt_file__,
                                                       reco_type=self.__reco_type__,
                                                       user_cuts_dict=self.__user_cuts__,
                                                       event_classes=self.__event_classes__,
                                                       store_msw_msl=self.__store_msw_msl__
                                                       )
        self.__gti__ = gti
        # This is an array of dicts for each event class (array of one if not using event classes.)
        self.__evt__ = evts
        self.__azimuth__ = ea_config["azimuth"]
        self.__zenith__ = ea_config["zenith"]
        self.__noise__ = ea_config["noise"]


    def __fill_gti__(self):
        pass


    def __fill_response__(self):
        az = self.__azimuth__
        ze = self.__zenith__
        nn = self.__noise__
        response_dicts = []
        # Event class mode
        if self.__event_classes__ is not None:
            # Setup
            multi_eclass = False
            if len(self.__event_classes__) > 1:
                multi_eclass = True
            eclass_idx = 0
            # Fill responses per event class
            for ec in self.__event_classes__:
                ea = ec.effective_area_IO
                response_dicts.append(
                    __fillRESPONSE_not_safe__(ea, az, ze, nn, self.__irf_to_store__,
                                              event_class_idx=eclass_idx, 
                                              multi_eclass=multi_eclass,
                                              msw_range=(ec.msw_lower, ec.msw_upper),
                                              )
                )
                eclass_idx += 1
        # When not using event classes
        else:
            response_dicts.append( 
                __fillRESPONSE_not_safe__(self.__ea_file__, az, ze, 
                                          nn, self.__irf_to_store__)
            )
        self.__response__ = response_dicts
