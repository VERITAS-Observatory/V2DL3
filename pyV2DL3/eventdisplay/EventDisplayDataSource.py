from pyV2DL3.VtsDataSource import VtsDataSource
#from pyV2DL3.eventdisplay.load_eventDisplay import EDStatus
from pyV2DL3.eventdisplay.fillEVENTS_not_safe import __fillEVENTS_not_safe__
from pyV2DL3.eventdisplay.fillRESPONSE_not_safe import __fillRESPONSE_not_safe__
import ROOT


class EventDisplayDataSource(VtsDataSource):
    def __init__(self, etv_file, ea_file):
        print("filepath? or name?",etv_file)
        super(EventDisplayDataSource, self).__init__('EventDisplay', etv_file, ea_file)

        # Loading eventDisplay if not already done so
        #self.ed_status = EDStatus()
        #self.ed_status.load_ed()
        #ROOT not needed to open evt file anymore
        #self.__evt_file__ = ROOT.TFile.Open(etv_file)
        self.__evt_file__ = etv_file
        self.__ea_file__ = ROOT.TFile.Open(ea_file)

        # Auxiliary storage 
        self.__azimuth__ = 0
        self.__zenith__ = 0
        self.__noise__ = 0
        self.__offset__= 0 #new

    def __fill_evt__(self):
        #can be simplified further:
        gti, ea_config, events = __fillEVENTS_not_safe__(self.__evt_file__)
        self.__gti__ = gti
        self.__evt__ = events
        self.__azimuth__ = ea_config['azimuth']
        self.__zenith__ = ea_config['zenith']
        self.__noise__ = ea_config['noise']
    def __fill_gti__(self):
        pass

    def __fill_response__(self):
        pass
        az = self.__azimuth__
        ze = self.__zenith__
        nn = self.__noise__
        oo = self.__offset__ #new
        print ("coordinates to pass to fillresponse:",az,ze,nn,oo)
        # TODO: for now, just on the 0.5 deg offset. Improve once we have all IRFs available
        self.__response__ = __fillRESPONSE_not_safe__(self.__ea_file__, az, ze, nn, oo,
                                                      self.__irf_to_store__)
