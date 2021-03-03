import ROOT
from pyV2DL3.VtsDataSource import VtsDataSource
from pyV2DL3.eventdisplay.fillEVENTS import __fillEVENTS__
from pyV2DL3.eventdisplay.fillRESPONSE import __fillRESPONSE__


class EventDisplayDataSource(VtsDataSource):
    def __init__(self, etv_file, ea_file):
        super(EventDisplayDataSource, self).__init__('EventDisplay', etv_file, ea_file)

        self.__evt_file__ = etv_file
        self.__ea_file__ = ROOT.TFile.Open(ea_file)

        # Auxiliary storage 
        self.__azimuth__ = 0
        self.__zenith__ = 0
        self.__noise__ = 0
        self.__offset__= 0 #new

    def __fill_evt__(self):
        # can be simplified further:
        gti, ea_config, events = __fillEVENTS__(self.__evt_file__)
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
        print("Coordinates to fillresponse:", az, ze, nn, oo)
        self.__response__ = __fillRESPONSE__(self.__evt_file__,self.__ea_file__, az, ze, nn, oo,
                                                      self.__irf_to_store__)
