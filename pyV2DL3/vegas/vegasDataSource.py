from pyV2DL3.vtsDataSource import vtsDataSource
from pyV2DL3.vegas.load_vegas import VEGASStatus    
from pyV2DL3.vegas.fillEVENTS_not_safe import __fillEVENTS_not_safe__
from pyV2DL3.vegas.fillRESPONSE_not_safe import __fillRESPONSE_not_safe__
import ROOT

class vegasDataSource(vtsDataSource):
    def __init__(self,etv_file,ea_file):
        super(vegasDataSource, self).__init__(etv_file,ea_file)

        # Loading VEGAS if not already done so
        self.vegas_status = VEGASStatus()        
        self.vegas_status.loadVEGAS()
        self.__evt_file__ = ROOT.VARootIO(etv_file,True)
        self.__ea_file__  = ROOT.VARootIO(ea_file,True)

        # Auxiliary storage 
        self.__azimuth__ = 0
        self.__zenith__  = 0
        self.__noise__   = 0 

    def __fill_evt__(self):
       gti,ea_config,evts = __fillEVENTS_not_safe__(self.__evt_file__)
       self.__gti__ = gti
       self.__evt__ = evts
       self.__azimuth__ = ea_config['azimuth'] 
       self.__zenith__ = ea_config['zenith'] 
       self.__noise__   = ea_config['noise']
   
    def __fill_gti__(self):
        pass

    def __fill_response__(self):
        az = self.__azimuth__ 
        ze = self.__zenith__
        nn = self.__noise__  
        self.__response__ = __fillRESPONSE_not_safe__(self.__ea_file__,az,ze,nn,0.5)    
        
