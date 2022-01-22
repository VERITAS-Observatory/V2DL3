import ROOT

from pyV2DL3.VtsDataSource import VtsDataSource, VtsDataSourceMeta
from pyV2DL3.vegas.load_vegas import VEGASStatus


class VegasDataSource(VtsDataSource, metaclass=VtsDataSourceMeta, file_type="VEGAS"):
    def __init__(self, etv_file, ea_file):
        super(VegasDataSource, self).__init__('VEGAS')
        
        # Loading VEGAS if not already done so
        self.vegas_status = VEGASStatus()
        self.vegas_status.loadVEGAS()
        self.__evt_file__ = ROOT.VARootIO(etv_file, True)
        self.__ea_file__ = ROOT.VARootIO(ea_file, True)
    
    def __del__(self):
        # Close the root files
        self.__evt_file__.closeTheRootFile()
        self.__ea_file__.closeTheRootFile()
