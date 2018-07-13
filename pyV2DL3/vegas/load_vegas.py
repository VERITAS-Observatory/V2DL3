from pyV2DL3.root_lib_util import SingletonDecorator
from ROOT import gSystem
import logging
logger = logging.getLogger(__name__)
class VEGASStatus:
    __vegas_loaded = False
    def __init__(self):
        pass

    def loadVEGAS(self):
        if(not self.__vegas_loaded):
            logger.debug('Load VEGAS lib')
            if gSystem.Load("libVEGAScommon"):
                raise Exception("Problem loading VEGAS Common libraries - please check this before proceeding")
            if gSystem.Load("libVEGASstage5"):
                raise Exception("Problem loading VEGAS Stage 5 libraries - please check this before proceeding")
            if gSystem.Load("libVEGASstage6"):
                raise Exception("Problem loading VEGAS Stage 6 libraries - please check this before proceeding") 
            self.__vegas_loaded = True

VEGASStatus = SingletonDecorator(VEGASStatus)
