from pyV2DL3.root_lib_util import SingletonDecorator
from ROOT import gSystem
import logging
logger = logging.getLogger(__name__)


class EDStatus:
    __ed_loaded = False

    def __init__(self):
        pass

    def loadED(self):
        if(not self.__ed_loaded):
            logger.debug('Here we should load eventDisplay libVAnaSum.so')
            logger.debug('Unfortunately, I\'m unable to compile ED with a working root+python+root_numpy...')
            logger.debug('For now, will just use ROOT.')
            # if gSystem.Load("libVEGASCommon.dylib"):
            #     raise Exception("Problem loading VEGAS Common libraries - please check this before proceeding")
            # if gSystem.Load("libVEGASStage6.dylib"):
            #     raise Exception("Problem loading VEGAS Stage 6 libraries - please check this before proceeding")
            # if gSystem.Load("libVEGASStage5.dylib"):
            #     raise Exception("Problem loading VEGAS Stage 5 libraries - please check this before proceeding")
            self.__ed_loaded = True


EDStatus = SingletonDecorator(EDStatus)
