from pyV2DL3.root_lib_util import SingletonDecorator
from ROOT import gSystem
import logging
logger = logging.getLogger(__name__)


class EDStatus:
    __ed_loaded = False

    def __init__(self):
        pass

    def load_ed(self):
        if not self.__ed_loaded:
            logger.debug('Load eventDisplay libVAnaSum.so')
            if gSystem.Load("$EVNDISPSYS/lib/libVAnaSum.so"):
                raise Exception("Problem loading eventDisplay libraries - please check this before proceeding")
            self.__ed_loaded = True


EDStatus = SingletonDecorator(EDStatus)