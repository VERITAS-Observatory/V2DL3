import ROOT
from ROOT import gSystem,gROOT
import contextlib
@contextlib.contextmanager
def CppPrintContext(verbose=True):
    if(not verbose):
        gROOT.ProcessLine("std::cout.setstate(std::ios_base::failbit)")
    else:
        pass
    yield
    if(not verbose):
        gROOT.ProcessLine("std::cout.clear()")

class SingletonDecorator:
    def __init__(self,klass):
        self.klass = klass
        self.instance = None
    def __call__(self,*args,**kwds):
        if self.instance == None:
            self.instance = self.klass(*args,**kwds)
        return self.instance

class VEGASStatus:
    __vegas_loaded = False
    def __init__(self):
        pass

    def loadVEGAS(self):
        if(not self.__vegas_loaded):
            if gSystem.Load("libVEGASCommon.dylib"):
                raise Exception("Problem loading VEGAS Common libraries - please check this before proceeding")
            if gSystem.Load("libVEGASStage6.dylib"):
                raise Exception("Problem loading VEGAS Stage 6 libraries - please check this before proceeding")
            if gSystem.Load("libVEGASStage5.dylib"):
                raise Exception("Problem loading VEGAS Stage 5 libraries - please check this before proceeding") 
            self.__vegas_loaded = True

VEGASStatus = SingletonDecorator(VEGASStatus)
