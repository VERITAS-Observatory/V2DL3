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

