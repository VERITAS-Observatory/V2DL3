import contextlib

from ROOT import gROOT


@contextlib.contextmanager
def cpp_print_context(verbose=True):
    if not verbose:
        gROOT.ProcessLine("std::cout.setstate(std::ios_base::failbit)")
    yield
    if not verbose:
        gROOT.ProcessLine("std::cout.clear()")


class SingletonDecorator:
    def __init__(self, klass):
        self.klass = klass
        self.instance = None

    def __call__(self, *args, **kwds):
        if self.instance is None:
            self.instance = self.klass(*args, **kwds)
        return self.instance
