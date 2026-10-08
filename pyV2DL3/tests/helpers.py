"""Shared test doubles for ROOT objects."""


class RootFile(dict):
    """Minimal context-manager-compatible stand-in for an uproot file."""

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class RootLog:
    """Minimal stand-in for a ROOT log object exposing ``fLines``."""

    def __init__(self, lines):
        self.lines = lines

    def member(self, name):
        if name != "fLines":
            raise KeyError(name)
        return self.lines
