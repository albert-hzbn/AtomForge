"""Optional third-party packages: tests needing one skip when it is missing.

    @requires("pymatgen")
    def test_against_pymatgen(self): ...

    if missing("seekpath"): self.skipTest(...)
"""
import importlib.util
import unittest


def missing(*modules):
    """The modules of `modules` that cannot be imported (top-level names)."""
    return [m for m in modules if importlib.util.find_spec(m.split(".")[0]) is None]


def requires(*modules):
    """Skip the decorated test (or class) unless every module is installed."""
    absent = missing(*modules)
    return unittest.skipIf(bool(absent), "needs " + ", ".join(absent) + " (optional)")
