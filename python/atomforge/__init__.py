"""
atomforge — Python interface for AtomForge structures.

Quick-start
-----------
>>> import atomforge as af
>>> s = af.load("my_crystal.cif")
>>> s.repeat(2, 2, 1).view()          # open in AtomForge GUI

>>> s2 = af.Structure()
>>> s2.add_atom("Fe", 0, 0, 0)
>>> s2.add_atom("C",  1.8, 0, 0)
>>> s2.save("dimer.xyz")
"""

from ._structure import Atom, Structure
from ._io import load, save
from ._viewer import view

__all__ = ["Atom", "Structure", "load", "save", "view"]
__version__ = "0.1.0"
