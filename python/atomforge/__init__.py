"""
atomforge — Python interface for AtomForge structures.

No external dependencies required.  Supported formats (read + write):
XYZ, extended-XYZ, VASP POSCAR/CONTCAR, PDB, CIF (P1), LAMMPS data.

Quick-start
-----------
>>> import atomforge as af
>>> s = af.load("crystal.cif")
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
