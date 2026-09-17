"""Optional scientific workflows backed by ASE and pymatgen.

Install ``atomforge-py[science]`` for these integrations. Calculators, VASP
datasets and external population-analysis executables remain user-supplied.
All public coordinates are Angstrom and energies are eV unless stated otherwise.
"""

from .simulation import to_ase, from_ase, relax, molecular_dynamics, phonons, mode_frames, load_calculator
from .spectra import read_bands, read_dos, read_projections, plot_bands, plot_dos
from .diffraction import powder_diffraction, single_crystal, plot_diffraction
from .population import bader, ddec, hirshfeld, orbital_density
from .trajectory import read_trajectory, write_trajectory, export_animation
from .reciprocal import BandGrid, read_bandgrid, plot_surface

__all__ = ["to_ase", "from_ase", "relax", "molecular_dynamics", "phonons", "mode_frames", "load_calculator",
           "read_bands", "read_dos", "read_projections", "plot_bands", "plot_dos",
           "powder_diffraction", "single_crystal", "plot_diffraction", "bader", "ddec",
           "hirshfeld", "orbital_density", "read_trajectory", "write_trajectory", "export_animation",
           "BandGrid", "read_bandgrid", "plot_surface"]
