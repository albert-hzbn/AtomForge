"""Structure relaxation, dynamics and finite-displacement phonons with ASE."""

from dataclasses import dataclass
import math
from pathlib import Path


def load_calculator(module, attribute, *args, **kwargs):
    """Import ``module`` and return ``getattr(module, attribute)(*args, **kwargs)``.

    A dotted-path convenience for handing relax()/molecular_dynamics()/phonons()
    any pip-installed ASE calculator by name -- a DFT calculator, or a machine-
    learned interatomic potential such as MACE (module="mace.calculators",
    attribute="mace_mp"), CHGNet ("chgnet.model.dynamics", "CHGNetCalculator"),
    MatterSim ("mattersim.forcefield", "MatterSimCalculator") or ORB
    ("orb_models.forcefield.calculator", "ORBCalculator", model). AtomForge does
    not install, select or validate the underlying package, model or weights;
    consult that package's own documentation for arguments, accuracy and
    hardware requirements. Import and construction errors propagate unchanged.
    """
    import importlib
    if not isinstance(module, str) or not module or not isinstance(attribute, str) or not attribute:
        raise ValueError("module and attribute must be nonempty strings")
    factory = getattr(importlib.import_module(module), attribute)
    return factory(*args, **kwargs)


def to_ase(structure):
    """Copy species, Cartesian positions and cell into an ASE Atoms object."""
    from ase import Atoms
    return Atoms(symbols=[a.symbol for a in structure.atoms],
                 positions=[(a.x, a.y, a.z) for a in structure.atoms],
                 cell=structure.cell, pbc=structure.cell is not None)


def from_ase(atoms):
    """Copy species, positions and a fully periodic cell into an AtomForge Structure.

    Mixed periodic/nonperiodic axes cannot be represented by Structure and are
    rejected instead of silently converting a slab or wire to 3D periodicity.
    Nonperiodic bounding boxes are not periodic unit cells.
    """
    from .._structure import Structure
    result = Structure()
    for symbol, point in zip(atoms.get_chemical_symbols(), atoms.positions):
        result.add_atom(symbol, *point)
    if atoms.pbc.any() and (not atoms.pbc.all() or atoms.cell.rank != 3):
        raise ValueError("AtomForge Structure requires either no periodic axes or a full 3D periodic cell")
    if atoms.pbc.all():
        result.cell = atoms.cell.tolist()
    return result


@dataclass
class RelaxationResult:
    structure: object
    converged: bool
    energy: float
    forces: object
    steps: int


def relax(structure, calculator, *, fmax=0.02, steps=200, relax_cell=False, trajectory=None):
    """Run BFGS with an explicitly supplied ASE calculator (including DFT engines).

    ``converged`` must be checked: reaching the step limit is not convergence.
    Calculator setup, pseudopotentials and electronic convergence are controlled
    by the caller. The source Structure is never modified.
    """
    if not math.isfinite(fmax) or fmax <= 0 or not isinstance(steps, int) or steps < 1:
        raise ValueError("Use positive fmax and step count")
    from ase.optimize import BFGS
    atoms = to_ase(structure)
    atoms.calc = calculator
    target = atoms
    if relax_cell:
        if not structure.cell:
            raise ValueError("Cell relaxation requires a periodic cell")
        from ase.filters import FrechetCellFilter
        target = FrechetCellFilter(atoms)
    optimizer = BFGS(target, trajectory=str(trajectory) if trajectory else None, logfile=None)
    converged = optimizer.run(fmax=fmax, steps=steps)
    return RelaxationResult(from_ase(atoms), bool(converged), float(atoms.get_potential_energy()),
                            atoms.get_forces().tolist(), optimizer.nsteps)


def molecular_dynamics(structure, calculator, *, steps=100, timestep_fs=1, temperature_K=300,
                       seed=0, sample_interval=1, velocities=None):
    """Run NVE velocity-Verlet; return frames and total energies in eV.

    With velocities omitted, initialize a Maxwell-Boltzmann distribution using
    seed. Supplied velocities use ASE's internal velocity units. Temperature is
    an initialization target, not a thermostat. Inspect energy drift to select dt.
    """
    if not isinstance(steps, int) or steps < 1 or not isinstance(sample_interval, int) or sample_interval < 1:
        raise ValueError("Steps and sampling interval must be positive integers")
    if not math.isfinite(timestep_fs) or timestep_fs <= 0 or not math.isfinite(temperature_K) or temperature_K < 0:
        raise ValueError("Invalid timestep or temperature")
    import numpy as np
    from ase import units
    from ase.md.verlet import VelocityVerlet
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    atoms = to_ase(structure)
    atoms.calc = calculator
    if velocities is None:
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K, rng=np.random.RandomState(seed))
    else:
        values = np.asarray(velocities, dtype=float)
        if values.shape != (len(atoms), 3) or not np.isfinite(values).all():
            raise ValueError("Velocities must have shape (atoms, 3) and be finite")
        atoms.set_velocities(values)
    frames, energies = [], []
    dynamics = VelocityVerlet(atoms, timestep=timestep_fs * units.fs, logfile=None)
    def record():
        frames.append(from_ase(atoms))
        energies.append(float(atoms.get_potential_energy() + atoms.get_kinetic_energy()))
    dynamics.attach(record, interval=sample_interval)
    dynamics.run(steps)
    return {"frames": frames, "total_energies_eV": energies,
            "sample_interval_fs": sample_interval * timestep_fs}


def phonons(structure, calculator, qpoints, directory, *, supercell=(2,2,2), displacement=0.01):
    """Finite-displacement phonons; energies in eV and complex eigenvectors.

    qpoints are fractional reciprocal coordinates. Negative energies denote
    imaginary frequencies. A fresh cache directory prevents stale force reuse.
    """
    import numpy as np
    from ase.phonons import Phonons
    if not structure.cell:
        raise ValueError("Phonons require a periodic cell")
    if len(supercell) != 3 or any(not isinstance(n,int) or n < 1 for n in supercell):
        raise ValueError("Supercell must contain three positive integers")
    if not math.isfinite(displacement) or displacement <= 0:
        raise ValueError("Displacement must be positive")
    qpoints = np.asarray(qpoints, dtype=float)
    if qpoints.ndim != 2 or qpoints.shape[1] != 3 or not np.isfinite(qpoints).all():
        raise ValueError("qpoints must have shape (n,3)")
    folder = Path(directory)
    folder.mkdir(parents=True, exist_ok=False)
    calculation = Phonons(to_ase(structure), calculator, supercell=supercell,
                          delta=displacement, name=str(folder / "forces"))
    calculation.run()
    calculation.read(acoustic=True)
    energies, modes = calculation.band_structure(qpoints, modes=True)
    return {"qpoints": qpoints, "energies_eV": energies, "modes": modes}


def mode_frames(structure, mode, *, amplitude=0.2, frames=40):
    """Animate one Gamma-point mode, with maximum displacement amplitude in A."""
    import numpy as np
    values = np.asarray(mode, dtype=complex)
    if values.shape != (len(structure),3) or not np.isfinite(values).all():
        raise ValueError("Mode must contain one finite Cartesian vector per atom")
    if not isinstance(frames,int) or frames < 2 or not math.isfinite(amplitude) or amplitude <= 0:
        raise ValueError("Invalid frame count or amplitude")
    norm = float(np.max(np.linalg.norm(values, axis=1)))
    if norm == 0:
        raise ValueError("Cannot animate a zero mode")
    result = []
    for phase in np.linspace(0,2*np.pi,frames,endpoint=False):
        frame = structure.copy()
        displacement = (values * np.exp(1j*phase)).real * amplitude / norm
        for atom, vector in zip(frame.atoms, displacement):
            atom.x += vector[0]; atom.y += vector[1]; atom.z += vector[2]
        result.append(frame)
    return result
