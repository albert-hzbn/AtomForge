"""Explicit-calculator migration paths and controlled MD ensembles."""

from ._validation import array, boolean, integer, positive
from .simulation import from_ase, to_ase


def migration_path(initial, final, calculator_factory, *, images=7, fmax=.03,
                   steps=300, spring_eV_per_A2=.1, climb=True, mic=False):
    """Optimize a fixed-cell climbing-image NEB using ASE's improved tangent.

    Endpoints must already be relaxed, have identical cell/species/order, and
    describe the intended atom mapping. Factory is called separately for EVERY
    image (external calculators must also use distinct working directories).
    Barriers are sampled image energies; convergence and endpoint forces are
    reported rather than assumed. mic=True chooses shortest periodic paths.
    """
    import numpy as np
    from ase.mep import NEB
    from ase.optimize import FIRE
    climb = boolean(climb, "climb")
    mic = boolean(mic, "mic")
    integer(images, "images", 3); integer(steps, "steps")
    positive(fmax, "fmax"); positive(spring_eV_per_A2, "spring_eV_per_A2")
    start, end = to_ase(initial), to_ase(final)
    if not len(start) or start.get_chemical_symbols() != end.get_chemical_symbols():
        raise ValueError("NEB requires identical ordered species at both endpoints")
    if not np.array_equal(start.pbc, end.pbc) or not np.allclose(start.cell, end.cell, rtol=0, atol=1e-10):
        raise ValueError("Variable-cell NEB is not supported; endpoint cells must match")
    chain = [start] + [start.copy() for _ in range(images - 2)] + [end]
    calculators = []
    for image in chain:
        calculator = calculator_factory()
        if calculator is None or any(calculator is previous for previous in calculators):
            raise ValueError("calculator_factory must return an independent calculator for each image")
        image.calc = calculator
        calculators.append(calculator)
    path = NEB(chain, k=spring_eV_per_A2, climb=bool(climb), method="improvedtangent")
    path.interpolate(mic=bool(mic))
    optimizer = FIRE(path, logfile=None)
    converged = optimizer.run(fmax=fmax, steps=steps)
    energies = np.array([image.get_potential_energy() for image in chain])
    endpoint_forces = [float(np.linalg.norm(image.get_forces(), axis=1).max()) for image in (start, end)]
    return {"images": [from_ase(image) for image in chain], "energies_eV": energies,
            "forward_barrier_eV": float(energies.max() - energies[0]),
            "reverse_barrier_eV": float(energies.max() - energies[-1]),
            "reaction_energy_eV": float(energies[-1] - energies[0]),
            "converged": bool(converged), "steps": optimizer.nsteps,
            "endpoint_max_forces_eV_per_A": endpoint_forces}


def _dynamics(structure, calculator, *, ensemble, steps, timestep_fs, temperature_K,
              seed, sample_interval, thermostat_fs, pressure_GPa=0., barostat_fs=1000.):
    import numpy as np
    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    integer(steps, "steps"); integer(sample_interval, "sample_interval")
    integer(seed, "seed", 0)
    dt = positive(timestep_fs, "timestep_fs")
    temperature = positive(temperature_K, "temperature_K")
    damping = positive(thermostat_fs, "thermostat_fs")
    atoms = to_ase(structure)
    if len(atoms) < 2:
        raise ValueError("At least two atoms are required")
    atoms.calc = calculator
    rng = np.random.RandomState(seed)
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature, rng=rng)
    if ensemble == "nvt":
        from ase.md.langevin import Langevin
        dynamics = Langevin(atoms, dt * units.fs, temperature_K=temperature,
                            friction=1/(damping * units.fs), rng=rng, fixcm=False, logfile=None)
    else:
        from ase.md.nose_hoover_chain import IsotropicMTKNPT
        if not atoms.pbc.all() or atoms.cell.rank != 3:
            raise ValueError("NPT requires a full periodic cell")
        positive(barostat_fs, "barostat_fs")
        pressure = float(array([pressure_GPa], "pressure_GPa", 1)[0])
        # Early capability check gives an actionable error for calculators
        # providing energies and forces but no stress.
        atoms.get_stress()
        dynamics = IsotropicMTKNPT(atoms, dt * units.fs, temperature_K=temperature,
            pressure_au=pressure * units.GPa, tdamp=damping * units.fs,
            pdamp=barostat_fs * units.fs, logfile=None)
    frames, times, temperatures, energies, volumes, velocities, pressures = [], [], [], [], [], [], []
    def record():
        energy = float(atoms.get_total_energy())
        if not np.isfinite(energy) or not np.isfinite(atoms.positions).all():
            raise RuntimeError("Dynamics became nonfinite; inspect timestep and calculator")
        frames.append(from_ase(atoms)); times.append(dynamics.nsteps * dt)
        temperatures.append(float(atoms.get_temperature())); energies.append(energy)
        volumes.append(float(atoms.get_volume()) if atoms.cell.rank == 3 else None)
        # ASE internal velocities -> A/fs, directly usable by VACF and spectra.
        velocities.append(atoms.get_velocities() * units.fs)
        if ensemble == "npt":
            pressures.append(float(-np.trace(atoms.get_stress(voigt=False, include_ideal_gas=True))/3/units.GPa))
    dynamics.attach(record, interval=sample_interval)
    dynamics.run(steps)
    if dynamics.nsteps % sample_interval:
        record()
    return {"frames": frames, "time_fs": np.asarray(times), "temperature_K": np.asarray(temperatures),
            "total_energy_eV": np.asarray(energies), "volume_A3": volumes,
            "velocities_A_per_fs": np.asarray(velocities), "pressure_GPa": pressures,
            "ensemble": ensemble, "target_temperature_K": temperature}


def nvt_dynamics(structure, calculator, *, steps=1000, timestep_fs=1., temperature_K=300.,
                 thermostat_fs=100., seed=0, sample_interval=10):
    """Langevin NVT with a thermostat relaxation time in fs and reproducible seed.

    Returns temperatures, physical energies, frames and velocities. Physical
    energy is not conserved under a thermostat. COM is not constrained; remove
    any unwanted drift explicitly during analysis. Equilibration is not inferred.
    """
    return _dynamics(structure, calculator, ensemble="nvt", steps=steps, timestep_fs=timestep_fs,
        temperature_K=temperature_K, thermostat_fs=thermostat_fs, seed=seed, sample_interval=sample_interval)


def npt_dynamics(structure, calculator, *, steps=1000, timestep_fs=1., temperature_K=300.,
                 pressure_GPa=0., thermostat_fs=100., barostat_fs=1000., seed=0, sample_interval=10):
    """Isotropic Martyna-Tobias-Klein NPT via ASE, pressure positive in compression.

    Cell volume changes while shape remains fixed. A stress-capable calculator
    is required. Finite runs must be checked for equilibration, fluctuations and
    integration convergence; instantaneous T and P need not equal targets.
    """
    return _dynamics(structure, calculator, ensemble="npt", steps=steps, timestep_fs=timestep_fs,
        temperature_K=temperature_K, pressure_GPa=pressure_GPa, thermostat_fs=thermostat_fs,
        barostat_fs=barostat_fs, seed=seed, sample_interval=sample_interval)


def reciprocal_path(structure, *, spacing_inv_A=.025, symprec_A=1e-5, time_reversal=True):
    """SeeK-path HPKOT high-symmetry path in the standardized primitive cell.

    Returns the primitive structure, reciprocal basis and segment endpoints so
    fractional k coordinates cannot be confused with the input supercell basis.
    Nonmagnetic spatial symmetry is used; time_reversal controls path completion,
    not a magnetic-space-group analysis. Disconnected segments remain explicit.
    """
    import numpy as np
    import seekpath
    time_reversal = boolean(time_reversal, "time_reversal")
    from ase import Atoms
    positive(spacing_inv_A, "spacing_inv_A"); positive(symprec_A, "symprec_A")
    atoms = to_ase(structure)
    if not atoms.pbc.all() or atoms.cell.rank != 3 or not len(atoms):
        raise ValueError("Reciprocal paths require a nonempty periodic structure")
    result = seekpath.get_explicit_k_path((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers),
        with_time_reversal=bool(time_reversal), reference_distance=spacing_inv_A, symprec=symprec_A)
    primitive = Atoms(numbers=result["primitive_types"], scaled_positions=result["primitive_positions"],
                      cell=result["primitive_lattice"], pbc=True)
    return {"primitive_structure": from_ase(primitive), "spacegroup_number": result["spacegroup_number"],
            "spacegroup_symbol": result["spacegroup_international"], "path": result["path"],
            "special_points_fractional": result["point_coords"],
            "kpoints_fractional": result["explicit_kpoints_rel"],
            "kpoints_inv_A": result["explicit_kpoints_abs"],
            "distance_inv_A": result["explicit_kpoints_linearcoord"],
            "labels": result["explicit_kpoints_labels"], "segments": result["explicit_segments"],
            "reciprocal_lattice_inv_A": np.asarray(result["reciprocal_primitive_lattice"])}
