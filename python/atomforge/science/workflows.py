"""Shared, allowlisted science workflows for Python, CLI and desktop clients.

Data arguments accept inline arrays or {"file": "path", "field": "positions"}
references. File references are resolved relative to the request, never evaluated
as Python. A calculator explicitly names an ASE-compatible factory and kwargs.
"""

from pathlib import Path
import inspect
import json
import math

from . import dynamics_analysis as dynamics
from . import local_structure as local
from . import electronic_properties as electronic
from . import thermomechanics as thermal
from . import advanced_simulation as simulation


# Public order is the implementation priority agreed for this feature set.
TOOLS = {
    "msd": ("Mean-square displacement", dynamics.mean_square_displacement),
    "diffusion": ("Diffusion coefficient", dynamics.diffusion_coefficient),
    "vacf": ("Velocity autocorrelation", dynamics.velocity_autocorrelation),
    "vibrational-spectrum": ("Trajectory vibrational spectrum", dynamics.vibrational_spectrum),
    "local-strain": ("Local strain and D2min", local.local_strain),
    "centrosymmetry": ("Centrosymmetry parameter", local.centrosymmetry),
    "bond-order": ("Steinhardt bond order", local.bond_order),
    "wigner-seitz": ("Wigner-Seitz defects", local.wigner_seitz),
    "structure-factor": ("Static structure factor", local.static_structure_factor),
    "band-gap": ("Band gap and band edges", electronic.band_gap),
    "effective-mass": ("Effective-mass tensor", electronic.effective_mass),
    "work-function": ("Work function", electronic.work_function),
    "equation-of-state": ("Equation of state", thermal.equation_of_state),
    "elastic-tensor": ("Elastic tensor and moduli", thermal.elastic_tensor),
    "phonon-dos": ("Phonon density of states", thermal.phonon_dos),
    "harmonic-thermodynamics": ("Harmonic thermodynamics", thermal.harmonic_thermodynamics),
    "neb": ("NEB migration path", simulation.migration_path),
    "nvt": ("NVT dynamics", simulation.nvt_dynamics),
    "npt": ("NPT dynamics", simulation.npt_dynamics),
    "reciprocal-path": ("Symmetry reciprocal-space path", simulation.reciprocal_path),
}


def tool_catalog():
    """Return names, documentation and parameter defaults for all twenty tools."""
    catalog = []
    for index, (key, (name, function)) in enumerate(TOOLS.items(), 1):
        parameters = []
        for parameter in inspect.signature(function).parameters.values():
            required = parameter.default is inspect.Parameter.empty
            parameters.append({"name": parameter.name, "required": required,
                               "default": None if required else parameter.default})
        catalog.append({"priority": index, "id": key, "name": name,
                        "description": inspect.getdoc(function), "parameters": parameters})
    return catalog


def _read_frames(path):
    if path.name.upper().startswith("XDATCAR"):
        # VASP can repeat the header/cell between configurations. ASE's
        # XDATCAR reader does not support that valid output variant.
        from pymatgen.io.vasp.outputs import Xdatcar
        from pymatgen.io.ase import AseAtomsAdaptor
        return [AseAtomsAdaptor.get_atoms(frame) for frame in Xdatcar(path).structures]
    from ase.io import read
    return read(str(path), index=":")


def _reference(value, name, base, *, timestep_fs=None):
    import numpy as np
    if not isinstance(value, dict) or "file" not in value:
        return value
    if set(value) - {"file", "field", "column"}:
        raise ValueError("Unknown data-reference options for " + name)
    path = (base / value["file"]).resolve()
    field = value.get("field")
    column = value.get("column")
    if column is not None and (isinstance(column, bool) or not isinstance(column, int) or column < 0):
        raise ValueError("Column must be a nonnegative integer")
    def select_column(result):
        if column is None:
            return result
        table = np.asarray(result)
        if table.ndim != 2 or column >= table.shape[1]:
            raise ValueError("Column is outside the two-dimensional data table")
        return table[:, column]
    if path.suffix.lower() == ".json":
        result = json.loads(path.read_text(encoding="utf-8"))
        if field is not None:
            for component in field.split("."):
                result = result[int(component)] if isinstance(result, list) else result[component]
        return select_column(result)
    if path.name.upper().startswith("EIGENVAL") and name == "energies_eV":
        if field is not None or column is not None:
            raise ValueError("EIGENVAL supplies all spin/k-point/band energies; field/column selection is unsupported")
        from .spectra import read_bands
        bands = read_bands(path)
        return np.asarray(list(bands["energies_eV"].values()))
    if path.suffix.lower() == ".npy":
        if field is not None:
            raise ValueError("NPY arrays do not have named fields")
        result = np.load(path, allow_pickle=False)
    elif path.suffix.lower() in (".csv", ".txt", ".dat"):
        if field is not None:
            raise ValueError("Numeric text tables do not have named fields; select a column instead")
        # Numerical CSV: no guessing about units, header rows or column identity.
        result = np.loadtxt(path, delimiter="," if path.suffix.lower() == ".csv" else None, ndmin=2)
        if column is None and result.shape[1] == 1:
            result = result[:, 0]
    else:
        from ase import units
        if column is not None:
            raise ValueError("Column selection is supported only for JSON/NPY/numeric tables")
        if field is None:
            field = {"velocities": "velocities", "masses": "masses", "cell": "cell",
                     "reference_cell": "cell", "current_cell": "cell"}.get(name, "positions")
        frames = _read_frames(path)
        if not frames:
            raise ValueError("Structure/trajectory file contains no frames")
        if timestep_fs is not None and field in ("positions", "velocities"):
            times = [frame.info.get("time_fs") for frame in frames]
            if any(time is not None for time in times):
                if any(time is None for time in times):
                    raise ValueError("Trajectory has incomplete time_fs metadata")
                intervals = np.diff(np.asarray(times, dtype=float))
                if not np.isfinite(times).all() or np.any(intervals <= 0) or not np.allclose(intervals, timestep_fs, rtol=1e-8, atol=1e-10):
                    raise ValueError("Trajectory sampling must be uniform and match timestep_fs; resample or trim the final partial interval")
        symbols = frames[0].get_chemical_symbols()
        if any(frame.get_chemical_symbols() != symbols for frame in frames):
            raise ValueError("Trajectory atom species/order must remain consistent")
        if field == "velocities":
            if any(not frame.has("momenta") for frame in frames):
                raise ValueError("Trajectory has no velocities; supply explicit A/fs data")
            result = np.asarray([frame.get_velocities()*units.fs for frame in frames])
        elif field == "cell":
            if any(not np.allclose(frame.cell, frames[0].cell, rtol=0, atol=1e-10) for frame in frames):
                raise ValueError("This analysis requires a fixed cell")
            result = frames[0].cell.array
        elif field == "masses":
            result = frames[0].get_masses()
        elif field in (None, "positions"):
            result = np.asarray([frame.positions for frame in frames])
            if len(frames) == 1:
                result = result[0]
        else:
            raise ValueError("Structure field must be positions, velocities, cell or masses")
    return select_column(result)


def _structure(value, base):
    from .._structure import Structure
    if isinstance(value, Structure):
        return value
    from ase import Atoms
    from ase.io import read
    from .simulation import from_ase
    if isinstance(value, dict) and set(value) == {"file"}:
        path = (base / value["file"]).resolve()
        if path.suffix.lower() == ".json":
            return _structure(json.loads(path.read_text(encoding="utf-8")), path.parent)
        return from_ase(read(str(path)))
    if not isinstance(value, dict) or set(value) - {"symbols", "positions", "cell"}:
        raise ValueError("Structure requires a file reference or symbols, positions and optional cell")
    return from_ase(Atoms(symbols=value["symbols"], positions=value["positions"], cell=value.get("cell"),
                         pbc=value.get("cell") is not None))


def run_tool(tool, parameters, *, base_directory="."):
    """Execute one catalog tool; preserve numerical arrays in the returned result."""
    if tool not in TOOLS:
        raise ValueError("Unknown scientific tool: " + str(tool))
    if not isinstance(parameters, dict):
        raise ValueError("parameters must be an object")
    function = TOOLS[tool][1]
    # Validate keys before importing calculators or reading data.
    inspect.signature(function).bind(**parameters)
    base = Path(base_directory).resolve()
    resolved = {}
    for name, value in parameters.items():
        if name in ("structure", "initial", "final"):
            resolved[name] = _structure(value, base)
        elif name in ("calculator", "calculator_factory"):
            if callable(value) or not isinstance(value, dict):
                resolved[name] = value
                continue
            if set(value) - {"module", "attribute", "kwargs"} or not {"module", "attribute"} <= set(value):
                raise ValueError("Calculator requires module, attribute and optional kwargs")
            from .simulation import load_calculator
            import copy
            configuration = copy.deepcopy(value)
            image_index = [0]
            def factory(config=configuration, is_factory=name == "calculator_factory", counter=image_index):
                kwargs = copy.deepcopy(config.get("kwargs", {}))
                if is_factory:
                    if "directory" in kwargs:
                        kwargs["directory"] = str((base / kwargs["directory"]).resolve() / ("image-" + str(counter[0])))
                    counter[0] += 1
                return load_calculator(config["module"], config["attribute"], **kwargs)
            resolved[name] = factory if name == "calculator_factory" else factory()
        else:
            timestep = parameters.get("timestep_fs") if tool in ("msd", "vacf", "vibrational-spectrum") else None
            resolved[name] = _reference(value, name, base, timestep_fs=timestep)
    return function(**resolved)


def json_result(value):
    """Strict JSON representation: invalid local environments become null, with validity masks retained."""
    from .._structure import Structure
    if isinstance(value, Structure):
        return {"symbols": [a.symbol for a in value.atoms], "positions": [[a.x, a.y, a.z] for a in value.atoms],
                "cell": value.cell}
    if isinstance(value, dict):
        return {str(key): json_result(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_result(item) for item in value]
    if hasattr(value, "tolist"):
        return json_result(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def result_report(tool, result):
    """Compact desktop summary; full arrays and structures remain in result JSON."""
    import numpy as np
    from .._structure import Structure
    lines = [TOOLS[tool][0], "", "Full numerical data are available through Save results.", ""]
    def describe(key, value):
        if isinstance(value, Structure):
            lines.append(key + ": " + str(len(value)) + " atoms")
        elif isinstance(value, dict):
            for child, item in value.items():
                describe(key + "." + child if key else child, item)
        elif isinstance(value, (list, tuple)) and value and isinstance(value[0], Structure):
            lines.append(key + ": " + str(len(value)) + " structural frames")
        elif isinstance(value, (list, tuple, np.ndarray)):
            data = np.asarray(value)
            if data.dtype.kind in "biuf" and data.size:
                finite = data[np.isfinite(data)]
                lines.append(key + " (shape " + str(data.shape) + "):")
                if data.size <= 36:
                    lines.append(np.array2string(data, precision=7, max_line_width=90))
                elif finite.size:
                    lines.append("  min={:.7g}, max={:.7g}, mean={:.7g}; {} valid / {} values".format(
                        finite.min(), finite.max(), finite.mean(), finite.size, data.size))
                else:
                    lines.append("  No valid values")
            else:
                lines.append(key + ": " + str(json_result(value))[:1000])
        else:
            lines.append(key + ": " + str(value))
    describe("", result)
    return "\n".join(lines) + "\n"


def export_structures(result, path, *, overwrite=False):
    """Export returned simulation frames/images or a primitive structure as extXYZ.

    Returns False for numerical analyses with no structural output. Velocities
    from NVT/NPT are retained as momenta through ASE's extXYZ convention.
    """
    from ase import units
    from ase.io import write
    from .simulation import to_ase
    frames = result.get("frames", result.get("images"))
    if frames is None and "primitive_structure" in result:
        frames = [result["primitive_structure"]]
    if frames is None:
        return False
    destination = Path(path)
    if destination.exists() and not overwrite:
        raise FileExistsError("Structure output exists; choose another file or use --overwrite")
    atoms = [to_ase(frame) for frame in frames]
    velocities = result.get("velocities_A_per_fs")
    for index, frame in enumerate(atoms):
        if velocities is not None:
            frame.set_velocities(velocities[index] / units.fs)
        if "time_fs" in result:
            frame.info["time_fs"] = float(result["time_fs"][index])
    destination.parent.mkdir(parents=True, exist_ok=True)
    import os
    import tempfile
    descriptor, temporary = tempfile.mkstemp(prefix=".atomforge-frames-", suffix=".xyz", dir=destination.parent)
    os.close(descriptor)
    try:
        write(temporary, atoms, format="extxyz")
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return True
