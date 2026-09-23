"""Structural analyses using the same numerical engine as the desktop.

Coordinates and distances are Angstrom. Results are lists of column dictionaries;
native errors propagate as CalledProcessError with diagnostic stderr.
"""

import csv
import os
from pathlib import Path
import subprocess
import tempfile

from ._structure import Structure
from .builders import _executable


def analyze(source, method, *, executable=None, timeout=300, **parameters):
    """Run cna, rdf, adf, sro or interstitial and return named CSV rows.

    Parameter names match CLI flags with underscores replacing hyphens.
    ``no_pbc=True`` and ``raw=True`` are switches. SRO currently uses the
    native shell convention and does not accept a no_pbc override.
    """
    if method not in ("cna", "rdf", "adf", "sro", "interstitial"):
        raise ValueError("Unknown analysis method")
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    allowed = {
        "cna": {"cutoff", "scale", "no_pbc"},
        "rdf": {"cutoff", "minimum", "bins", "centre_z", "neighbor_z", "smooth", "no_pbc", "raw"},
        "adf": {"cutoff", "bins", "smooth", "centre", "neighbor_a", "neighbor_b", "no_pbc", "raw"},
        "sro": {"shells", "tolerance"},
        "interstitial": {"resolution", "limit", "clearance", "separation"},
    }
    if set(parameters) - allowed[method]:
        raise ValueError("Unsupported parameters: " + ", ".join(sorted(set(parameters) - allowed[method])))
    with tempfile.TemporaryDirectory(prefix="atomforge_analysis_") as directory:
        root = Path(directory)
        if isinstance(source, Structure):
            source_path = root / ("source.vasp" if source.cell else "source.xyz")
            source.save(source_path)
        else:
            source_path = Path(source).expanduser().resolve()
            if not source_path.is_file():
                raise FileNotFoundError(str(source_path))
        result = root / "result.csv"
        command = [_executable(executable), "--analyze", method, "--input", str(source_path),
                   "--output", str(result)]
        for key, value in parameters.items():
            if key in ("no_pbc", "raw"):
                if not isinstance(value, bool):
                    raise TypeError(key + " must be boolean")
                if value:
                    command.append("--" + key.replace("_", "-"))
            else:
                command.extend(("--" + key.replace("_", "-"), str(value)))
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        with result.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            for key, value in row.items():
                try:
                    row[key] = float(value)
                except ValueError:
                    pass
        return rows


def cna(source, **parameters):
    """Common-neighbour indices and full-environment classification per atom."""
    return analyze(source, "cna", **parameters)


def rdf(source, **parameters):
    """Radial distribution, directed raw pair counts and coordination."""
    return analyze(source, "rdf", **parameters)


def adf(source, **parameters):
    """Angular distribution of unordered neighbour pairs, in degrees."""
    return analyze(source, "adf", **parameters)


def sro(source, **parameters):
    """Warren-Cowley and Rao-Curtin short-range order by neighbour shell."""
    return analyze(source, "sro", **parameters)


def interstitial_sites(source, **parameters):
    """Interstitial candidates, clearance, volume and coordination."""
    return analyze(source, "interstitial", **parameters)


def nye_tensor(deformed, reference, *, cutoff=3.0, no_pbc=False, executable=None, timeout=300):
    """Nye (dislocation density) tensor per atom, comparing `deformed` against `reference`.

    Both structures need the same atom count and ordering (atom n in one is
    the same physical atom as atom n in the other), e.g. a perfect crystal
    and that same crystal after inserting a dislocation. Returns one row per
    atom: index, symbol, the nine alpha_ij tensor components (1/Angstrom),
    and their Frobenius norm. See algorithms/NyeTensor.h for the exact
    algorithm (a native reimplementation of BABEL's documented Hartley &
    Mishin lattice-correspondence method).
    """
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    with tempfile.TemporaryDirectory(prefix="atomforge_nye_") as directory:
        root = Path(directory)

        def resolve(source, name):
            if isinstance(source, Structure):
                path = root / (name + (".vasp" if source.cell else ".xyz"))
                source.save(path)
                return path
            path = Path(source).expanduser().resolve()
            if not path.is_file():
                raise FileNotFoundError(str(path))
            return path

        deformed_path = resolve(deformed, "deformed")
        reference_path = resolve(reference, "reference")
        result = root / "result.csv"
        command = [_executable(executable), "--analyze", "nye",
                   "--input", str(deformed_path), "--reference", str(reference_path),
                   "--cutoff", str(cutoff), "--output", str(result)]
        if no_pbc:
            command.append("--no-pbc")
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        with result.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            for key, value in row.items():
                if key in ("index", "symbol"):
                    continue
                row[key] = float(value)
            row["index"] = int(row["index"])
        return rows


def vitek_map(deformed, reference, line, *, burgers=0.0, cutoff=0.0, no_pbc=False,
              executable=None, timeout=300):
    """Differential-displacement (Vitek) map between `deformed` and `reference`.

    Both structures need the same atom count and ordering. `line` is the
    dislocation line/Burgers direction in Cartesian coordinates (e.g.
    `dislocation.lineDirection`, or (1,1,1) for a <111> bcc screw); `burgers`
    is its magnitude in Angstrom (0 disables wrapping the screw component
    onto the elastic part of the displacement). `cutoff` <= 0 auto-detects
    the first/second-neighbor-shell midpoint. Returns one row per neighbor
    pair: atom indices, both atoms' reference positions, the screw
    (along-line) component, and the two in-plane edge components -- the same
    quantities BABEL's vitek.f90 writes out for gnuplot differential-
    displacement arrows. See algorithms/VitekMap.h.
    """
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    if len(line) != 3:
        raise ValueError("line must have three components")
    with tempfile.TemporaryDirectory(prefix="atomforge_vitek_") as directory:
        root = Path(directory)

        def resolve(source, name):
            if isinstance(source, Structure):
                path = root / (name + (".vasp" if source.cell else ".xyz"))
                source.save(path)
                return path
            path = Path(source).expanduser().resolve()
            if not path.is_file():
                raise FileNotFoundError(str(path))
            return path

        deformed_path = resolve(deformed, "deformed")
        reference_path = resolve(reference, "reference")
        result = root / "result.csv"
        command = [_executable(executable), "--analyze", "vitek",
                   "--input", str(deformed_path), "--reference", str(reference_path),
                   "--line", "{} {} {}".format(*line), "--burgers", str(burgers),
                   "--cutoff", str(cutoff), "--output", str(result)]
        if no_pbc:
            command.append("--no-pbc")
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        with result.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            for key, value in row.items():
                if key in ("index_i", "index_j"):
                    row[key] = int(value)
                else:
                    row[key] = float(value)
        return rows


def build_pattern(reference, *, cutoff=3.2, no_pbc=False, executable=None, timeout=300):
    """Extract a crystallographic neighbor-direction pattern from a perfect structure.

    Returns {"directions": [(dx,dy,dz), ...]} (unit vectors from a
    representative atom to each first-neighbor), for use with
    `detect_pattern`. `cutoff` should enclose exactly the first coordination
    shell. Follows the same two-stage design as BABEL's patternInit
    program, generalized to work with any reference structure rather than a
    fixed bcc/fcc/hcp table. See algorithms/PatternMatch.h.
    """
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    with tempfile.TemporaryDirectory(prefix="atomforge_pattern_") as directory:
        root = Path(directory)
        if isinstance(reference, Structure):
            path = root / ("reference.vasp" if reference.cell else "reference.xyz")
            reference.save(path)
        else:
            path = Path(reference).expanduser().resolve()
            if not path.is_file():
                raise FileNotFoundError(str(path))
        result = root / "pattern.csv"
        command = [_executable(executable), "--analyze", "pattern-init",
                   "--input", str(path), "--cutoff", str(cutoff), "--output", str(result)]
        if no_pbc:
            command.append("--no-pbc")
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        with result.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        directions = [(float(row["dx"]), float(row["dy"]), float(row["dz"])) for row in rows]
        return {"directions": directions, "cutoff": cutoff}


def detect_pattern(structure, pattern, *, cutoff=None, angle_threshold=10.0, no_pbc=False,
                   executable=None, timeout=300):
    """Match each atom's local neighbor environment against `pattern` (from `build_pattern`).

    Returns one row per atom: index, matched (bool), neighbor_count, and
    max_angle_deviation_deg (for matched atoms). Atoms not matching the
    pattern are near a defect (dislocation core, stacking fault, free
    surface, ...). `cutoff` defaults to the cutoff `pattern` was built with.
    """
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    if cutoff is None:
        cutoff = pattern.get("cutoff", 3.2)
    with tempfile.TemporaryDirectory(prefix="atomforge_pattern_") as directory:
        root = Path(directory)
        if isinstance(structure, Structure):
            path = root / ("structure.vasp" if structure.cell else "structure.xyz")
            structure.save(path)
        else:
            path = Path(structure).expanduser().resolve()
            if not path.is_file():
                raise FileNotFoundError(str(path))
        pattern_path = root / "pattern.csv"
        with pattern_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["dx", "dy", "dz"])
            writer.writerows(pattern["directions"])
        result = root / "result.csv"
        command = [_executable(executable), "--analyze", "pattern-detect",
                   "--input", str(path), "--pattern", str(pattern_path),
                   "--cutoff", str(cutoff), "--angle-threshold", str(angle_threshold),
                   "--output", str(result)]
        if no_pbc:
            command.append("--no-pbc")
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        with result.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            row["index"] = int(row["index"])
            row["matched"] = row["matched"] == "1"
            row["neighbor_count"] = int(row["neighbor_count"])
            row["max_angle_deviation_deg"] = float(row["max_angle_deviation_deg"])
        return rows


def fit_dislocation(deformed, reference, line, *, cutoff=3.0, area=1.0, no_pbc=False,
                    executable=None, timeout=300):
    """Recover a dislocation's position and Burgers vector from a measured field.

    Serves the same purpose as BABEL's displacementFit/vitekFit programs
    (self-described in BABEL's own docs as "not documented yet"):
    extracting dislocation properties from data, the reverse of `af.build`
    inserting a known dislocation. Uses the Nye tensor's own defining
    properties directly (its |alpha|-weighted centroid is the best-fit line
    position; its cross-sectional integral is the Burgers vector, by Nye's
    theorem) rather than BABEL's iterative nonlinear least-squares fit.
    `line` is the assumed Cartesian line direction; `area` is the
    cross-sectional area per atom (Angstrom^2) used to scale the Burgers
    vector estimate (leave at 1.0 for an uncalibrated, direction-only
    result). Returns {"line_position": (x,y,z), "burgers_vector": (x,y,z)}.
    See algorithms/DislocationFit.h.
    """
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    if len(line) != 3:
        raise ValueError("line must have three components")
    with tempfile.TemporaryDirectory(prefix="atomforge_fit_") as directory:
        root = Path(directory)

        def resolve(source, name):
            if isinstance(source, Structure):
                path = root / (name + (".vasp" if source.cell else ".xyz"))
                source.save(path)
                return path
            path = Path(source).expanduser().resolve()
            if not path.is_file():
                raise FileNotFoundError(str(path))
            return path

        deformed_path = resolve(deformed, "deformed")
        reference_path = resolve(reference, "reference")
        result = root / "result.csv"
        command = [_executable(executable), "--analyze", "fit-dislocation",
                   "--input", str(deformed_path), "--reference", str(reference_path),
                   "--line", "{} {} {}".format(*line), "--cutoff", str(cutoff),
                   "--area", str(area), "--output", str(result)]
        if no_pbc:
            command.append("--no-pbc")
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        with result.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        values = {row["quantity"]: (float(row["x"]), float(row["y"]), float(row["z"])) for row in rows}
        return values


def prepare_drag(initial, final, *, zeta=0.5, clip_displacement=False, no_pbc=False,
                 executable=None, timeout=300):
    """Prepare a constrained-minimization ("drag") migration-barrier calculation.

    Linearly interpolates atom positions between `initial` and `final`
    (same atom count/order in both) at reaction coordinate `zeta` in
    [0, 1], e.g. to sample points along a dislocation's migration path
    between two neighboring equilibrium (Peierls-valley) positions.
    `clip_displacement` wraps each atom's final-minus-initial vector to its
    minimum-image representation first, for migrations that would otherwise
    appear to cross most of the periodic cell. AtomForge does not perform
    the constrained minimization itself (like BABEL's prepareDrag program);
    this only prepares its inputs. Returns
    (interpolated_structure, constraint_directions), where
    constraint_directions is a list of (dx, dy, dz) per atom -- the
    reaction-coordinate direction an external relaxation/MD code would fix
    the projection onto while relaxing every perpendicular degree of
    freedom. See algorithms/DragPrep.h.
    """
    from ._io import load
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    if not (0.0 <= zeta <= 1.0):
        raise ValueError("zeta must be in [0, 1]")
    with tempfile.TemporaryDirectory(prefix="atomforge_drag_") as directory:
        root = Path(directory)

        def resolve(source, name):
            if isinstance(source, Structure):
                path = root / (name + (".vasp" if source.cell else ".xyz"))
                source.save(path)
                return path
            path = Path(source).expanduser().resolve()
            if not path.is_file():
                raise FileNotFoundError(str(path))
            return path

        initial_path = resolve(initial, "initial")
        final_path = resolve(final, "final")
        output = root / "interpolated.vasp"
        command = [_executable(executable), "--analyze", "drag",
                   "--input", str(initial_path), "--final", str(final_path),
                   "--zeta", str(zeta), "--output", str(output)]
        if clip_displacement:
            command.append("--clip-displacement")
        if no_pbc:
            command.append("--no-pbc")
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        interpolated = load(output)
        constraint_path = Path(str(output) + ".constraint.csv")
        with constraint_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        directions = [(float(row["dx"]), float(row["dy"]), float(row["dz"])) for row in rows]
        return interpolated, directions


def sculpt(source, slabs, *, repeat=(1,1,1), periodic=False, executable=None, timeout=300):
    """Apply intersecting Miller slabs with the desktop sculptor engine.

    Each slab is (h,k,l,lower,upper). Manual bounds are A relative to the
    supercell centroid; periodic bounds mean integer start plane and period
    count, using the native structural-period convention.
    """
    from ._io import load
    if len(repeat)!=3 or any(not isinstance(n,int) or n<1 for n in repeat):
        raise ValueError("Repeat must contain three positive integers")
    rows=[]
    for slab in slabs:
        if len(slab)!=5:
            raise ValueError("Each slab needs h,k,l,lower,upper")
        rows.append(" ".join(str(v) for v in (*slab,int(periodic))))
    if not rows:
        raise ValueError("Provide at least one slab")
    with tempfile.TemporaryDirectory(prefix="atomforge_sculpt_") as directory:
        root=Path(directory)
        if isinstance(source,Structure):
            path=root/"source.vasp"; source.save(path)
        else:
            path=Path(source).resolve()
        output=root/"result.vasp"
        command=[_executable(executable),"--analyze","sculpt","--input",str(path),
                 "--output",str(output),"--slabs",";".join(rows)]
        for flag,value in zip(("--nx","--ny","--nz"),repeat):
            command.extend((flag,str(value)))
        subprocess.run(command,check=True,capture_output=True,text=True,timeout=timeout,
                       creationflags=subprocess.CREATE_NO_WINDOW if os.name=="nt" else 0)
        return load(output)
