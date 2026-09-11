"""Headless access to the builders in an installed AtomForge executable."""

from dataclasses import dataclass
import os
from pathlib import Path
import subprocess
import tempfile

from ._io import load
from ._structure import Structure
from ._viewer import _find_atomforge


_MODES = ("bulk", "gb", "poly", "nano", "amorphous", "sss", "dislocation", "custom")
_RESERVED = {"--build", "--help", "--output", "--input"}


@dataclass
class BuildResult:
    """Built structure and native diagnostics; output is a saved path or None.

    Inspect stdout and stderr for algorithm warnings. Explicit output files and
    native orientation sidecars are retained. Structure holds coordinates and
    cell data; it does not import the sidecar's grain metadata.
    """

    structure: Structure
    stdout: str
    stderr: str
    output: object = None


def _executable(executable):
    path = os.fspath(executable) if executable is not None else _find_atomforge()
    if not path:
        raise FileNotFoundError("Install AtomForge and set ATOMFORGE_PATH to its executable")
    return path


def builder_help(mode, *, executable=None):
    """Return installed native help for one of the eight builder modes."""
    if mode not in _MODES:
        raise ValueError("Unknown builder mode: {!r}".format(mode))
    return subprocess.run(
        [_executable(executable), "--help", mode], check=True, capture_output=True,
        text=True, encoding="utf-8", errors="replace", timeout=30,
        creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
    ).stdout


def build(mode, options=(), *, source=None, output=None, executable=None, timeout=300):
    """Run a native builder and return BuildResult without opening the GUI.

    mode is bulk, gb, poly, nano, amorphous, sss, dislocation, or custom.
    options is an argument sequence, e.g. ['--a', '3.61', '--atom', 'Cu 0 0 0'].
    Repeated flags are retained; multi-component values occupy one argument.
    source accepts a Structure or an input path. output optionally retains a
    CIF, XYZ/extXYZ, VASP, PDB, or LAMMPS result and native sidecars; otherwise
    temporary VASP files are cleaned up after loading. Coordinates are Angstrom.
    Native failures raise subprocess.CalledProcessError with stdout/stderr;
    timeout raises subprocess.TimeoutExpired. No command shell is used.
    """
    if mode not in _MODES:
        raise ValueError("Unknown builder mode: {!r}".format(mode))
    if isinstance(options, (str, bytes)):
        raise TypeError("options must be an argument sequence, not a command string")
    arguments = [str(value) for value in options]
    if any(arg in _RESERVED for arg in arguments):
        raise ValueError("Use mode, source, and output arguments instead of reserved CLI flags")
    if timeout is not None and timeout <= 0:
        raise ValueError("timeout must be positive or None")
    exe = _executable(executable)
    destination = Path(output).expanduser().resolve() if output is not None else None
    if destination is not None and destination.suffix.lower() not in (
        ".cif", ".xyz", ".extxyz", ".vasp", ".poscar", ".contcar", ".pdb", ".lmp", ".lammps"
    ):
        raise ValueError("output must use a format readable by the Python package")
    with tempfile.TemporaryDirectory(prefix="atomforge_build_") as directory:
        root = Path(directory)
        command = [exe, "--build", mode]
        if source is not None:
            if isinstance(source, Structure):
                input_path = root / ("source.vasp" if source.cell else "source.xyz")
                source.save(str(input_path))
            else:
                input_path = Path(source).expanduser().resolve()
                if not input_path.is_file():
                    raise FileNotFoundError(str(input_path))
            command += ["--input", str(input_path)]
        target = destination or root / "result.vasp"
        command += arguments + ["--output", str(target)]
        result = subprocess.run(
            command, check=True, capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=timeout,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
        )
        structure = load(str(target))
        return BuildResult(structure, result.stdout, result.stderr, destination)


def wulff(source, facets, *, radius=20, vacuum=5, output=None, executable=None, timeout=300):
    """Build a Wulff particle using (h, k, l, energy) facet families.

    Positive relative surface energies set plane distances; radius is the
    maximum plane distance in Angstrom, not the farthest vertex radius.
    Uses symmetry expansion and the same native engine as the desktop.
    Returns BuildResult. Requires an executable supporting --shape wulff.
    """
    options = ["--shape", "wulff", "--radius", str(radius), "--vacuum", str(vacuum)]
    rows = list(facets)
    if not rows:
        raise ValueError("Provide at least one facet family")
    for row in rows:
        row = tuple(row)
        if len(row) != 4:
            raise ValueError("Each facet must contain h, k, l, and energy")
        options += ["--facet", " ".join(str(value) for value in row)]
    return build("nano", options, source=source, output=output,
                 executable=executable, timeout=timeout)
