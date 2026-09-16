"""Additional structure formats through the installed desktop I/O engine."""

import os
from pathlib import Path
import subprocess
import tempfile


def convert(source, destination, *, format=None, executable=None, timeout=300):
    """Convert with native Open Babel support without opening a GUI."""
    from ..builders import _executable
    command = [_executable(executable), "--convert", "--input", str(Path(source).resolve()),
               "--output", str(Path(destination).resolve())]
    if format:
        command.extend(("--format", format))
    subprocess.run(command, check=True, capture_output=True, text=True, timeout=timeout,
                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)


def _load_native(path):
    from .xyz import _load_xyz
    extension = Path(path).suffix.lower()
    if extension in (".pwi", ".gjf", ".com"):
        from ase.io import read
        from ..science.simulation import from_ase
        return from_ase(read(str(path), format="espresso-in" if extension==".pwi" else "gaussian-in"))
    with tempfile.TemporaryDirectory(prefix="atomforge_format_") as directory:
        target = Path(directory) / "converted.xyz"
        convert(path, target, format="extxyz")
        return _load_xyz(target)


def _save_native(structure, path):
    from .cif import _save_cif
    from .xyz import _save_xyz
    extension = Path(path).suffix.lower()
    if extension in (".pwi", ".gjf", ".com"):
        from ase.io import write
        from ..science.simulation import to_ase
        atoms = to_ase(structure)
        if extension==".pwi":
            if structure.cell is None:
                raise ValueError("Quantum ESPRESSO input requires a unit cell")
            # Structural template: users must supply actual pseudopotential files.
            write(str(path), atoms, format="espresso-in",
                  pseudopotentials={a.symbol:a.symbol+".UPF" for a in structure.atoms})
        else:
            write(str(path), atoms, format="gaussian-in")
        return
    with tempfile.TemporaryDirectory(prefix="atomforge_format_") as directory:
        source = Path(directory) / ("source.cif" if structure.cell else "source.xyz")
        (_save_cif if structure.cell else _save_xyz)(structure, source)
        convert(source, path)
