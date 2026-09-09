"""Format dispatch and compatibility imports for native structure I/O.

Add a codec under _formats and register its extensions in _FORMAT_MAP.
Existing private helper imports remain available for downstream callers.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Callable, Dict, Tuple

from ._elements import _COLORS, _MASSES, _default_color, _clean_symbol
from ._cell import _cell_from_params, _frac_to_cart, _cart_to_frac
from ._formats.xyz import _load_xyz, _save_xyz
from ._formats.vasp import _load_vasp, _save_vasp
from ._formats.pdb import _load_pdb, _save_pdb
from ._formats.cif import _load_cif, _save_cif
from ._formats.lammps import _load_lammps, _save_lammps

if TYPE_CHECKING:
    from ._structure import Structure

_FORMAT_MAP: Dict[str, Tuple[Callable, Callable]] = {
    ".xyz":     (_load_xyz,    _save_xyz),
    ".extxyz":  (_load_xyz,    _save_xyz),
    ".vasp":    (_load_vasp,   _save_vasp),
    ".poscar":  (_load_vasp,   _save_vasp),
    ".contcar": (_load_vasp,   _save_vasp),
    ".pdb":     (_load_pdb,    _save_pdb),
    ".ent":     (_load_pdb,    _save_pdb),
    ".cif":     (_load_cif,    _save_cif),
    ".lmp":     (_load_lammps, _save_lammps),
    ".lammps":  (_load_lammps, _save_lammps),
    ".data":    (_load_lammps, _save_lammps),
    ".dump":    (_load_lammps, _save_lammps),
}


def _resolve_format(path: str) -> Tuple[Callable, Callable]:
    """Resolve both directions together so read/write filename rules agree."""
    filename = Path(path)
    extension = filename.suffix.lower()
    if not extension and filename.name.upper() in ("POSCAR", "CONTCAR"):
        return _load_vasp, _save_vasp
    return _FORMAT_MAP.get(extension, (_load_xyz, _save_xyz))


def load(path: str) -> "Structure":
    """Load a structure; format is inferred from the filename."""
    loader, _ = _resolve_format(path)
    return loader(path)


def save(s: "Structure", path: str) -> None:
    """Save a structure; format is inferred from the filename."""
    _, saver = _resolve_format(path)
    saver(s, path)
