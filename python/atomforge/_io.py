"""File I/O: native XYZ/extxyz, and optional ASE fallback for CIF/VASP/PDB/etc."""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING, Tuple

if TYPE_CHECKING:
    from ._structure import Structure

# ── CPK / Jmol colour table (RGB 0-1) ──────────────────────────────────────────

_COLORS: dict[str, Tuple[float, float, float]] = {
    "H":  (1.000, 1.000, 1.000), "He": (0.851, 1.000, 1.000),
    "Li": (0.800, 0.502, 1.000), "Be": (0.761, 1.000, 0.000),
    "B":  (1.000, 0.710, 0.710), "C":  (0.565, 0.565, 0.565),
    "N":  (0.188, 0.314, 0.973), "O":  (1.000, 0.051, 0.051),
    "F":  (0.565, 0.878, 0.314), "Ne": (0.702, 0.890, 0.961),
    "Na": (0.671, 0.361, 0.949), "Mg": (0.541, 1.000, 0.000),
    "Al": (0.749, 0.651, 0.651), "Si": (0.941, 0.784, 0.627),
    "P":  (1.000, 0.502, 0.000), "S":  (1.000, 1.000, 0.188),
    "Cl": (0.122, 0.941, 0.122), "Ar": (0.502, 0.820, 0.890),
    "K":  (0.561, 0.251, 0.831), "Ca": (0.239, 1.000, 0.000),
    "Sc": (0.902, 0.902, 0.902), "Ti": (0.749, 0.761, 0.780),
    "V":  (0.651, 0.651, 0.671), "Cr": (0.541, 0.600, 0.780),
    "Mn": (0.612, 0.478, 0.780), "Fe": (0.878, 0.400, 0.200),
    "Co": (0.941, 0.565, 0.627), "Ni": (0.314, 0.816, 0.314),
    "Cu": (0.784, 0.502, 0.200), "Zn": (0.490, 0.502, 0.690),
    "Ga": (0.761, 0.561, 0.561), "Ge": (0.400, 0.561, 0.561),
    "As": (0.741, 0.502, 0.890), "Se": (1.000, 0.631, 0.000),
    "Br": (0.651, 0.161, 0.161), "Kr": (0.361, 0.722, 0.820),
    "Rb": (0.439, 0.180, 0.690), "Sr": (0.000, 1.000, 0.000),
    "Y":  (0.580, 1.000, 1.000), "Zr": (0.580, 0.878, 0.878),
    "Nb": (0.451, 0.761, 0.788), "Mo": (0.329, 0.710, 0.710),
    "Tc": (0.231, 0.620, 0.620), "Ru": (0.141, 0.561, 0.561),
    "Rh": (0.039, 0.490, 0.549), "Pd": (0.000, 0.412, 0.522),
    "Ag": (0.753, 0.753, 0.753), "Cd": (1.000, 0.851, 0.561),
    "In": (0.651, 0.459, 0.451), "Sn": (0.400, 0.502, 0.502),
    "Sb": (0.620, 0.388, 0.710), "Te": (0.831, 0.478, 0.000),
    "I":  (0.580, 0.000, 0.580), "Xe": (0.259, 0.620, 0.690),
    "Cs": (0.341, 0.090, 0.561), "Ba": (0.000, 0.788, 0.000),
    "La": (0.439, 0.831, 1.000), "Ce": (1.000, 1.000, 0.780),
    "Pr": (0.851, 1.000, 0.780), "Nd": (0.780, 1.000, 0.780),
    "Pm": (0.639, 1.000, 0.780), "Sm": (0.561, 1.000, 0.780),
    "Eu": (0.380, 1.000, 0.780), "Gd": (0.271, 1.000, 0.780),
    "Tb": (0.188, 1.000, 0.780), "Dy": (0.122, 1.000, 0.780),
    "Ho": (0.000, 1.000, 0.612), "Er": (0.000, 0.902, 0.459),
    "Tm": (0.000, 0.831, 0.322), "Yb": (0.000, 0.749, 0.220),
    "Lu": (0.000, 0.671, 0.141), "Hf": (0.302, 0.761, 1.000),
    "Ta": (0.302, 0.651, 1.000), "W":  (0.129, 0.580, 0.839),
    "Re": (0.149, 0.490, 0.671), "Os": (0.149, 0.400, 0.588),
    "Ir": (0.090, 0.329, 0.529), "Pt": (0.816, 0.816, 0.878),
    "Au": (1.000, 0.820, 0.137), "Hg": (0.722, 0.722, 0.816),
    "Tl": (0.651, 0.329, 0.302), "Pb": (0.341, 0.349, 0.380),
    "Bi": (0.620, 0.310, 0.710), "Po": (0.671, 0.361, 0.000),
    "At": (0.459, 0.310, 0.271), "Rn": (0.259, 0.510, 0.588),
    "Fr": (0.259, 0.000, 0.400), "Ra": (0.000, 0.490, 0.000),
    "Ac": (0.439, 0.671, 0.980), "Th": (0.000, 0.729, 1.000),
    "Pa": (0.000, 0.631, 1.000), "U":  (0.000, 0.561, 1.000),
    "Np": (0.000, 0.502, 1.000), "Pu": (0.000, 0.420, 1.000),
}


def _default_color(symbol: str) -> Tuple[float, float, float]:
    return _COLORS.get(symbol, (0.8, 0.8, 0.8))


# ── XYZ / extended-XYZ ─────────────────────────────────────────────────────────

def _load_xyz(path: str) -> "Structure":
    from ._structure import Atom, Structure
    s = Structure()
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()
    if not lines:
        return s
    n = int(lines[0].strip())
    comment = lines[1] if len(lines) > 1 else ""

    # Parse extended-XYZ lattice from comment line
    if 'Lattice=' in comment or 'lattice=' in comment:
        import re
        m = re.search(r'[Ll]attice="([^"]+)"', comment)
        if m:
            vals = list(map(float, m.group(1).split()))
            if len(vals) == 9:
                s.cell = [vals[0:3], vals[3:6], vals[6:9]]

    atom_lines = lines[2:2 + n]
    properties_order = ["species", "pos", "pos", "pos"]
    if 'Properties=' in comment or 'properties=' in comment:
        import re
        m = re.search(r'[Pp]roperties=(\S+)', comment)
        if m:
            properties_order = m.group(1).split(":")

    for line in atom_lines:
        parts = line.split()
        if not parts:
            continue
        sym = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        cr, cg, cb = _default_color(sym)
        s.atoms.append(Atom(sym, x, y, z, cr, cg, cb))
    return s


def _save_xyz(s: "Structure", path: str) -> None:
    lines = [str(len(s.atoms)) + "\n"]
    if s.cell:
        a, b, c = s.cell
        lat = " ".join(f"{v:.6f}" for row in [a, b, c] for v in row)
        lines.append(f'Lattice="{lat}" Properties=species:S:1:pos:R:3\n')
    else:
        lines.append("AtomForge structure\n")
    for atom in s.atoms:
        lines.append(f"{atom.symbol} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)


# ── ASE fallback ────────────────────────────────────────────────────────────────

def _load_ase(path: str) -> "Structure":
    try:
        import ase.io
        import numpy as np
    except ImportError:
        raise ImportError(
            "Install the 'ase' package to load this file format:\n"
            "    pip install ase"
        )
    from ._structure import Atom, Structure
    atoms_ase = ase.io.read(path)
    s = Structure()
    if atoms_ase.cell is not None and any(atoms_ase.pbc):
        cell = np.array(atoms_ase.cell)
        s.cell = [list(cell[0]), list(cell[1]), list(cell[2])]
    for sym, (x, y, z) in zip(atoms_ase.get_chemical_symbols(),
                               atoms_ase.get_positions()):
        cr, cg, cb = _default_color(sym)
        s.atoms.append(Atom(sym, float(x), float(y), float(z), cr, cg, cb))
    return s


def _save_ase(s: "Structure", path: str) -> None:
    try:
        import ase
        import ase.io
        import numpy as np
    except ImportError:
        raise ImportError(
            "Install the 'ase' package to save this file format:\n"
            "    pip install ase"
        )
    syms = [a.symbol for a in s.atoms]
    pos  = np.array([[a.x, a.y, a.z] for a in s.atoms])
    cell = np.array(s.cell) if s.cell else None
    pbc  = s.cell is not None
    ase_atoms = ase.Atoms(symbols=syms, positions=pos, cell=cell, pbc=pbc)
    ase.io.write(path, ase_atoms)


# ── Format dispatch ─────────────────────────────────────────────────────────────

_XYZ_EXTS = {".xyz", ".extxyz"}
_ASE_EXTS = {".cif", ".vasp", ".poscar", ".pdb", ".mol", ".mol2", ".lmp",
             ".lammps", ".dump", ".gen", ".sdf", ".cfg"}


def load(path: str) -> "Structure":
    """Load a structure from file.  XYZ/extxyz is handled natively; other formats require ASE."""
    ext = Path(path).suffix.lower()
    if ext in _XYZ_EXTS:
        return _load_xyz(path)
    if ext in _ASE_EXTS or not ext:
        return _load_ase(path)
    # Try XYZ first, fall back to ASE
    try:
        return _load_xyz(path)
    except Exception:
        return _load_ase(path)


def save(s: "Structure", path: str) -> None:
    """Save a structure to file.  XYZ/extxyz is handled natively; other formats require ASE."""
    ext = Path(path).suffix.lower()
    if ext in _XYZ_EXTS:
        _save_xyz(s, path)
    elif ext in _ASE_EXTS:
        _save_ase(s, path)
    else:
        # Default to XYZ
        _save_xyz(s, path)
