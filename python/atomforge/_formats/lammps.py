"""Native LAMMPS structure reader and writer."""
from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

from .._elements import _default_color, _clean_symbol, _MASSES

if TYPE_CHECKING:
    from .._structure import Structure


def _load_lammps(path: str) -> "Structure":
    """
    Load a LAMMPS data file (atom_style atomic or full).
    Element names are read from '# elem1 elem2 ...' comment on the Masses line
    or inferred from atomic mass if the comment is absent.
    """
    from .._structure import Atom, Structure

    # Common atomic masses → symbol
    _MASS_TO_SYM = {
        1.008: "H",   4.003: "He",  6.941: "Li",  9.012: "Be", 10.811: "B",
        12.011: "C",  14.007: "N",  15.999: "O",  18.998: "F",  20.180: "Ne",
        22.990: "Na", 24.305: "Mg", 26.982: "Al", 28.086: "Si", 30.974: "P",
        32.065: "S",  35.453: "Cl", 39.948: "Ar", 39.098: "K",  40.078: "Ca",
        47.867: "Ti", 51.996: "Cr", 54.938: "Mn", 55.845: "Fe", 58.933: "Co",
        58.693: "Ni", 63.546: "Cu", 65.38:  "Zn", 107.868: "Ag",183.84: "W",
        195.084: "Pt", 196.967: "Au",
    }

    def _closest_mass(m: float):
        """The element whose standard mass is within 0.5 u, else None."""
        best, sym = 1e9, None
        for mass, s in _MASS_TO_SYM.items():
            if abs(m - mass) < best:
                best, sym = abs(m - mass), s
        return sym if best < 0.5 else None

    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    type_to_sym: dict[int, str] = {}
    xlo, xhi, ylo, yhi, zlo, zhi = 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
    xy, xz, yz = 0.0, 0.0, 0.0
    atoms: list[tuple[int, float, float, float]] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if "xlo xhi" in line:
            parts = line.split()
            xlo, xhi = float(parts[0]), float(parts[1])
        elif "ylo yhi" in line:
            parts = line.split()
            ylo, yhi = float(parts[0]), float(parts[1])
        elif "zlo zhi" in line:
            parts = line.split()
            zlo, zhi = float(parts[0]), float(parts[1])
        elif "xy xz yz" in line:
            parts = line.split()
            xy, xz, yz = float(parts[0]), float(parts[1]), float(parts[2])
        elif line.startswith("Masses"):
            i += 1
            while i < len(lines):
                mline = lines[i].strip()
                i += 1
                if not mline or mline.startswith("#"):
                    continue  # skip blank separator lines
                if mline[0].isalpha():
                    i -= 1   # put the section header back
                    break
                parts = mline.split()
                if not parts or not parts[0].isdigit():
                    break
                tid = int(parts[0])
                mass = float(parts[1])
                cm = re.search(r"#\s*([A-Za-z]+)", mline)
                sym = _clean_symbol(cm.group(1)) if cm else _closest_mass(mass)
                if sym:
                    type_to_sym[tid] = sym
            continue

        elif line.startswith("Atoms"):
            style = line.partition("#")[2].strip().split()
            style = style[0].lower() if style else None
            if style not in (None, "atomic", "full", "charge", "molecular", "bond", "angle"):
                raise ValueError(f"unsupported LAMMPS atom style: {style}")
            i += 1
            while i < len(lines):
                aline = lines[i].strip()
                i += 1
                if not aline:
                    continue  # skip blank separator lines
                if aline[0].isalpha():
                    i -= 1
                    break
                parts = aline.partition("#")[0].split()
                # atomic: id type x y z;  charge: id type q x y z
                # full:   id mol type q x y z;  molecular: id mol type x y z
                # (image flags may follow). Without a style comment the
                # column count decides, as in the desktop.
                row_style = style or {5: "atomic", 8: "atomic", 6: "charge", 9: "charge",
                                      7: "full", 10: "full"}.get(len(parts))
                if row_style is None:
                    raise ValueError(f"cannot tell the LAMMPS atom style of {aline!r}; "
                                     "name it after the section, e.g. 'Atoms # atomic'")
                type_col, x_col = {"atomic": (1, 2), "charge": (1, 3), "full": (2, 4)}.get(row_style, (2, 3))
                try:
                    atoms.append((int(parts[0]), int(parts[type_col]),
                                  float(parts[x_col]), float(parts[x_col + 1]), float(parts[x_col + 2])))
                except (ValueError, IndexError) as error:
                    raise ValueError(f"cannot read LAMMPS atom line {aline!r} as {row_style} style") from error
            continue
        i += 1

    s = Structure()
    s.cell = [
        [xhi - xlo,  0.0,       0.0],
        [xy,         yhi - ylo, 0.0],
        [xz,         yz,        zhi - zlo],
    ]
    unknown = sorted({tid for _, tid, *_ in atoms if tid not in type_to_sym})
    if unknown:
        raise ValueError(f"the element of LAMMPS atom type {', '.join(map(str, unknown))} is unknown: "
                         "give its mass, or its symbol after the mass ('1 63.546 # Cu')")
    atoms.sort(key=lambda row: row[0])  # LAMMPS order is by atom id
    for _, tid, x, y, z in atoms:
        sym = type_to_sym[tid]
        cr, cg, cb = _default_color(sym)
        s.atoms.append(Atom(sym, x - xlo, y - ylo, z - zlo, cr, cg, cb))
    return s


def _save_lammps(s: "Structure", path: str) -> None:
    seen: list[str] = []
    for a in s.atoms:
        if a.symbol not in seen:
            seen.append(a.symbol)
    sym_to_tid = {sym: i + 1 for i, sym in enumerate(seen)}

    if s.cell:
        a, b, c = s.cell
        # Express the cell and atoms in the restricted triclinic basis.
        lx = math.sqrt(sum(v*v for v in a))
        if lx <= 1e-12:
            raise ValueError("LAMMPS export requires a non-degenerate cell")
        ex = [v / lx for v in a]
        xy = sum(b[i]*ex[i] for i in range(3))
        by = [b[i] - xy*ex[i] for i in range(3)]
        ly = math.sqrt(sum(v*v for v in by))
        if ly <= 1e-12:
            raise ValueError("LAMMPS export requires a non-degenerate cell")
        ey = [v / ly for v in by]
        ez = [ex[1]*ey[2] - ex[2]*ey[1], ex[2]*ey[0] - ex[0]*ey[2],
              ex[0]*ey[1] - ex[1]*ey[0]]
        xz = sum(c[i]*ex[i] for i in range(3))
        yz = sum(c[i]*ey[i] for i in range(3))
        lz = sum(c[i]*ez[i] for i in range(3))
        if not all(math.isfinite(v) for v in (lx, ly, lz, xy, xz, yz)) or lz <= 1e-12:
            raise ValueError("LAMMPS export requires a finite right-handed cell")
        positions = [tuple(sum(v*d for v, d in zip((atom.x, atom.y, atom.z), axis))
                           for axis in (ex, ey, ez)) for atom in s.atoms]
        xlo, xhi, ylo, yhi, zlo, zhi = 0.0, lx, 0.0, ly, 0.0, lz
    else:
        positions = [(a.x, a.y, a.z) for a in s.atoms]
        bounds = [(min((p[i] for p in positions), default=0.0) - 1.0,
                   max((p[i] for p in positions), default=0.0) + 1.0)
                  for i in range(3)]
        (xlo, xhi), (ylo, yhi), (zlo, zhi) = bounds
        xy = xz = yz = 0.0

    lines = [
        "AtomForge structure\n\n",
        f"{len(s.atoms)} atoms\n",
        f"{len(seen)} atom types\n\n",
        f"{xlo:.10f} {xhi:.10f} xlo xhi\n",
        f"{ylo:.10f} {yhi:.10f} ylo yhi\n",
        f"{zlo:.10f} {zhi:.10f} zlo zhi\n",
        f"{xy:.10f} {xz:.10f} {yz:.10f} xy xz yz\n\n",
        "Masses\n\n",
    ]
    for sym in seen:
        if sym not in _MASSES:
            raise ValueError(f"unknown atomic mass for element {sym!r}")
        mass = _MASSES[sym]
        lines.append(f"  {sym_to_tid[sym]}  {mass:.3f}  # {sym}\n")
    lines.append("\nAtoms  # atomic\n\n")
    for i, (a, (x, y, z)) in enumerate(zip(s.atoms, positions), 1):
        lines.append(f"{i} {sym_to_tid[a.symbol]} {x:.10f} {y:.10f} {z:.10f}\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
