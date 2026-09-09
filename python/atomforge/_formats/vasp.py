"""Native VASP structure reader and writer."""
from __future__ import annotations

from typing import TYPE_CHECKING

from .._elements import _default_color, _clean_symbol
from .._cell import _frac_to_cart, _cart_to_frac

if TYPE_CHECKING:
    from .._structure import Structure


def _load_vasp(path: str) -> "Structure":
    from .._structure import Atom, Structure
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [l.rstrip("\n") for l in fh]

    if len(lines) < 8:
        raise ValueError("invalid POSCAR: file is incomplete")

    scale = float(lines[1].split()[0])
    raw_cell = []
    for i in range(2, 5):
        row = [float(v) for v in lines[i].split()[:3]]
        if len(row) != 3:
            raise ValueError("invalid POSCAR: each lattice vector needs three values")
        raw_cell.append(row)

    if scale == 0.0:
        raise ValueError("invalid POSCAR: scale factor cannot be zero")
    if scale < 0.0:
        # A negative VASP scale specifies the desired cell volume.
        a, b, c = raw_cell
        volume = abs(
            a[0] * (b[1]*c[2] - b[2]*c[1])
            - a[1] * (b[0]*c[2] - b[2]*c[0])
            + a[2] * (b[0]*c[1] - b[1]*c[0])
        )
        if volume <= 1e-15:
            raise ValueError("invalid POSCAR: lattice vectors are singular")
        scale = (-scale / volume) ** (1.0 / 3.0)
    cell = [[value * scale for value in row] for row in raw_cell]

    # VASP5: species names on line 5, counts on line 6
    # VASP4: counts on line 5 (no species names)
    tok5 = lines[5].split()
    tok6 = lines[6].split() if len(lines) > 6 else []
    if tok5 and not tok5[0][0].isdigit():
        species = tok5
        counts  = [int(x) for x in tok6]
        coord_line = 7
    else:
        species = [f"X{i+1}" for i in range(len(tok5))]
        counts  = [int(x) for x in tok5]
        coord_line = 6

    if len(species) != len(counts) or not counts or any(n < 0 for n in counts):
        raise ValueError("invalid POSCAR: species and non-negative counts must match")
    if coord_line >= len(lines):
        raise ValueError("invalid POSCAR: missing coordinate mode")
    if lines[coord_line].strip().lower().startswith("s"):
        coord_line += 1
    if coord_line >= len(lines):
        raise ValueError("invalid POSCAR: missing coordinate mode")
    mode_line = lines[coord_line].strip().lower()
    if not (mode_line.startswith("d") or mode_line.startswith("c") or
            mode_line.startswith("k")):
        raise ValueError(f"invalid POSCAR coordinate mode: {lines[coord_line]!r}")
    direct = mode_line.startswith("d")
    coord_line += 1

    s = Structure()
    s.cell = cell
    for sym, count in zip(species, counts):
        sym = _clean_symbol(sym)
        for _ in range(count):
            if coord_line >= len(lines):
                raise ValueError("invalid POSCAR: fewer atom coordinates than declared")
            row = lines[coord_line].split()
            coord_line += 1
            if len(row) < 3:
                raise ValueError("invalid POSCAR: atom coordinate needs three values")
            fx, fy, fz = float(row[0]), float(row[1]), float(row[2])
            if direct:
                x, y, z = _frac_to_cart(fx, fy, fz, cell)
            else:
                x, y, z = fx * scale, fy * scale, fz * scale
            cr, cg, cb = _default_color(sym)
            s.atoms.append(Atom(sym, x, y, z, cr, cg, cb))
    return s


def _save_vasp(s: "Structure", path: str) -> None:
    # Collect species order (preserving first appearance)
    seen: list[str] = []
    for a in s.atoms:
        if a.symbol not in seen:
            seen.append(a.symbol)

    groups = {sym: [] for sym in seen}
    for a in s.atoms:
        groups[a.symbol].append(a)

    cell = s.cell or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    lines = ["AtomForge structure\n", "1.0\n"]
    for row in cell:
        lines.append(f"  {row[0]:16.10f}  {row[1]:16.10f}  {row[2]:16.10f}\n")
    lines.append("  " + "  ".join(seen) + "\n")
    lines.append("  " + "  ".join(str(len(groups[sym])) for sym in seen) + "\n")
    lines.append("Direct\n")
    for sym in seen:
        for a in groups[sym]:
            if s.cell:
                fx, fy, fz = _cart_to_frac(a.x, a.y, a.z, cell)
            else:
                fx, fy, fz = a.x, a.y, a.z
            lines.append(f"  {fx:16.10f}  {fy:16.10f}  {fz:16.10f}\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
