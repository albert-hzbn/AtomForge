"""Native XYZ structure reader and writer."""
from __future__ import annotations

import re
from typing import TYPE_CHECKING

from .._elements import _default_color

if TYPE_CHECKING:
    from .._structure import Structure


def _load_xyz(path: str) -> "Structure":
    from .._structure import Atom, Structure
    s = Structure()
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    if not lines:
        return s
    n = int(lines[0].strip())
    if n < 0 or len(lines) < n + 2:
        raise ValueError("invalid XYZ: incomplete file or negative atom count")
    comment = lines[1] if len(lines) > 1 else ""
    m = re.search(r'[Ll]attice="([^"]+)"', comment)
    if m:
        vals = list(map(float, m.group(1).split()))
        if len(vals) == 9:
            s.cell = [vals[0:3], vals[3:6], vals[6:9]]
    for line in lines[2:2 + n]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError("invalid XYZ: atom needs a symbol and three coordinates")
        sym = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        cr, cg, cb = _default_color(sym)
        s.atoms.append(Atom(sym, x, y, z, cr, cg, cb))
    return s


def _save_xyz(s: "Structure", path: str) -> None:
    lines = [str(len(s.atoms)) + "\n"]
    if s.cell:
        lat = " ".join(f"{v:.6f}" for row in s.cell for v in row)
        lines.append(f'Lattice="{lat}" Properties=species:S:1:pos:R:3\n')
    else:
        lines.append("AtomForge structure\n")
    for a in s.atoms:
        lines.append(f"{a.symbol} {a.x:.6f} {a.y:.6f} {a.z:.6f}\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
