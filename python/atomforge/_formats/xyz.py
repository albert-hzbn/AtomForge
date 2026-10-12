"""Native XYZ structure reader and writer."""
from __future__ import annotations

import re
from typing import TYPE_CHECKING

from .._elements import _default_color

if TYPE_CHECKING:
    from .._structure import Structure


def _element_of(label: str) -> str:
    """An XYZ label's element, as the desktop reads it: a symbol in any case
    (CU, cu), a label starting with one (Cu1, OW) or an atomic number (29).
    Other labels are kept as written."""
    from .._elements import _MASSES
    if label.isdigit():
        symbols = list(_MASSES)
        number = int(label)
        return symbols[number - 1] if 1 <= number <= len(symbols) else label
    letters = re.match(r"[A-Za-z]{1,2}", label)
    if letters:
        symbol = letters.group(0).capitalize()
        if symbol in _MASSES:
            return symbol
        if symbol[0] in _MASSES:  # OW, HW1: a one-letter element and a tag
            return symbol[0]
    return label


def _columns(comment: str) -> tuple[int, int]:
    """Species and position columns from an extended XYZ Properties key
    (e.g. id:I:1:species:S:1:pos:R:3; an atomic-number column Z:I:1 stands
    in for a missing species column); plain XYZ is (0, 1)."""
    m = re.search(r'[Pp]roperties=(\S+)', comment)
    if not m:
        return 0, 1
    fields = m.group(1).split(":")
    if len(fields) % 3:
        raise ValueError("invalid extended XYZ Properties")
    species = pos = number = None
    column = 0
    for name, _kind, count in zip(fields[0::3], fields[1::3], fields[2::3]):
        width = int(count)
        if width < 1:
            raise ValueError("invalid extended XYZ Properties")
        if name.lower() == "species" and width == 1:
            species = column
        if name == "Z" and width == 1:
            number = column
        if name.lower() in ("pos", "positions") and width == 3:
            pos = column
        column += width
    if species is None:
        species = number
    if species is None or pos is None:
        raise ValueError("extended XYZ needs species (or Z) and pos columns")
    return species, pos


def _load_xyz(path: str) -> "Structure":
    from .._structure import Atom, Structure
    s = Structure()
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    while lines and not lines[0].strip():  # leading blank lines
        lines.pop(0)
    if not lines:
        return s
    n = int(lines[0].strip())
    if n < 0 or len(lines) < n + 2:
        raise ValueError("invalid XYZ: incomplete file or negative atom count")
    comment = lines[1] if len(lines) > 1 else ""
    m = re.search(r"""[Ll]attice=(?:"([^"]+)"|'([^']+)')""", comment)
    if m:
        vals = list(map(float, (m.group(1) or m.group(2)).split()))
        if len(vals) == 9:
            s.cell = [vals[0:3], vals[3:6], vals[6:9]]
    species_col, pos_col = _columns(comment)
    for line in lines[2:2 + n]:
        parts = line.split()
        if len(parts) <= max(species_col, pos_col + 2):
            raise ValueError("invalid XYZ: atom needs a symbol and three coordinates")
        # Element symbols in any case ("CU", "cu" -> "Cu"); other labels as written.
        sym = _element_of(parts[species_col])
        x, y, z = (float(v) for v in parts[pos_col:pos_col + 3])
        cr, cg, cb = _default_color(sym)
        s.atoms.append(Atom(sym, x, y, z, cr, cg, cb))
    return s


def _save_xyz(s: "Structure", path: str) -> None:
    lines = [str(len(s.atoms)) + "\n"]
    if s.cell:
        lat = " ".join(f"{v:.10g}" for row in s.cell for v in row)
        lines.append(f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
    else:
        lines.append("AtomForge structure\n")
    for a in s.atoms:
        lines.append(f"{a.symbol} {a.x:.10f} {a.y:.10f} {a.z:.10f}\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
