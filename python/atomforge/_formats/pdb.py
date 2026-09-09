"""Native PDB structure reader and writer."""
from __future__ import annotations

from typing import TYPE_CHECKING

from .._elements import _default_color, _clean_symbol
from .._cell import _params_from_cell, _cell_from_params

if TYPE_CHECKING:
    from .._structure import Structure


def _load_pdb(path: str) -> "Structure":
    from .._structure import Atom, Structure
    s = Structure()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "CRYST1":
                try:
                    a  = float(line[6:15])
                    b  = float(line[15:24])
                    c  = float(line[24:33])
                    al = float(line[33:40])
                    be = float(line[40:47])
                    ga = float(line[47:54])
                    s.cell = _cell_from_params(a, b, c, al, be, ga)
                except (ValueError, IndexError):
                    pass
            elif rec in ("ATOM", "HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    # Element column (cols 76-78) preferred, fall back to name col
                    elem = line[76:78].strip() if len(line) > 76 else ""
                    if not elem:
                        elem = line[12:16].strip()
                    sym = _clean_symbol(elem)
                    cr, cg, cb = _default_color(sym)
                    s.atoms.append(Atom(sym, x, y, z, cr, cg, cb))
                except (ValueError, IndexError):
                    pass
    return s


def _save_pdb(s: "Structure", path: str) -> None:
    lines = []
    if s.cell:
        a, b, c, alpha, beta, gamma = _params_from_cell(s.cell)
        lines.append(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1\n")
    for i, atom in enumerate(s.atoms, 1):
        name = f"{atom.symbol:<4}"
        lines.append(
            f"HETATM{i:5d} {name} UNK A{1:4d}    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
            f"  1.00  0.00          {atom.symbol:>2}\n"
        )
    lines.append("END\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
