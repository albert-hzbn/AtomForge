"""Native CIF structure reader and writer."""
from __future__ import annotations

import re
from typing import TYPE_CHECKING, Optional

from .._elements import _default_color, _clean_symbol
from .._cell import _params_from_cell, _cell_from_params, _frac_to_cart, _cart_to_frac

if TYPE_CHECKING:
    from .._structure import Structure


def _load_cif(path: str) -> "Structure":
    """
    Load a CIF file.  Reads cell parameters and the first _atom_site loop.
    Handles fractional (fract_x/y/z) and Cartesian (Cartn_x/y/z) coords.
    Symmetry expansion is NOT performed — the file must already be in P1
    (all symmetry-equivalent atoms listed explicitly), which is what
    AtomForge and most MD/DFT export tools produce.
    """
    from .._structure import Atom, Structure

    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()

    # Remove CIF comments
    text = re.sub(r"#[^\n]*", "", text)

    def _get_scalar(tag: str) -> Optional[str]:
        m = re.search(rf"(?i){re.escape(tag)}\s+(\S+)", text)
        return m.group(1) if m else None

    def _float(tag: str) -> Optional[float]:
        v = _get_scalar(tag)
        if v is None:
            return None
        # Strip uncertainty like "3.456(7)"
        v = re.sub(r"\(\d+\)$", "", v)
        try:
            return float(v)
        except ValueError:
            return None

    a     = _float("_cell_length_a")
    b     = _float("_cell_length_b")
    c     = _float("_cell_length_c")
    alpha = _float("_cell_angle_alpha")
    beta  = _float("_cell_angle_beta")
    gamma = _float("_cell_angle_gamma")
    alpha = 90.0 if alpha is None else alpha
    beta = 90.0 if beta is None else beta
    gamma = 90.0 if gamma is None else gamma

    cell = _cell_from_params(a, b, c, alpha, beta, gamma) if all(v is not None for v in (a, b, c)) else None

    # Find the _atom_site loop
    loop_m = re.search(r"(?i)loop_\s*((?:_atom_site_\S+\s*)+)", text)
    if not loop_m:
        s = Structure()
        s.cell = cell
        return s

    header_block = loop_m.group(1)
    tags = [t.lower() for t in re.findall(r"_atom_site_\S+", header_block)]

    # Everything after the header until the next loop_ or data_ block
    data_start = loop_m.end()
    data_end   = re.search(r"(?i)(loop_|data_)", text[data_start:])
    data_text  = text[data_start: data_start + data_end.start()] if data_end else text[data_start:]

    tokens = []
    for line in data_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("_"):
            break
        # Tokenise, respecting quoted strings
        tokens.extend(re.findall(r"'[^']*'|\"[^\"]*\"|\S+", line))
    if len(tokens) % len(tags):
        raise ValueError("invalid CIF: incomplete atom-site loop row")
    rows = [tokens[i:i + len(tags)] for i in range(0, len(tokens), len(tags))]

    def col(tag_suffix: str) -> int:
        for i, t in enumerate(tags):
            if t.endswith(tag_suffix):
                return i
        return -1

    idx_sym    = col("type_symbol")
    idx_label  = col("label")
    idx_fx     = col("fract_x")
    idx_fy     = col("fract_y")
    idx_fz     = col("fract_z")
    idx_cx     = col("cartn_x")
    idx_cy     = col("cartn_y")
    idx_cz     = col("cartn_z")

    s = Structure()
    s.cell = cell

    def val(row: list, idx: int) -> str:
        if idx < 0 or idx >= len(row):
            return ""
        return re.sub(r"\(\d+\)$", "", row[idx].strip("'\""))

    for row in rows:
        raw_sym = val(row, idx_sym) or val(row, idx_label)
        if not raw_sym:
            continue
        sym = _clean_symbol(raw_sym)

        if idx_fx >= 0 and cell:
            try:
                fx, fy, fz = float(val(row, idx_fx)), float(val(row, idx_fy)), float(val(row, idx_fz))
                x, y, z = _frac_to_cart(fx, fy, fz, cell)
            except ValueError:
                continue
        elif idx_cx >= 0:
            try:
                x, y, z = float(val(row, idx_cx)), float(val(row, idx_cy)), float(val(row, idx_cz))
            except ValueError:
                continue
        else:
            continue

        cr, cg, cb = _default_color(sym)
        s.atoms.append(Atom(sym, x, y, z, cr, cg, cb))
    return s


def _save_cif(s: "Structure", path: str) -> None:
    lines = ["data_atomforge\n\n"]
    if s.cell:
        a, b, c, alpha, beta, gamma = _params_from_cell(s.cell)
        lines += [
            f"_cell_length_a   {a:.6f}\n",
            f"_cell_length_b   {b:.6f}\n",
            f"_cell_length_c   {c:.6f}\n",
            f"_cell_angle_alpha   {alpha:.4f}\n",
            f"_cell_angle_beta    {beta:.4f}\n",
            f"_cell_angle_gamma   {gamma:.4f}\n",
            "_symmetry_space_group_name_H-M  'P 1'\n\n",
        ]
    lines += [
        "loop_\n",
        "_atom_site_type_symbol\n",
        "_atom_site_label\n",
        "_atom_site_fract_x\n" if s.cell else "_atom_site_Cartn_x\n",
        "_atom_site_fract_y\n" if s.cell else "_atom_site_Cartn_y\n",
        "_atom_site_fract_z\n" if s.cell else "_atom_site_Cartn_z\n",
    ]
    counts: dict[str, int] = {}
    for a in s.atoms:
        counts[a.symbol] = counts.get(a.symbol, 0) + 1
        label = f"{a.symbol}{counts[a.symbol]}"
        if s.cell:
            fx, fy, fz = _cart_to_frac(a.x, a.y, a.z, s.cell)
        else:
            fx, fy, fz = a.x, a.y, a.z
        lines.append(f"{a.symbol}  {label}  {fx:.6f}  {fy:.6f}  {fz:.6f}\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
