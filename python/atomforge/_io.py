"""
File I/O — pure Python, no external dependencies.

Supported formats
-----------------
Read + Write : XYZ, extended-XYZ, VASP POSCAR/CONTCAR, PDB, CIF (P1), LAMMPS data
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Tuple

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


def _clean_symbol(raw: str) -> str:
    """Normalise an element label: strip trailing digits/signs, title-case."""
    s = re.sub(r"[^A-Za-z]", "", raw)
    return s.capitalize() if s else raw


def _cell_from_params(a: float, b: float, c: float,
                      alpha: float, beta: float, gamma: float
                      ) -> List[List[float]]:
    """Build a 3×3 cell matrix from lattice parameters (lengths Å, angles °)."""
    rad = math.pi / 180.0
    ca, cb, cg = math.cos(alpha * rad), math.cos(beta * rad), math.cos(gamma * rad)
    sg = math.sin(gamma * rad)
    ax = a
    bx = b * cg
    by = b * sg
    cx = c * cb
    cy = c * (ca - cb * cg) / sg if sg > 1e-10 else 0.0
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return [[ax, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]]


def _frac_to_cart(fx: float, fy: float, fz: float,
                  cell: List[List[float]]) -> Tuple[float, float, float]:
    a, b, c = cell
    x = fx * a[0] + fy * b[0] + fz * c[0]
    y = fx * a[1] + fy * b[1] + fz * c[1]
    z = fx * a[2] + fy * b[2] + fz * c[2]
    return x, y, z


def _cart_to_frac(x: float, y: float, z: float,
                  cell: List[List[float]]) -> Tuple[float, float, float]:
    """Invert a (generally non-orthogonal) cell matrix to get fractional coords."""
    a, b, c = cell
    # Build 3×3 and invert
    m = [list(col) for col in zip(a, b, c)]
    # Determinant
    det = (m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
         - m[0][1] * (m[1][0]*m[2][2] - m[1][2]*m[2][0])
         + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]))
    if abs(det) < 1e-12:
        return x, y, z
    inv = [
        [(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det,
         (m[0][2]*m[2][1]-m[0][1]*m[2][2])/det,
         (m[0][1]*m[1][2]-m[0][2]*m[1][1])/det],
        [(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det,
         (m[0][0]*m[2][2]-m[0][2]*m[2][0])/det,
         (m[0][2]*m[1][0]-m[0][0]*m[1][2])/det],
        [(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det,
         (m[0][1]*m[2][0]-m[0][0]*m[2][1])/det,
         (m[0][0]*m[1][1]-m[0][1]*m[1][0])/det],
    ]
    fx = inv[0][0]*x + inv[0][1]*y + inv[0][2]*z
    fy = inv[1][0]*x + inv[1][1]*y + inv[1][2]*z
    fz = inv[2][0]*x + inv[2][1]*y + inv[2][2]*z
    return fx, fy, fz


# ── XYZ / extended-XYZ ─────────────────────────────────────────────────────────

def _load_xyz(path: str) -> "Structure":
    from ._structure import Atom, Structure
    s = Structure()
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    if not lines:
        return s
    n = int(lines[0].strip())
    comment = lines[1] if len(lines) > 1 else ""
    m = re.search(r'[Ll]attice="([^"]+)"', comment)
    if m:
        vals = list(map(float, m.group(1).split()))
        if len(vals) == 9:
            s.cell = [vals[0:3], vals[3:6], vals[6:9]]
    for line in lines[2:2 + n]:
        parts = line.split()
        if len(parts) < 4:
            continue
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


# ── VASP POSCAR / CONTCAR ───────────────────────────────────────────────────────

def _load_vasp(path: str) -> "Structure":
    from ._structure import Atom, Structure
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [l.rstrip("\n") for l in fh]

    scale = float(lines[1].split()[0])
    cell = []
    for i in range(2, 5):
        row = [float(v) * scale for v in lines[i].split()]
        cell.append(row)

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

    mode_line = lines[coord_line].strip().lower()
    direct = mode_line.startswith("d")
    coord_line += 1

    s = Structure()
    s.cell = cell
    for sym, count in zip(species, counts):
        sym = _clean_symbol(sym)
        for _ in range(count):
            row = lines[coord_line].split()
            coord_line += 1
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


# ── PDB ────────────────────────────────────────────────────────────────────────

def _load_pdb(path: str) -> "Structure":
    from ._structure import Atom, Structure
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
        a_vec, b_vec, c_vec = s.cell
        a = math.sqrt(sum(v*v for v in a_vec))
        b = math.sqrt(sum(v*v for v in b_vec))
        c = math.sqrt(sum(v*v for v in c_vec))
        cos_alpha = sum(b_vec[i]*c_vec[i] for i in range(3)) / (b * c) if b*c else 0
        cos_beta  = sum(a_vec[i]*c_vec[i] for i in range(3)) / (a * c) if a*c else 0
        cos_gamma = sum(a_vec[i]*b_vec[i] for i in range(3)) / (a * b) if a*b else 0
        alpha = math.degrees(math.acos(max(-1.0, min(1.0, cos_alpha))))
        beta  = math.degrees(math.acos(max(-1.0, min(1.0, cos_beta))))
        gamma = math.degrees(math.acos(max(-1.0, min(1.0, cos_gamma))))
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


# ── CIF (P1 / pre-expanded) ────────────────────────────────────────────────────

def _load_cif(path: str) -> "Structure":
    """
    Load a CIF file.  Reads cell parameters and the first _atom_site loop.
    Handles fractional (fract_x/y/z) and Cartesian (Cartn_x/y/z) coords.
    Symmetry expansion is NOT performed — the file must already be in P1
    (all symmetry-equivalent atoms listed explicitly), which is what
    AtomForge and most MD/DFT export tools produce.
    """
    from ._structure import Atom, Structure

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
    alpha = _float("_cell_angle_alpha") or 90.0
    beta  = _float("_cell_angle_beta")  or 90.0
    gamma = _float("_cell_angle_gamma") or 90.0

    cell = _cell_from_params(a, b, c, alpha, beta, gamma) if (a and b and c) else None

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

    rows = []
    for line in data_text.splitlines():
        line = line.strip()
        if not line or line.startswith("_"):
            break
        # Tokenise, respecting quoted strings
        tokens = re.findall(r"'[^']*'|\"[^\"]*\"|\S+", line)
        if tokens:
            rows.append(tokens)

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
        a_vec, b_vec, c_vec = s.cell
        a = math.sqrt(sum(v*v for v in a_vec))
        b = math.sqrt(sum(v*v for v in b_vec))
        c = math.sqrt(sum(v*v for v in c_vec))
        cos_al = sum(b_vec[i]*c_vec[i] for i in range(3)) / (b*c) if b*c else 0
        cos_be = sum(a_vec[i]*c_vec[i] for i in range(3)) / (a*c) if a*c else 0
        cos_ga = sum(a_vec[i]*b_vec[i] for i in range(3)) / (a*b) if a*b else 0
        alpha = math.degrees(math.acos(max(-1.0, min(1.0, cos_al))))
        beta  = math.degrees(math.acos(max(-1.0, min(1.0, cos_be))))
        gamma = math.degrees(math.acos(max(-1.0, min(1.0, cos_ga))))
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
        "_atom_site_fract_x\n",
        "_atom_site_fract_y\n",
        "_atom_site_fract_z\n",
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


# ── LAMMPS data ────────────────────────────────────────────────────────────────

def _load_lammps(path: str) -> "Structure":
    """
    Load a LAMMPS data file (atom_style atomic or full).
    Element names are read from '# elem1 elem2 ...' comment on the Masses line
    or inferred from atomic mass if the comment is absent.
    """
    from ._structure import Atom, Structure

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

    def _closest_mass(m: float) -> str:
        best, sym = 1e9, "X"
        for mass, s in _MASS_TO_SYM.items():
            if abs(m - mass) < best:
                best, sym = abs(m - mass), s
        return sym if best < 2.0 else "X"

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
                if not mline:
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
                type_to_sym[tid] = sym
            continue

        elif line.startswith("Atoms"):
            i += 1
            while i < len(lines):
                aline = lines[i].strip()
                i += 1
                if not aline:
                    continue  # skip blank separator lines
                if aline[0].isalpha():
                    i -= 1
                    break
                parts = aline.split()
                if len(parts) < 5:
                    continue
                # atomic: id type x y z
                # full:   id mol type charge x y z
                try:
                    if len(parts) >= 7:
                        tid = int(parts[2]); x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                    else:
                        tid = int(parts[1]); x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atoms.append((tid, x, y, z))
                except (ValueError, IndexError):
                    pass
            continue
        i += 1

    s = Structure()
    s.cell = [
        [xhi - xlo,  0.0,       0.0],
        [xy,         yhi - ylo, 0.0],
        [xz,         yz,        zhi - zlo],
    ]
    for tid, x, y, z in atoms:
        sym = type_to_sym.get(tid, f"X{tid}")
        cr, cg, cb = _default_color(sym)
        s.atoms.append(Atom(sym, x, y, z, cr, cg, cb))
    return s


def _save_lammps(s: "Structure", path: str) -> None:
    seen: list[str] = []
    for a in s.atoms:
        if a.symbol not in seen:
            seen.append(a.symbol)
    sym_to_tid = {sym: i + 1 for i, sym in enumerate(seen)}

    cell = s.cell or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    # Only orthorhombic box supported for simplicity
    xlo, xhi = 0.0, math.sqrt(sum(v*v for v in cell[0]))
    ylo, yhi = 0.0, math.sqrt(sum(v*v for v in cell[1]))
    zlo, zhi = 0.0, math.sqrt(sum(v*v for v in cell[2]))

    lines = [
        "AtomForge structure\n\n",
        f"{len(s.atoms)} atoms\n",
        f"{len(seen)} atom types\n\n",
        f"{xlo:.6f} {xhi:.6f} xlo xhi\n",
        f"{ylo:.6f} {yhi:.6f} ylo yhi\n",
        f"{zlo:.6f} {zhi:.6f} zlo zhi\n\n",
        "Masses\n\n",
    ]
    # Approximate masses
    _MASSES = {
        "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
        "Al": 26.982, "Si": 28.086, "Fe": 55.845, "Cu": 63.546,
        "Ni": 58.693, "Ti": 47.867, "W": 183.84, "Au": 196.967,
    }
    for sym in seen:
        mass = _MASSES.get(sym, 1.0)
        lines.append(f"  {sym_to_tid[sym]}  {mass:.3f}  # {sym}\n")
    lines.append("\nAtoms  # atomic\n\n")
    for i, a in enumerate(s.atoms, 1):
        lines.append(f"{i} {sym_to_tid[a.symbol]} {a.x:.6f} {a.y:.6f} {a.z:.6f}\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)


# ── Format dispatch ─────────────────────────────────────────────────────────────

_FORMAT_MAP: dict[str, tuple] = {
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


def load(path: str) -> "Structure":
    """Load a structure from file.  Format is inferred from the file extension."""
    ext = Path(path).suffix.lower()
    loader, _ = _FORMAT_MAP.get(ext, (_load_xyz, _save_xyz))
    # If the filename has no extension or is POSCAR/CONTCAR, try VASP
    stem = Path(path).stem.upper()
    if not ext and stem in ("POSCAR", "CONTCAR"):
        loader = _load_vasp
    return loader(path)


def save(s: "Structure", path: str) -> None:
    """Save a structure to file.  Format is inferred from the file extension."""
    ext = Path(path).suffix.lower()
    _, saver = _FORMAT_MAP.get(ext, (_load_xyz, _save_xyz))
    saver(s, path)
