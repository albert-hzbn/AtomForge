"""Comprehensive tests for the atomforge package — no external dependencies."""

import math
import os
import sys
import tempfile
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import atomforge as af
from atomforge._io import (
    _load_xyz, _save_xyz,
    _load_vasp, _save_vasp,
    _load_pdb, _save_pdb,
    _load_cif, _save_cif,
    _load_lammps, _save_lammps,
    _cell_from_params, _frac_to_cart, _cart_to_frac,
)

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"
_failures = []


def check(name: str, cond: bool, detail: str = ""):
    if cond:
        print(f"  {PASS}  {name}")
    else:
        msg = f"  {FAIL}  {name}" + (f": {detail}" if detail else "")
        print(msg)
        _failures.append(name)


def approx(a: float, b: float, tol: float = 1e-4) -> bool:
    return abs(a - b) < tol


def _bcc_fe() -> af.Structure:
    s = af.Structure()
    s.set_cell(2.87, 2.87, 2.87)
    s.add_atom("Fe", 0.0,   0.0,   0.0)
    s.add_atom("Fe", 1.435, 1.435, 1.435)
    return s


def _water() -> af.Structure:
    s = af.Structure()
    s.add_atom("O", 0.0,  0.0,  0.0)
    s.add_atom("H", 0.96, 0.0,  0.0)
    s.add_atom("H", 0.0,  0.96, 0.0)
    return s


# ── Structure / Atom API ────────────────────────────────────────────────────────

def test_structure_api():
    print("\n[Structure API]")
    s = _bcc_fe()
    check("len", len(s) == 2)
    check("repr periodic", "periodic" in repr(s))
    check("atom symbol", s.atoms[0].symbol == "Fe")
    check("atom coords", approx(s.atoms[1].x, 1.435))

    # copy
    c = s.copy()
    c.atoms[0].x = 99.0
    check("copy is independent", approx(s.atoms[0].x, 0.0))

    # translate
    s2 = _water()
    s2.translate(1, 2, 3)
    check("translate x", approx(s2.atoms[0].x, 1.0))
    check("translate y", approx(s2.atoms[0].y, 2.0))

    # scale
    s3 = _bcc_fe()
    s3.scale(2.0)
    check("scale atom", approx(s3.atoms[1].x, 2.87))
    check("scale cell", approx(s3.cell[0][0], 5.74))

    # remove_atom
    s4 = _water()
    removed = s4.remove_atom(0)
    check("remove_atom symbol", removed.symbol == "O")
    check("remove_atom len", len(s4) == 2)

    # filter_species
    s5 = _bcc_fe().repeat(2, 2, 2)
    fe = s5.filter_species("Fe")
    check("filter_species", len(fe) == 16)

    # repeat
    sup = _bcc_fe().repeat(3, 3, 3)
    check("repeat atom count", len(sup) == 54)
    check("repeat cell", approx(sup.cell[0][0], 2.87 * 3))

    # non-periodic repr
    check("repr non-periodic", "non-periodic" in repr(_water()))


# ── Math helpers ────────────────────────────────────────────────────────────────

def test_math():
    print("\n[Math helpers]")
    cell = _cell_from_params(4.0, 4.0, 4.0, 90, 90, 90)
    check("cubic cell a", approx(cell[0][0], 4.0))
    check("cubic cell off-diag", approx(cell[0][1], 0.0))

    # frac→cart→frac round-trip
    cell_hex = _cell_from_params(3.0, 3.0, 5.0, 90, 90, 120)
    fx, fy, fz = 0.333, 0.667, 0.5
    x, y, z = _frac_to_cart(fx, fy, fz, cell_hex)
    rfx, rfy, rfz = _cart_to_frac(x, y, z, cell_hex)
    check("frac->cart->frac x", approx(rfx, fx, 1e-6))
    check("frac->cart->frac y", approx(rfy, fy, 1e-6))
    check("frac->cart->frac z", approx(rfz, fz, 1e-6))


# ── XYZ / extXYZ ───────────────────────────────────────────────────────────────

def test_xyz():
    print("\n[XYZ / extXYZ]")
    with tempfile.NamedTemporaryFile(suffix=".xyz", delete=False, mode="w") as f:
        name = f.name

    try:
        # Cluster (no cell)
        s = _water()
        _save_xyz(s, name)
        r = _load_xyz(name)
        check("xyz cluster atom count", len(r) == 3)
        check("xyz cluster symbol", r.atoms[1].symbol == "H")
        check("xyz cluster x", approx(r.atoms[1].x, 0.96))
        check("xyz cluster no cell", r.cell is None)

        # Periodic (extXYZ)
        s2 = _bcc_fe()
        _save_xyz(s2, name)
        r2 = _load_xyz(name)
        check("extxyz cell present", r2.cell is not None)
        check("extxyz cell a", approx(r2.cell[0][0], 2.87))
        check("extxyz atom count", len(r2) == 2)

        # Round-trip via af.load / af.save
        af.save(s2, name)
        r3 = af.load(name)
        check("af.load/save xyz", len(r3) == 2)
    finally:
        os.unlink(name)


# ── VASP POSCAR ─────────────────────────────────────────────────────────────────

_POSCAR_BCC = """\
BCC Fe
1.0
  2.870000  0.000000  0.000000
  0.000000  2.870000  0.000000
  0.000000  0.000000  2.870000
Fe
2
Direct
  0.000000  0.000000  0.000000
  0.500000  0.500000  0.500000
"""

def test_vasp():
    print("\n[VASP POSCAR]")
    with tempfile.NamedTemporaryFile(suffix=".vasp", delete=False, mode="w") as f:
        f.write(_POSCAR_BCC)
        name = f.name
    try:
        s = _load_vasp(name)
        check("vasp atom count", len(s) == 2)
        check("vasp symbol", s.atoms[0].symbol == "Fe")
        check("vasp cell a", approx(s.cell[0][0], 2.87))
        check("vasp frac->cart", approx(s.atoms[1].x, 1.435))

        # Round-trip
        with tempfile.NamedTemporaryFile(suffix=".vasp", delete=False) as f2:
            name2 = f2.name
        _save_vasp(s, name2)
        s2 = _load_vasp(name2)
        check("vasp round-trip count", len(s2) == 2)
        check("vasp round-trip x1", approx(s2.atoms[1].x, 1.435, 1e-3))
        os.unlink(name2)

        # af.load dispatch
        s3 = af.load(name)
        check("af.load .vasp", s3.atoms[0].symbol == "Fe")
    finally:
        os.unlink(name)


# ── PDB ────────────────────────────────────────────────────────────────────────

_PDB_WATER = """\
CRYST1    5.000    5.000    5.000  90.00  90.00  90.00 P 1           1
HETATM    1 O    UNK A   1       0.000   0.000   0.000  1.00  0.00           O
HETATM    2 H    UNK A   1       0.960   0.000   0.000  1.00  0.00           H
HETATM    3 H    UNK A   1       0.000   0.960   0.000  1.00  0.00           H
END
"""

def test_pdb():
    print("\n[PDB]")
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w") as f:
        f.write(_PDB_WATER)
        name = f.name
    try:
        s = _load_pdb(name)
        check("pdb atom count", len(s) == 3)
        check("pdb O symbol", s.atoms[0].symbol == "O")
        check("pdb H x", approx(s.atoms[1].x, 0.96))
        check("pdb cell a", approx(s.cell[0][0], 5.0))

        # Round-trip
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f2:
            name2 = f2.name
        _save_pdb(s, name2)
        s2 = _load_pdb(name2)
        check("pdb round-trip count", len(s2) == 3)
        check("pdb round-trip cell", approx(s2.cell[0][0], 5.0, 0.01))
        check("pdb round-trip H x", approx(s2.atoms[1].x, 0.96, 0.01))
        os.unlink(name2)

        # af.load dispatch
        s3 = af.load(name)
        check("af.load .pdb", s3.atoms[0].symbol == "O")
    finally:
        os.unlink(name)


# ── CIF ────────────────────────────────────────────────────────────────────────

_CIF_BCC = """\
data_bcc_fe

_cell_length_a   2.8700
_cell_length_b   2.8700
_cell_length_c   2.8700
_cell_angle_alpha   90.0000
_cell_angle_beta    90.0000
_cell_angle_gamma   90.0000
_symmetry_space_group_name_H-M  'P 1'

loop_
_atom_site_type_symbol
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe  Fe1  0.000000  0.000000  0.000000
Fe  Fe2  0.500000  0.500000  0.500000
"""

def test_cif():
    print("\n[CIF]")
    with tempfile.NamedTemporaryFile(suffix=".cif", delete=False, mode="w") as f:
        f.write(_CIF_BCC)
        name = f.name
    try:
        s = _load_cif(name)
        check("cif atom count", len(s) == 2)
        check("cif symbol", s.atoms[0].symbol == "Fe")
        check("cif cell a", approx(s.cell[0][0], 2.87))
        # fractional 0.5,0.5,0.5 in cubic cell → Cartesian 1.435,1.435,1.435
        check("cif frac->cart x", approx(s.atoms[1].x, 1.435))

        # Round-trip
        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as f2:
            name2 = f2.name
        _save_cif(s, name2)
        s2 = _load_cif(name2)
        check("cif round-trip count", len(s2) == 2)
        check("cif round-trip x1", approx(s2.atoms[1].x, 1.435, 1e-3))
        os.unlink(name2)

        # af.load dispatch
        s3 = af.load(name)
        check("af.load .cif", s3.atoms[0].symbol == "Fe")
    finally:
        os.unlink(name)


# ── LAMMPS ─────────────────────────────────────────────────────────────────────

_LAMMPS_BCC = """\
BCC Fe

2 atoms
1 atom types

0.000000 2.870000 xlo xhi
0.000000 2.870000 ylo yhi
0.000000 2.870000 zlo zhi

Masses

  1  55.845  # Fe

Atoms  # atomic

1 1 0.000000 0.000000 0.000000
2 1 1.435000 1.435000 1.435000
"""

def test_lammps():
    print("\n[LAMMPS data]")
    with tempfile.NamedTemporaryFile(suffix=".lmp", delete=False, mode="w") as f:
        f.write(_LAMMPS_BCC)
        name = f.name
    try:
        s = _load_lammps(name)
        check("lammps atom count", len(s) == 2)
        check("lammps symbol from mass", s.atoms[0].symbol == "Fe")
        check("lammps x1", approx(s.atoms[1].x, 1.435))
        check("lammps cell a", approx(s.cell[0][0], 2.87))

        # Round-trip
        with tempfile.NamedTemporaryFile(suffix=".lmp", delete=False) as f2:
            name2 = f2.name
        _save_lammps(s, name2)
        s2 = _load_lammps(name2)
        check("lammps round-trip count", len(s2) == 2)
        check("lammps round-trip x1", approx(s2.atoms[1].x, 1.435, 1e-3))
        os.unlink(name2)

        # af.load dispatch
        s3 = af.load(name)
        check("af.load .lmp", s3.atoms[0].symbol == "Fe")
    finally:
        os.unlink(name)


# ── Cross-format round-trips ────────────────────────────────────────────────────

def test_cross_format():
    print("\n[Cross-format round-trips]")
    s = _bcc_fe().repeat(2, 2, 2)   # 16-atom supercell
    fmts = [".xyz", ".vasp", ".pdb", ".cif", ".lmp"]
    for ext in fmts:
        with tempfile.NamedTemporaryFile(suffix=ext, delete=False) as f:
            name = f.name
        try:
            af.save(s, name)
            r = af.load(name)
            check(f"cross {ext} count", len(r) == 16,
                  f"got {len(r)}")
            check(f"cross {ext} symbol", r.atoms[0].symbol == "Fe")
        finally:
            os.unlink(name)


# ── No external imports ─────────────────────────────────────────────────────────

def test_no_external_deps():
    print("\n[No external dependencies]")
    import importlib, pkgutil
    # Verify ase is not imported anywhere in the package
    import atomforge
    import atomforge._io
    import atomforge._structure
    import atomforge._viewer
    for mod_name in ["atomforge", "atomforge._io", "atomforge._structure", "atomforge._viewer"]:
        mod = sys.modules.get(mod_name)
        if mod is None:
            continue
        src = getattr(mod, "__file__", "") or ""
        if src.endswith(".py"):
            with open(src, encoding="utf-8") as fh:
                text = fh.read()
            has_ase = "import ase" in text or "from ase" in text
            check(f"no ase import in {mod_name}", not has_ase)


# ── Run all ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("atomforge test suite")
    print("=" * 60)
    for fn in [
        test_structure_api,
        test_math,
        test_xyz,
        test_vasp,
        test_pdb,
        test_cif,
        test_lammps,
        test_cross_format,
        test_no_external_deps,
    ]:
        try:
            fn()
        except Exception:
            print(f"\n  {FAIL}  {fn.__name__} raised an exception:")
            traceback.print_exc()
            _failures.append(fn.__name__)

    print("\n" + "=" * 60)
    if _failures:
        print(f"FAILED  ({len(_failures)} failures): {', '.join(_failures)}")
        sys.exit(1)
    else:
        print(f"All tests passed.")
    print("=" * 60)
