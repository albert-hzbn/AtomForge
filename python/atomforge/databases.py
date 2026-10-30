"""Structures from online databases: the Crystallography Open Database (COD)
and the Materials Project.

    import atomforge as af
    from atomforge.databases import fetch_cod, fetch_materials_project, search_cod, search_materials_project
    rutile = fetch_cod(9015662)                      # open access
    hits = search_cod(elements=["Ti", "O"], strict=True)    # entries with exactly Ti and O
    iron = fetch_materials_project("mp-13")           # needs an API key (MP_API_KEY)
    stable = search_materials_project(formula="TiO2", stable_only=True)

Both download over HTTPS when called; nothing is fetched on import.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Callable, Optional

from ._structure import Atom, Structure
from ._elements import _default_color

COD_URL = "https://www.crystallography.net/cod/{id}.cif"
COD_SEARCH_URL = "https://www.crystallography.net/cod/result"
MP_URL = "https://api.materialsproject.org/materials/summary/"


class DatabaseError(RuntimeError):
    """A database request failed or returned no usable structure."""


def _get(url: str, headers: Optional[dict] = None, timeout: float = 30.0, opener: Optional[Callable] = None) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "AtomForge", **(headers or {})})
    try:
        with (opener or urllib.request.urlopen)(request, timeout=timeout) as response:
            return response.read()
    except Exception as error:  # network errors, HTTP errors
        raise DatabaseError(f"Could not download {url}: {error}") from error


def _structure_from_cif_text(text: str) -> Structure:
    """A CIF (with symmetry) as explicit atoms: the Python reader, or the
    AtomForge executable's reader when the crystallography extra is missing."""
    from ._io import load
    with tempfile.TemporaryDirectory(prefix="atomforge_db_") as folder:
        path = Path(folder) / "structure.cif"
        path.write_text(text, encoding="utf-8")
        try:
            return load(str(path))
        except ImportError:
            from .builders import _executable
            converted = Path(folder) / "structure.xyz"
            done = subprocess.run([_executable(None), "--convert", "--input", str(path), "--output", str(converted)],
                                  capture_output=True, text=True, timeout=120)
            if done.returncode != 0 or not converted.exists():
                raise DatabaseError("Reading a CIF with symmetry needs atomforge-py[crystallography] or the AtomForge executable "
                                    f"(ATOMFORGE_PATH): {done.stderr.strip()}")
            return load(str(converted))


def fetch_cod(cod_id: int, *, timeout: float = 30.0, opener: Optional[Callable] = None) -> Structure:
    """The structure of a COD entry (e.g. 9008460) from its CIF, symmetry expanded."""
    if not str(cod_id).isdigit():
        raise ValueError(f"COD ids are numbers, not {cod_id!r}")
    text = _get(COD_URL.format(id=int(cod_id)), timeout=timeout, opener=opener).decode("utf-8", errors="replace")
    if "_atom_site" not in text:
        raise DatabaseError(f"COD entry {cod_id} has no atom sites")
    return _structure_from_cif_text(text)


def _hill_formula(formula: str) -> str:
    """COD's formula search wants Hill order with spaces ("O2 Ti"); accept "TiO2" too."""
    import re
    parts = re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula.replace(" ", ""))
    if not parts or "".join(e + n for e, n in parts) != formula.replace(" ", ""):
        raise ValueError(f"Not a chemical formula: {formula!r}")
    counts: dict = {}
    for element, number in parts:
        counts[element] = counts.get(element, 0) + (float(number) if number else 1)
    order = (["C", "H"] if "C" in counts else []) + sorted(e for e in counts if not ("C" in counts and e in ("C", "H")))
    text = lambda n: str(int(n)) if n == int(n) else str(n)
    return " ".join(e + (text(counts[e]) if counts[e] != 1 else "") for e in order)


def search_cod(formula: Optional[str] = None, elements: Optional[list] = None, *, strict: bool = False, text: Optional[str] = None,
               max_results: int = 100, timeout: float = 60.0, opener: Optional[Callable] = None) -> list:
    """COD entries matching a formula ("TiO2" or Hill "O2 Ti"), containing the
    elements (exactly those with strict=True) and/or a word in the title or
    mineral name. Each hit is a dict: id, formula, spacegroup, a, b, c, alpha,
    beta, gamma, volume, mineral, title, year; pass the id to fetch_cod."""
    query = [("format", "json")]
    if formula:
        query.append(("formula", _hill_formula(formula)))
    for k, element in enumerate(elements or [], start=1):
        if k > 8:
            raise ValueError("COD searches take at most 8 elements")
        query.append((f"el{k}", element))
    if elements and strict:
        query += [("strictmin", str(len(elements))), ("strictmax", str(len(elements)))]
    if text:
        query.append(("text", text))
    if len(query) == 1:
        raise ValueError("Give a formula, elements or text to search for")
    raw = _get(f"{COD_SEARCH_URL}?{urllib.parse.urlencode(query)}", timeout=timeout, opener=opener).decode("utf-8", errors="replace")
    try:
        rows = json.loads(raw) if raw.strip() else []
    except ValueError as error:
        raise DatabaseError(f"Unexpected COD search reply: {raw[:200]}") from error
    number = lambda v: float(v) if v not in (None, "") else None
    hits = []
    for row in rows[:max_results]:
        hits.append({
            "id": int(row.get("file")),
            "formula": (row.get("formula") or "").strip("- "),
            "spacegroup": row.get("sg"),
            "a": number(row.get("a")), "b": number(row.get("b")), "c": number(row.get("c")),
            "alpha": number(row.get("alpha")), "beta": number(row.get("beta")), "gamma": number(row.get("gamma")),
            "volume": number(row.get("vol")),
            "mineral": row.get("mineral"), "title": row.get("title"), "year": row.get("year"),
        })
    return hits


def search_materials_project(formula: Optional[str] = None, elements: Optional[list] = None, *, chemsys: Optional[str] = None,
                             stable_only: bool = False, max_energy_above_hull: Optional[float] = None, max_results: int = 100,
                             api_key: Optional[str] = None, timeout: float = 60.0, opener: Optional[Callable] = None) -> list:
    """Materials Project entries by formula ("TiO2"), by elements they must
    contain, or by chemical system ("Li-Fe-O": only those elements); optionally
    only stable ones or within an energy above the convex hull (eV/atom).
    Each hit is a dict: id, formula, spacegroup, crystal_system, nsites,
    energy_above_hull, band_gap, is_stable, is_magnetic; sorted by stability."""
    key = api_key or os.environ.get("MP_API_KEY")
    if not key:
        raise DatabaseError("The Materials Project needs an API key: pass api_key= or set MP_API_KEY")
    fields = "material_id,formula_pretty,symmetry,nsites,energy_above_hull,band_gap,is_stable,is_magnetic"
    query = {"_fields": fields, "_limit": str(max_results)}
    if formula:
        query["formula"] = formula
    if elements:
        query["elements"] = ",".join(elements)
    if chemsys:
        query["chemsys"] = chemsys
    if len(query) == 2:
        raise ValueError("Give a formula, elements or a chemical system to search for")
    if stable_only:
        query["is_stable"] = "true"
    if max_energy_above_hull is not None:
        query["energy_above_hull_max"] = str(max_energy_above_hull)
    payload = json.loads(_get(f"{MP_URL}?{urllib.parse.urlencode(query)}", headers={"X-API-KEY": key}, timeout=timeout,
                              opener=opener).decode("utf-8"))
    hits = []
    for entry in payload.get("data") or []:
        symmetry = entry.get("symmetry") or {}
        hits.append({
            "id": entry.get("material_id"), "formula": entry.get("formula_pretty"),
            "spacegroup": symmetry.get("symbol"), "crystal_system": symmetry.get("crystal_system"),
            "nsites": entry.get("nsites"), "energy_above_hull": entry.get("energy_above_hull"),
            "band_gap": entry.get("band_gap"), "is_stable": entry.get("is_stable"), "is_magnetic": entry.get("is_magnetic"),
        })
    hits.sort(key=lambda h: (h["energy_above_hull"] is None, h["energy_above_hull"] or 0.0))
    return hits


def conventional_cell(structure: Structure, symprec: float = 1e-3) -> Structure:
    """The conventional standard cell (spglib), keeping magnetic moments."""
    try:
        import spglib
    except ImportError as error:
        raise DatabaseError("The conventional cell needs spglib (pip install spglib)") from error
    # spglib needs distinct integer types only.
    kinds = sorted({a.symbol for a in structure.atoms})
    numbers = [kinds.index(a.symbol) + 1 for a in structure.atoms]
    lattice = [list(row) for row in structure.cell]
    import numpy as np
    inverse = np.linalg.inv(np.array(lattice))
    positions = [list(np.array([a.x, a.y, a.z]) @ inverse) for a in structure.atoms]
    dataset = spglib.get_symmetry_dataset((lattice, positions, numbers), symprec=symprec)
    if dataset is None:
        raise DatabaseError("spglib found no symmetry for this structure")
    get = (lambda name: getattr(dataset, name)) if hasattr(dataset, "std_lattice") else (lambda name: dataset[name])
    moment_of_primitive = {}
    for i, p in enumerate(get("mapping_to_primitive")):
        moment_of_primitive.setdefault(int(p), structure.atoms[i].magmom)
    out = Structure()
    std = np.array(get("std_lattice"))
    out.cell = std.tolist()
    for j, (f, z) in enumerate(zip(get("std_positions"), get("std_types"))):
        x, y, zc = (np.array(f) % 1.0) @ std
        symbol = kinds[int(z) - 1]
        atom = Atom(symbol, float(x), float(y), float(zc), *_default_color(symbol))
        atom.magmom = moment_of_primitive.get(int(get("std_mapping_to_primitive")[j]), 0.0)
        out.atoms.append(atom)
    return out


def _structure_from_pymatgen_dict(data: dict) -> Structure:
    """A structure from pymatgen's Structure.as_dict() layout (lattice matrix,
    sites with species and Cartesian xyz, optional magmom properties)."""
    structure = Structure()
    try:
        structure.cell = [[float(v) for v in row] for row in data["lattice"]["matrix"]]
        for site in data["sites"]:
            species = max(site["species"], key=lambda s: s.get("occu", 1.0))
            symbol = species["element"]
            x, y, z = (float(v) for v in site["xyz"])
            atom = Atom(symbol, x, y, z, *_default_color(symbol))
            magmom = (site.get("properties") or {}).get("magmom")
            if isinstance(magmom, (int, float)):
                atom.magmom = float(magmom)
            structure.atoms.append(atom)
    except (KeyError, TypeError, ValueError) as error:
        raise DatabaseError(f"Unexpected structure data: {error}") from error
    if not structure.atoms:
        raise DatabaseError("The structure has no sites")
    return structure


def main(argv=None) -> int:
    """python -m atomforge.databases cod 9008460 out.cif | mp mp-13 out.xyz [--conventional]
    | search cod TiO2 | search mp Fe2O3 --stable"""
    import argparse
    from ._io import save
    parser = argparse.ArgumentParser(prog="python -m atomforge.databases", description="Download or search structures in COD or the Materials Project.")
    parser.add_argument("database", choices=["cod", "mp", "search"])
    parser.add_argument("id", help="entry id; for search: cod or mp")
    parser.add_argument("output", help="structure file to write (format from the extension); for search: the formula")
    parser.add_argument("--api-key", help="Materials Project API key (default: MP_API_KEY)")
    parser.add_argument("--conventional", action="store_true", help="Materials Project: the conventional standard cell")
    parser.add_argument("--stable", action="store_true", help="search mp: stable entries only")
    args = parser.parse_args(argv)
    try:
        if args.database == "search":
            if args.id == "cod":
                for hit in search_cod(formula=args.output):
                    print(f"{hit['id']}  {hit['formula']:<20} {hit['spacegroup'] or '':<12} {hit['mineral'] or ''}")
            elif args.id == "mp":
                for hit in search_materials_project(formula=args.output, stable_only=args.stable, api_key=args.api_key):
                    print(f"{hit['id']:<12} {hit['formula']:<12} {hit['spacegroup'] or '':<10} E_hull {hit['energy_above_hull']} eV/atom")
            else:
                raise ValueError("search cod or search mp")
            return 0
        structure = (fetch_cod(args.id) if args.database == "cod"
                     else fetch_materials_project(args.id, api_key=args.api_key, conventional=args.conventional))
    except (DatabaseError, ValueError) as error:
        print(f"Error: {error}")
        return 1
    save(structure, args.output)
    print(f"Saved {len(structure.atoms)} atoms to {args.output}")
    return 0


def fetch_materials_project(material_id: str, *, api_key: Optional[str] = None, conventional: bool = False, timeout: float = 30.0,
                            opener: Optional[Callable] = None) -> Structure:
    """A Materials Project structure (e.g. "mp-149") with its magnetic moments
    where the database gives them, as stored (usually primitive) or as the
    conventional standard cell (needs spglib). The API key comes from
    ``api_key`` or the environment variable MP_API_KEY (free at materialsproject.org)."""
    key = api_key or os.environ.get("MP_API_KEY")
    if not key:
        raise DatabaseError("The Materials Project needs an API key: pass api_key= or set MP_API_KEY")
    query = urllib.parse.urlencode({"material_ids": material_id, "_fields": "material_id,structure"})
    payload = json.loads(_get(f"{MP_URL}?{query}", headers={"X-API-KEY": key}, timeout=timeout, opener=opener).decode("utf-8"))
    entries = payload.get("data") or []
    if not entries:
        raise DatabaseError(f"No Materials Project entry {material_id}")
    structure = _structure_from_pymatgen_dict(entries[0]["structure"])
    return conventional_cell(structure) if conventional else structure


if __name__ == "__main__":
    raise SystemExit(main())
