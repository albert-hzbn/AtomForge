"""Structures from online databases: the Crystallography Open Database (COD)
and the Materials Project.

    import atomforge as af
    from atomforge.databases import fetch_cod, fetch_materials_project
    rutile = fetch_cod(9015662)                      # open access
    iron = fetch_materials_project("mp-13")           # needs an API key (MP_API_KEY)

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
    """python -m atomforge.databases cod 9008460 out.cif | mp mp-13 out.xyz"""
    import argparse
    from ._io import save
    parser = argparse.ArgumentParser(prog="python -m atomforge.databases", description="Download a structure from COD or the Materials Project.")
    parser.add_argument("database", choices=["cod", "mp"])
    parser.add_argument("id")
    parser.add_argument("output", help="structure file to write (format from the extension)")
    parser.add_argument("--api-key", help="Materials Project API key (default: MP_API_KEY)")
    args = parser.parse_args(argv)
    try:
        structure = fetch_cod(args.id) if args.database == "cod" else fetch_materials_project(args.id, api_key=args.api_key)
    except (DatabaseError, ValueError) as error:
        print(f"Error: {error}")
        return 1
    save(structure, args.output)
    print(f"Saved {len(structure.atoms)} atoms to {args.output}")
    return 0


def fetch_materials_project(material_id: str, *, api_key: Optional[str] = None, timeout: float = 30.0,
                            opener: Optional[Callable] = None) -> Structure:
    """A Materials Project structure (e.g. "mp-149") with its magnetic moments
    where the database gives them. The API key comes from ``api_key`` or the
    environment variable MP_API_KEY (free at materialsproject.org)."""
    key = api_key or os.environ.get("MP_API_KEY")
    if not key:
        raise DatabaseError("The Materials Project needs an API key: pass api_key= or set MP_API_KEY")
    query = urllib.parse.urlencode({"material_ids": material_id, "_fields": "material_id,structure"})
    payload = json.loads(_get(f"{MP_URL}?{query}", headers={"X-API-KEY": key}, timeout=timeout, opener=opener).decode("utf-8"))
    entries = payload.get("data") or []
    if not entries:
        raise DatabaseError(f"No Materials Project entry {material_id}")
    return _structure_from_pymatgen_dict(entries[0]["structure"])


if __name__ == "__main__":
    raise SystemExit(main())
