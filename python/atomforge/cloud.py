"""Atom clouds: structures of up to billions of atoms for AtomForge's large-data viewer.

An ``.afcloud`` file stores atoms by space in chunks (8 bytes per atom, in random
order inside each chunk), which the desktop draws with level of detail and
streaming, so a structure's size is limited by the disk rather than by memory::

    from atomforge import cloud
    cloud.generate("cu1b.afcloud", "fcc", 3.615, (630, 630, 630), "Cu")   # 1.0 billion atoms
    cloud.build("run.dump", "run.afcloud", types=["Cu", "Ni"])            # stream a LAMMPS dump
    info = cloud.info("run.afcloud")
    positions, species = cloud.read("run.afcloud", chunk=0)
    cloud.from_structure(structure, "small.afcloud")

Building and generating run the AtomForge executable (``ATOMFORGE_PATH``);
reading needs only Python (:func:`read` gives NumPy arrays when NumPy is installed)."""

from __future__ import annotations

import json
import os
import struct
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

PathLike = Union[str, os.PathLike]

_MAGIC = b"AFCLOUD1"


@dataclass
class CloudInfo:
    """Description of an atom cloud: atom count, box, species and chunks."""

    atoms: int
    lower: Tuple[float, float, float]
    upper: Tuple[float, float, float]
    species: List[Dict[str, Any]]
    chunks: List[Dict[str, Any]]
    source: str = ""
    cell: Optional[List[List[float]]] = None
    cell_origin: Optional[Tuple[float, float, float]] = None
    path: str = field(default="", repr=False)


def _executable() -> str:
    from .pipeline import _executable as find
    return find()


def _run(arguments: Sequence[Any]) -> None:
    result = subprocess.run([_executable(), "--cloud", *map(str, arguments), "--quiet"], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError((result.stderr or result.stdout).strip().replace("Error: ", "", 1))


def info(path: PathLike) -> CloudInfo:
    """Reads the header of an ``.afcloud`` file (no executable needed)."""
    with open(path, "rb") as handle:
        if handle.read(8) != _MAGIC:
            raise ValueError(f"{path} is not an AtomForge atom cloud")
        offset, length = struct.unpack("<QQ", handle.read(16))
        if not offset:
            raise ValueError(f"{path} is incomplete (its build did not finish)")
        handle.seek(offset)
        header = json.loads(handle.read(length).decode("utf-8"))
    return CloudInfo(
        atoms=int(header["atoms"]), lower=tuple(header["lower"]), upper=tuple(header["upper"]),
        species=header["species"], chunks=header["chunks"], source=header.get("source", ""),
        cell=header.get("cell"), cell_origin=tuple(header["cell_origin"]) if "cell_origin" in header else None,
        path=str(path))


def read(path: PathLike, chunk: int, count: Optional[int] = None, first: int = 0):
    """Atoms of one chunk as ``(positions, species)``: positions in Angstrom and species indices.

    NumPy arrays when NumPy is installed, else lists. A prefix (``count`` from the
    start) is a uniform random subsample of the chunk."""
    from array import array
    import sys
    description = info(path)
    entry = description.chunks[chunk]
    total = int(entry["count"])
    count = total - first if count is None else max(0, min(count, total - first))
    values = array("H")
    with open(path, "rb") as handle:
        handle.seek(int(entry["offset"]) + 8 * first)
        values.frombytes(handle.read(8 * count))
    if sys.byteorder != "little":
        values.byteswap()
    lower = [float(v) for v in entry["lower"]]
    scale = [(float(u) - l) / 65535.0 for u, l in zip(entry["upper"], lower)]
    try:
        import numpy as np
    except ImportError:
        positions = [(lower[0] + values[i] * scale[0], lower[1] + values[i + 1] * scale[1], lower[2] + values[i + 2] * scale[2])
                     for i in range(0, len(values), 4)]
        return positions, list(values[3::4])
    packed = np.asarray(values, dtype=np.uint16).reshape(-1, 4)
    return np.asarray(lower) + packed[:, :3].astype(float) * np.asarray(scale), packed[:, 3].astype(int)


def generate(output: PathLike, lattice: str, a: float, cells: Sequence[int], element: str,
             second: Optional[str] = None, c: Optional[float] = None, chunk: Optional[int] = None) -> CloudInfo:
    """Writes a crystal of ``cells`` (nx, ny, nz) unit cells: lattice fcc, bcc, sc, diamond or hcp.

    ``second`` alternates a second element on the basis (ordered alloy test data)."""
    arguments: List[Any] = ["generate", "--lattice", lattice, "--a", a, "--cells", *cells, "--element", element, "--output", output]
    if second:
        arguments += ["--second", second]
    if c:
        arguments += ["--c", c]
    if chunk:
        arguments += ["--chunk", chunk]
    _run(arguments)
    return info(output)


def build(input: PathLike, output: PathLike, types: Optional[Sequence[str]] = None, frame: int = 0,
          chunk: Optional[int] = None, temporary: Optional[PathLike] = None) -> CloudInfo:
    """Converts a LAMMPS dump or (extended) XYZ file, streaming it, or any other structure file.

    ``types`` names LAMMPS atom types 1, 2, ...; ``frame`` picks a frame of a trajectory."""
    arguments: List[Any] = ["build", "--input", input, "--output", output, "--frame", frame]
    if types:
        arguments += ["--types", ",".join(types)]
    if chunk:
        arguments += ["--chunk", chunk]
    if temporary:
        arguments += ["--temp", temporary]
    _run(arguments)
    return info(output)


def from_structure(structure, output: PathLike, chunk: Optional[int] = None) -> CloudInfo:
    """Writes an :class:`atomforge.Structure` as an atom cloud (through extended XYZ)."""
    with tempfile.TemporaryDirectory(prefix="atomforge_cloud_") as folder:
        source = Path(folder) / "structure.extxyz"
        structure.save(str(source))
        return build(source, output, chunk=chunk)
