"""Native AtomForge scientific tools from Python, without extra dependencies.

Every tool of ``AtomForge --science`` (structure type, symmetry, dislocation
lines, clusters, voids, diffraction, VASP band structures, phonons,
relaxation, dynamics, NEB, LAMMPS export, ...) runs in the compiled engine
used by the desktop. Parameters are the JSON request keys listed by
:func:`catalog` (see also ``AtomForge --science --catalog``).

>>> import atomforge as af
>>> from atomforge import native
>>> s = af.load("POSCAR")
>>> sym = native.symmetry(structure=s)
>>> sym.result["international_symbol"], sym.result["wyckoff_letters"]
>>> conventional = sym.frames[-1]          # an atomforge.Structure

``native.run("structure-type", structure=s)`` is the explicit form; tool ids
with hyphens are also available as functions with underscores.

Structures are passed as :class:`atomforge.Structure` objects (or the JSON
dicts used in request files); files are referenced as ``{"file": path}``,
relative to ``base``. The engine is the shared library
``atomforge_science_native``; set ``ATOMFORGE_SCIENCE_LIBRARY`` to its path
when it is not found next to the package or the AtomForge executable.
"""

import ctypes as ct
import json
import os
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

from ._structure import Atom, Structure

__all__ = ["NativeResult", "catalog", "run", "tool_ids"]

_library = None
_dll_directories = []


def _load():
    global _library
    if _library is not None:
        return _library
    name = "atomforge_science_native" + (".dll" if sys.platform == "win32" else ".dylib" if sys.platform == "darwin" else ".so")
    package = Path(__file__).resolve().parent
    explicit = os.environ.get("ATOMFORGE_SCIENCE_LIBRARY")
    candidates = [Path(explicit)] if explicit else [
        package / name, package.parent / name, package.parent.parent / name,
        package.parent.parent / "lib" / name, package.parent.parent / "bin" / name,
        package.parent.parent / "build" / name, package.parent.parent / "build" / "Release" / name,
    ]
    if not explicit:
        from ._viewer import _find_atomforge
        viewer = _find_atomforge()
        if viewer:
            root = Path(viewer).resolve().parent
            candidates[:0] = [root / name, root / "lib" / name, root.parent / "lib" / name]
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        raise FileNotFoundError(
            "AtomForge native science library not found. Build target atomforge_science_native "
            "and set ATOMFORGE_SCIENCE_LIBRARY to its full path.")
    if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
        _dll_directories.append(os.add_dll_directory(str(path.resolve().parent)))
    lib = ct.CDLL(str(path.resolve()))
    lib.afs_version.restype = ct.c_int
    lib.afs_free.argtypes = [ct.c_void_p]
    lib.afs_catalog.restype = ct.c_void_p
    lib.afs_run.restype = ct.c_void_p
    lib.afs_run.argtypes = [ct.c_char_p, ct.c_char_p, ct.c_char_p, ct.POINTER(ct.c_int)]
    _library = lib
    return lib


def _take(lib, pointer) -> str:
    if not pointer:
        raise MemoryError("native science call returned no data")
    try:
        return ct.string_at(pointer).decode("utf-8")
    finally:
        lib.afs_free(pointer)


def catalog() -> List[Dict[str, Any]]:
    """Tools with their parameters: ``[{"id", "title", "category", "help", "parameters": [...]}]``."""
    lib = _load()
    return json.loads(_take(lib, lib.afs_catalog()))


def tool_ids() -> List[str]:
    """Identifiers of every native tool, such as ``"symmetry"`` or ``"relax"``."""
    return [tool["id"] for tool in catalog()]


def _structure_json(structure: Structure) -> Dict[str, Any]:
    data: Dict[str, Any] = {"symbols": [a.symbol for a in structure.atoms],
                            "positions": [[a.x, a.y, a.z] for a in structure.atoms]}
    if structure.cell:
        data["cell"] = [list(map(float, row)) for row in structure.cell]  # a cell makes it periodic
    return data


def _from_json(data: Dict[str, Any]) -> Structure:
    structure = Structure()
    for symbol, (x, y, z) in zip(data["symbols"], data["positions"]):
        structure.atoms.append(Atom(symbol, float(x), float(y), float(z)))
    if data.get("cell"):
        structure.cell = [list(map(float, row)) for row in data["cell"]]
    return structure


def _encode(value: Any) -> Any:
    if isinstance(value, Structure):
        return _structure_json(value)
    if isinstance(value, (list, tuple)):
        if value and all(isinstance(item, Structure) for item in value):
            return [_structure_json(item) for item in value]
        return [_encode(item) for item in value]
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    if hasattr(value, "tolist"):  # NumPy arrays and scalars
        return value.tolist()
    if isinstance(value, Path):
        return {"file": str(value)}
    return value


_catalog_cache: Optional[Dict[str, List[str]]] = None


def _parameter_names(tool: str) -> List[str]:
    global _catalog_cache
    if _catalog_cache is None:
        _catalog_cache = {t["id"]: [p["name"] for p in t["parameters"]] for t in catalog()}
    return _catalog_cache.get(tool, [])


def _expand_structures(tool: str, payload: Dict[str, Any]) -> None:
    """Position-based analyses take positions, cell and pbc; accept ``structure=`` for them
    too, and a list of structures as trajectory frames for ``positions``."""
    names = _parameter_names(tool)
    structure = payload.get("structure")
    if isinstance(structure, Structure) and "structure" not in names and "positions" in names:
        del payload["structure"]
        payload["positions"] = structure
    positions = payload.get("positions")
    frames = positions if isinstance(positions, (list, tuple)) and positions and all(isinstance(p, Structure) for p in positions) else None
    if isinstance(positions, Structure) or frames:
        first = frames[0] if frames else positions
        payload["positions"] = ([[[a.x, a.y, a.z] for a in f.atoms] for f in frames] if frames
                                else [[a.x, a.y, a.z] for a in positions.atoms])
        if first.cell and "cell" in names and "cell" not in payload:
            payload["cell"] = first.cell
        if first.cell and "pbc" in names and "pbc" not in payload:
            payload["pbc"] = [True, True, True]


class NativeResult:
    """Result of a native tool: ``result`` (dict), ``frames`` (atomforge structures),
    ``velocities`` (Angstrom/fs per frame) and ``times_fs`` for dynamics."""

    def __init__(self, tool: str, document: Dict[str, Any]):
        self.tool = tool
        self.result: Dict[str, Any] = document["result"]
        self.frames: List[Structure] = [_from_json(frame) for frame in document.get("frames", [])]
        self.velocities: List[List[List[float]]] = document.get("velocities", [])
        self.times_fs: List[float] = document.get("times_fs", [])

    @property
    def files(self) -> Dict[str, str]:
        """Generated text files (for example KPOINTS or in.lammps) by name."""
        return dict(self.result.get("files", {}))

    def write_files(self, directory) -> List[Path]:
        """Write the generated files into ``directory`` (created if needed); returns their paths."""
        folder = Path(directory)
        folder.mkdir(parents=True, exist_ok=True)
        written = []
        for name, text in self.files.items():
            target = folder / name
            target.write_text(text, encoding="utf-8")
            written.append(target)
        return written

    def __getitem__(self, key: str) -> Any:
        return self.result[key]

    def __repr__(self) -> str:
        return f"NativeResult({self.tool!r}, keys={sorted(self.result)[:6]}..., frames={len(self.frames)})"


def run(tool: str, request: Optional[Dict[str, Any]] = None, *, base: Optional[os.PathLike] = None,
        **parameters: Any) -> NativeResult:
    """Run native tool ``tool`` with a request dict and/or keyword parameters.

    Raises ``ValueError`` with the engine's message for invalid requests.
    """
    payload = dict(request or {})
    payload.update(parameters)
    _expand_structures(tool, payload)
    lib = _load()
    status = ct.c_int(1)
    text = _take(lib, lib.afs_run(tool.encode("utf-8"), json.dumps(_encode(payload)).encode("utf-8"),
                                  str(base or ".").encode("utf-8"), ct.byref(status)))
    if status.value != 0:
        raise ValueError(f"{tool}: {text}")
    return NativeResult(tool, json.loads(text))


def __getattr__(name: str):
    # native.structure_type(...) == native.run("structure-type", ...)
    if name.startswith("_"):
        raise AttributeError(name)
    tool = name.replace("_", "-")
    if tool not in tool_ids():
        raise AttributeError(f"atomforge.native has no tool {name!r}")

    def call(request: Optional[Dict[str, Any]] = None, **parameters: Any) -> NativeResult:
        return run(tool, request, **parameters)

    call.__name__ = name
    call.__doc__ = next(t["help"] for t in catalog() if t["id"] == tool)
    return call
