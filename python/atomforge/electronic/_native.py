"""Lazy ctypes boundary for the shared AtomForge electronic calculation engine."""

import ctypes as ct
import os
from pathlib import Path
import sys


_library = None
_dll_directories = []
Doubles = ct.POINTER(ct.c_double)
Ints = ct.POINTER(ct.c_int)


def library():
    global _library
    if _library is not None:
        return _library
    name = "atomforge_electronic" + (
        ".dll" if sys.platform == "win32" else ".dylib" if sys.platform == "darwin" else ".so"
    )
    package = Path(__file__).resolve().parents[1]
    explicit = os.environ.get("ATOMFORGE_ELECTRONIC_LIBRARY")
    candidates = [Path(explicit)] if explicit else [
        package / name, package.parent / name, package.parent.parent / name,
        package.parent.parent / "lib" / name,
        package.parent.parent / "bin" / name,
        package.parent.parent / "build" / name,
        package.parent.parent / "build" / "Release" / name,
    ]
    viewer = os.environ.get("ATOMFORGE_PATH")
    if viewer and not explicit:
        candidates.insert(0, Path(viewer).resolve().parent / name)
        candidates.insert(1, Path(viewer).resolve().parent / "lib" / name)
        candidates.insert(2, Path(viewer).resolve().parent.parent / "lib" / name)
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        raise FileNotFoundError(
            "AtomForge electronic library not found. Build target atomforge_electronic "
            "and set ATOMFORGE_ELECTRONIC_LIBRARY to its full path. "
            "The standard structure APIs do not require this library."
        )
    if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
        _dll_directories.append(os.add_dll_directory(str(path.resolve().parent)))
    lib = ct.CDLL(str(path.resolve()))
    signatures = {
        "version": (ct.c_int, []),
        "error": (ct.c_char_p, []),
        "free": (None, [ct.c_void_p]),
        "load": (ct.c_void_p, [ct.c_char_p, ct.c_char_p, ct.c_char_p]),
        "create": (ct.c_void_p, [Ints, Doubles, Doubles, ct.c_int, Doubles, ct.c_size_t, ct.c_char_p]),
        "select": (ct.c_void_p, [ct.c_void_p, ct.c_int]),
        "fields": (ct.c_int, [ct.c_void_p]),
        "label": (ct.c_char_p, [ct.c_void_p, ct.c_int, ct.c_int]),
        "data": (ct.c_size_t, [ct.c_void_p, Doubles]),
        "geometry": (ct.c_int, [ct.c_void_p, Ints, Doubles, Doubles]),
        "columns": (ct.c_int, [ct.c_void_p]),
        "sites": (ct.c_size_t, [ct.c_void_p, Doubles]),
        "with_sites": (ct.c_void_p, [ct.c_void_p, Doubles, ct.c_size_t]),
        "save": (ct.c_int, [ct.c_void_p, ct.c_char_p, ct.c_char_p]),
        "calculate": (ct.c_void_p, [ct.c_void_p, ct.c_void_p, ct.c_char_p, Doubles, ct.c_size_t]),
    }
    for key, (result, args) in signatures.items():
        function = getattr(lib, "af_electronic_" + key)
        function.restype = result
        function.argtypes = args
    if lib.af_electronic_version() != 1:
        raise RuntimeError("Incompatible AtomForge electronic library ABI; rebuild/update the library")
    _library = lib
    return lib


def doubles(values):
    values = list(values)
    return (ct.c_double * len(values))(*values)


class Handle:
    def __init__(self, pointer):
        self.lib = library()
        if not pointer:
            raise ValueError(self.lib.af_electronic_error().decode("utf-8", errors="replace"))
        self.pointer = pointer

    def __del__(self):
        pointer = getattr(self, "pointer", None)
        if pointer:
            self.lib.af_electronic_free(pointer)
            self.pointer = None

    def calculate(self, operation, parameters=(), other=None):
        data = doubles(parameters)
        return Handle(self.lib.af_electronic_calculate(
            self.pointer, other.pointer if other is not None else None,
            operation.encode("ascii"), data, len(data),
        ))

    def data(self):
        count = self.lib.af_electronic_data(self.pointer, None)
        data = (ct.c_double * count)()
        self.lib.af_electronic_data(self.pointer, data)
        return list(data)

    def rows(self):
        columns = self.lib.af_electronic_columns(self.pointer)
        if not columns:
            raise TypeError("Result is a grid, not a table")
        data = self.data()
        return [tuple(data[i:i + columns]) for i in range(0, len(data), columns)]

    def save(self, path, format):
        if not self.lib.af_electronic_save(self.pointer, str(path).encode("utf-8"), format.encode("ascii")):
            raise ValueError(self.lib.af_electronic_error().decode("utf-8", errors="replace"))
