"""Launch the AtomForge GUI to view a Structure."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from ._structure import Structure


def _find_atomforge() -> Optional[str]:
    """Return the path to the AtomForge executable, or None if not found."""
    # 1. Explicit env var
    env = os.environ.get("ATOMFORGE_PATH")
    if env and Path(env).is_file():
        return env

    # 2. Relative to this package file (handles unpacked release: python/ lives next to AtomForge.exe)
    pkg_dir = Path(__file__).resolve().parent
    for candidate in [
        pkg_dir / "AtomForge.exe",
        pkg_dir / "AtomForge",
        pkg_dir.parent / "AtomForge.exe",
        pkg_dir.parent / "AtomForge",
        pkg_dir.parent.parent / "AtomForge.exe",
        pkg_dir.parent.parent / "AtomForge",
        pkg_dir.parent.parent / "build" / "Release" / "AtomForge.exe",
        pkg_dir.parent.parent / "build" / "AtomForge.exe",
        pkg_dir.parent.parent / "build" / "AtomForge",
    ]:
        if candidate.is_file():
            return str(candidate)

    # 3. System PATH
    import shutil
    found = shutil.which("AtomForge") or shutil.which("atomforge")
    if found:
        return found

    # 4. Common Windows install paths
    if sys.platform == "win32":
        program_files = os.environ.get("ProgramFiles", r"C:\Program Files")
        for root in [program_files, os.environ.get("ProgramFiles(x86)", "")]:
            candidate = Path(root) / "AtomForge" / "AtomForge.exe"
            if candidate.is_file():
                return str(candidate)

    return None


def view(s: "Structure") -> None:
    """
    Open *s* in the AtomForge GUI (non-blocking).

    The structure is serialised to a temporary extXYZ file and AtomForge is
    launched as a detached subprocess.  The temp file is cleaned up after
    AtomForge exits (on POSIX) or left to the OS temp-file sweeper (Windows),
    because AtomForge must be able to read the file after this function returns.
    """
    exe = _find_atomforge()
    if exe is None:
        raise FileNotFoundError(
            "AtomForge executable not found.  Set the ATOMFORGE_PATH environment "
            "variable to the full path of AtomForge.exe, or ensure it is on PATH."
        )

    # Write to a named temp file that survives past this function
    suffix = ".xyz"
    tmp = tempfile.NamedTemporaryFile(
        prefix="atomforge_view_", suffix=suffix, delete=False
    )
    tmp_path = tmp.name
    tmp.close()

    from ._io import _save_xyz
    _save_xyz(s, tmp_path)

    if sys.platform == "win32":
        # DETACHED_PROCESS: new console, parent death does not kill child
        DETACHED_PROCESS = 0x00000008
        CREATE_NEW_PROCESS_GROUP = 0x00000200
        subprocess.Popen(
            [exe, tmp_path],
            creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
            close_fds=True,
        )
    else:
        subprocess.Popen(
            [exe, tmp_path],
            start_new_session=True,
            close_fds=True,
        )
