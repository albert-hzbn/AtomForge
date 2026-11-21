#!/usr/bin/env python3
"""Builds the portable AtomForge package (replaces the C++ CPack packages).

    python rust/tools/package.py [--no-build] [--out DIR]

Layout, as the C++ portable packages:

  Windows  AtomForge/AtomForge.exe and the native libraries (DLLs) beside it
  Linux    AtomForge/bin/AtomForge and the native libraries (.so) beside it

and in both AtomForge/python/atomforge (the Python package),
AtomForge/README.md, AtomForge/AtomForge-manual.pdf and AtomForge/examples.
The archive is AtomForge-<version>-win64.zip or AtomForge-<version>-linux-x86_64.tar.gz
in rust/target/package (or --out). AtomForge is one program: the desktop,
and the command line when given a mode such as --build or --pipe.
"""
import argparse
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUST = ROOT / "rust"
EXAMPLES = ["cu_fcc.cif", "fe_bcc.cif", "cu_ni_alloy.cif", "cu_sigma5.cif", "cu_polycrystal.cif",
            "cu_dislocation.cif", "cu_sphere.xyz", "sio2_pack.xyz", "cu_host.cif", "octahedron.obj"]
WINDOWS = sys.platform == "win32"


def version():
    text = (RUST / "crates" / "af-cli" / "src" / "lib.rs").read_text(encoding="utf-8")
    m = re.search(r'VERSION:\s*&str\s*=\s*"([^"]+)"', text)
    if m:
        return m.group(1)
    text = (RUST / "Cargo.toml").read_text(encoding="utf-8")
    return re.search(r'^version\s*=\s*"([^"]+)"', text, re.M).group(1)


def build():
    cargo = shutil.which("cargo") or str(Path.home() / ".cargo" / "bin" / ("cargo.exe" if WINDOWS else "cargo"))
    subprocess.run([cargo, "build", "--release", "-p", "af-app", "-p", "af-ffi-science", "-p", "af-ffi-electronic"], cwd=RUST, check=True)


def native_files(release: Path):
    """(source, name in the package) of the program and its libraries."""
    if WINDOWS:
        return [(release / "AtomForge.exe", "AtomForge.exe"),
                (release / "atomforge_science_native.dll", "atomforge_science_native.dll"),
                (release / "atomforge_electronic.dll", "atomforge_electronic.dll"),
                (RUST / "crates" / "af-spglib" / "vendor" / "lib" / "symspg.dll", "symspg.dll")]
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    # Cargo names cdylibs lib*.so; the Python package looks for the C++ names.
    return [(release / "AtomForge", "AtomForge"),
            (release / f"libatomforge_science_native{suffix}", f"atomforge_science_native{suffix}"),
            (release / f"libatomforge_electronic{suffix}", f"atomforge_electronic{suffix}")]


def system_libraries(binary: Path):
    """Linux: the non-system shared libraries a binary needs (spglib), via ldd."""
    if WINDOWS or sys.platform == "darwin" or not shutil.which("ldd"):
        return []
    out = subprocess.run(["ldd", str(binary)], capture_output=True, text=True).stdout
    keep = []
    for line in out.splitlines():
        m = re.match(r"\s*(\S+)\s*=>\s*(\S+)", line)
        if m and "symspg" in m.group(1):
            keep.append(Path(m.group(2)))
    return keep


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-build", action="store_true", help="package the existing release build")
    parser.add_argument("--out", default=str(RUST / "target" / "package"))
    args = parser.parse_args()
    if not args.no_build:
        build()
    release = RUST / "target" / "release"
    out = Path(args.out)
    stage = out / "AtomForge"
    if stage.exists():
        shutil.rmtree(stage)
    bin_dir = stage if WINDOWS else stage / "bin"
    bin_dir.mkdir(parents=True)
    for source, name in native_files(release):
        if not source.is_file():
            sys.exit(f"missing {source}: build first (or drop --no-build)")
        shutil.copy2(source, bin_dir / name)
    for library in system_libraries(release / "AtomForge"):
        shutil.copy2(library, bin_dir / library.name)
    if not WINDOWS:
        # The program finds its libraries beside it.
        (bin_dir / "AtomForge").chmod(0o755)
    shutil.copytree(ROOT / "python" / "atomforge", stage / "python" / "atomforge",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copy2(ROOT / "README.md", stage / "README.md")
    shutil.copy2(ROOT / "docs" / "manual" / "AtomForge-manual.pdf", stage / "AtomForge-manual.pdf")
    (stage / "examples").mkdir()
    for name in EXAMPLES:
        shutil.copy2(ROOT / "docs" / "manual" / "examples" / name, stage / "examples" / name)
    v = version()
    if WINDOWS:
        archive = out / f"AtomForge-{v}-win64.zip"
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
            for path in sorted(stage.rglob("*")):
                z.write(path, path.relative_to(out))
    else:
        machine = platform.machine() or "x86_64"
        archive = out / f"AtomForge-{v}-{'macos' if sys.platform == 'darwin' else 'linux'}-{machine}.tar.gz"
        with tarfile.open(archive, "w:gz") as t:
            t.add(stage, arcname="AtomForge")
    size = archive.stat().st_size / 1e6
    print(f"{archive} ({size:.1f} MB)")


if __name__ == "__main__":
    main()
