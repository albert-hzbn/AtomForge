"""Runs the C++ and the Rust command lines on the same cases and compares
exit codes, stdout, stderr and every file they write, byte for byte (timings
are masked). Run from the repository root:

    python rust/tools/compare_cli.py [CPP_EXE] [RUST_EXE] [-k FILTER]

Defaults: build-global-safety/AtomForge.exe and rust/target/release/atomforge.exe.
Each case runs in a fresh folder holding copies of the input files; a case may
chain several commands (later ones read what earlier ones wrote)."""
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EX = ROOT / "docs/manual/examples"
IO = ROOT / "rust/tests/fixtures/io"
INPUTS = {
    "cu.cif": EX / "cu_fcc.cif",
    "fe.cif": EX / "fe_bcc.cif",
    "oct.obj": EX / "octahedron.obj",
    "nacl.vasp": IO / "nacl.vasp",
    "hcp.cif": IO / "hcp_number.cif",
    "aspirin.mol": IO / "aspirin.mol",
    "benzene.pdb": IO / "benzene.pdb",
}

BATCH = ('{"tool": "harmonic-thermodynamics", "base": {"energies_eV": [[0.01, 0.02]]}, '
         '"sweep": {"temperatures_K": [[100], [300]]}, "collect": ["heat_capacity_eV_per_K.0"]}')
THERMO = '{"energies_eV": [[0.01, 0.02, 0.05]], "temperatures_K": [100, 300]}'
SYMMETRY = '{"structure": {"file": "cu.cif"}}'


def B(*a):
    return ["--build", *a]


CU8 = B("bulk", "--a", "3.615", "--atom", "Cu 0 0 0", "--output", "cu.vasp")
# (name, [commands], {file: text} written first, stdin text or None)
CASES = [
    ("version", [["--version"], ["-v"]], {}, None),
    ("help", [["--help"], ["-h"], ["--help", "nonsense"], ["--build"], ["--build", "warp"]], {}, None),
    *[("help " + m, [["--help", m]], {}, None) for m in
      ("bulk", "gb", "poly", "nano", "amorphous", "sss", "dislocation", "custom", "interface", "stacking-fault",
       "vacancy", "strain", "primitive", "surface", "sqs", "nanowire", "core-shell")],
    ("bulk fcc", [B("bulk", "--spacegroup", "225", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", "a.cif"),
                  B("bulk", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", "a.xyz"),
                  B("bulk", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", "POSCAR")], {}, None),
    ("bulk others", [B("bulk", "--spacegroup", "227", "--origin", "2", "--a", "5.43", "--atom", "Si 1/8 1/8 1/8", "--output", "si.vasp"),
                     B("bulk", "--spacegroup", "194", "--a", "3.21", "--c", "5.21", "--atom", "Mg 1/3 2/3 1/4", "--output", "mg.cif"),
                     B("bulk", "--system", "cubic", "--spacegroup", "221", "--a", "3.9", "--atom", "Sr 0 0 0", "--atom", "Ti 0.5 0.5 0.5",
                       "--atom", "O 0.5 0.5 0", "--output", "sto.xyz"),
                     B("bulk", "--spacegroup", "62", "--a", "5.2", "--b", "6.1", "--c", "4.4", "--atom", "Fe 0.1 0.25 0.3", "--output", "o.cif")],
     {}, None),
    ("bulk errors", [B("bulk", "--system", "cubik", "--atom", "Cu 0 0 0", "--output", "x.cif"),
                     B("bulk", "--spacegroup", "300", "--atom", "Cu 0 0 0", "--output", "x.cif"),
                     B("bulk", "--system", "hexagonal", "--spacegroup", "225", "--atom", "Cu 0 0 0", "--output", "x.cif"),
                     B("bulk", "--a", "3.6", "--b", "4", "--atom", "Cu 0 0 0", "--output", "x.cif"),
                     B("bulk", "--spacegroup", "166", "--a", "4", "--alpha", "60", "--beta", "60", "--gamma", "60", "--atom", "Bi 0 0 0",
                       "--output", "x.cif"),
                     B("bulk", "--atom", "Cu 0 0", "--output", "x.cif"),
                     B("bulk", "--atom", "Qq 0 0 0", "--output", "x.cif"),
                     B("bulk", "--output", "x.cif"),
                     B("bulk", "--atom", "Cu 0 0 0"),
                     B("bulk", "--spacegrup", "225", "--atom", "Cu 0 0 0", "--output", "x.cif"),
                     B("bulk", "-a", "3", "--atom", "Cu 0 0 0", "--output", "x.cif"),
                     B("bulk", "--a", "abc", "--atom", "Cu 0 0 0", "--output", "x.cif"),
                     B("bulk", "--a", "-3", "--atom", "Cu 0 0 0", "--output", "x.cif")], {}, None),
    ("custom", [B("custom", "--input", "cu.cif", "--mesh", "oct.obj", "--scale", "3", "--output", "c.xyz"),
                B("custom", "--input", "cu.cif", "--mesh", "oct.obj", "--scale", "3", "--miller", "1 1 0", "--no-cell", "--output", "c2.xyz"),
                B("custom", "--input", "cu.cif", "--mesh", "oct.txt", "--output", "c3.xyz"),
                B("custom", "--input", "cu.cif", "--output", "c3.xyz")], {}, None),
    ("sss", [B("sss", "--input", "cu.cif", "--frac", "Cu=0.7,Zn=0.3", "--seed", "42", "--output", "s.cif"),
             B("sss", "--input", "cu.cif", "--frac", "Cu=0.7,Zn=0.2", "--output", "s2.cif"),
             B("sss", "--input", "cu.cif", "--frac", "Cu:0.7", "--output", "s2.cif"),
             B("sss", "--input", "cu.cif", "--frac", "Cu=x", "--output", "s2.cif")], {}, None),
    ("lattice tools", [CU8,
                       B("primitive", "--input", "cu.vasp", "--output", "p.vasp"),
                       B("surface", "--input", "cu.vasp", "--h", "1", "--k", "1", "--l", "1", "--layers", "3", "--output", "s.xyz"),
                       B("sqs", "--input", "cu.vasp", "--element", "Cu 0.5", "--element", "Ni 0.5", "--steps", "200", "--output", "q.xyz"),
                       B("vacancy", "--input", "cu.vasp", "--count", "1", "--output", "v.xyz"),
                       B("vacancy", "--input", "cu.vasp", "--count", "99", "--output", "v2.xyz"),
                       B("strain", "--input", "cu.vasp", "--exx", "0.01", "--eyz", "0.02", "--output", "e.vasp"),
                       B("strain", "--input", "cu.vasp", "--matrix", "1 0 0 0 1 0 0 0 1.1", "--output", "e2.vasp"),
                       B("strain", "--input", "cu.vasp", "--matrix", "1 0 0", "--output", "e3.vasp"),
                       B("nanowire", "--input", "cu.vasp", "--radius", "6", "--output", "w.xyz"),
                       B("core-shell", "--input", "cu.vasp", "--core-radius", "2", "--core-element", "Au", "--shell-element", "Ag",
                         "--output", "cs.xyz"),
                       B("primitive", "--output", "p.vasp")], {}, None),
    ("gb", [B("gb", "--input", "cu.cif", "--sigma", "5", "--axis", "0 0 1", "--output", "gb.cif"),
            B("gb", "--input", "cu.cif", "--sigma", "7", "--axis", "1 1 1", "--plane", "1", "--conventional", "--output", "gb2.xyz"),
            B("gb", "--input", "cu.cif", "--sigma", "4", "--output", "gb3.cif"),
            B("gb", "--input", "cu.cif", "--sigma", "99", "--sigmamax", "20", "--output", "gb3.cif"),
            B("gb", "--input", "cu.cif", "--axis", "1 1", "--output", "gb3.cif")], {}, None),
    ("poly", [B("poly", "--input", "cu.cif", "--grains", "2", "--sizex", "12", "--sizey", "12", "--sizez", "12", "--output", "p.cif"),
              B("poly", "--input", "cu.cif", "--grains", "3", "--sizex", "12", "--sizey", "12", "--sizez", "12", "--euler", "0 0 0",
                "--grain-euler", "3 10 20 30", "--output", "p2.xyz")], {}, None),
    ("nano", [B("nano", "--input", "cu.cif", "--radius", "6", "--output", "n.xyz"),
              B("nano", "--input", "cu.cif", "--shape", "octahedron", "--radius", "7", "--no-cell", "--output", "n2.xyz"),
              B("nano", "--input", "cu.cif", "--shape", "wulff", "--facet", "1 1 1 1.0", "--facet", "1 0 0 1.2", "--radius", "7",
                "--output", "n3.xyz"),
              B("nano", "--input", "cu.cif", "--shape", "box", "--hx", "4", "--hy", "5", "--hz", "6", "--output", "n4.cif"),
              B("nano", "--input", "cu.cif", "--shape", "blob", "--output", "n5.xyz"),
              B("nano", "--input", "cu.cif", "--radius", "0", "--output", "n5.xyz"),
              B("nano", "--input", "cu.cif", "--facet", "1 1 1 1", "--output", "n5.xyz")], {}, None),
    ("amorphous", [B("amorphous", "--element", "Si 32", "--density", "2.33", "--seed", "3", "--output", "a.xyz"),
                   B("amorphous", "--element", "Si 8", "--element", "O 16", "--boxa", "8", "--mindist", "Si O 1.4", "--output", "b.cif"),
                   B("amorphous", "--element", "Si x", "--output", "c.xyz")], {}, None),
    ("dislocation", [B("bulk", "--a", "3.615", "--atom", "Cu 0 0 0", "--output", "cu.vasp"),
                     ["--pipe", "replicate 4 4 4", "--input", "cu.vasp", "--output", "big.vasp", "--quiet"],
                     B("dislocation", "--input", "big.vasp", "--character", "edge", "--output", "d.xyz"),
                     B("dislocation", "--input", "big.vasp", "--character", "screw", "--shape", "cylinder", "--cyl-radius", "6",
                       "--output", "d2.xyz"),
                     B("dislocation", "--input", "big.vasp", "--character", "twisty", "--output", "d3.xyz")], {}, None),
    ("interface", [B("interface", "--layerA", "cu.cif", "--layerB", "nacl.vasp", "--hklA", "1 0 0", "--hklB", "1 0 0", "--max-area", "40",
                     "--layersA", "2", "--layersB", "2", "--output", "i.xyz"),
                   B("interface", "--layerA", "cu.cif", "--layerB", "nacl.vasp", "--hklA", "1 0 0", "--hklB", "1 0 0", "--max-area", "1",
                     "--output", "i2.xyz")], {}, None),
    ("stacking fault", [B("stacking-fault", "--input", "cu.cif", "--plane", "1", "--layers", "6", "--frame", "2", "--sequence", "frames",
                          "--output", "sf.vasp"),
                        B("stacking-fault", "--input", "cu.cif", "--plane", "9", "--output", "sf2.vasp")], {}, None),
    ("analyze", [["--analyze", "cna", "--input", "cu.cif", "--output", "cna.csv"],
                 ["--analyze", "rdf", "--input", "cu.cif", "--output", "rdf.csv", "--bins", "50"],
                 ["--analyze", "adf", "--input", "nacl.vasp", "--output", "adf.csv", "--centre", "Na", "--bins", "30"],
                 ["--analyze", "sro", "--input", "nacl.vasp", "--output", "sro.csv"],
                 ["--analyze", "interstitial", "--input", "cu.cif", "--output", "voids.csv"],
                 ["--analyze", "pattern-init", "--input", "cu.cif", "--cutoff", "3", "--output", "pattern.csv"],
                 ["--analyze", "pattern-detect", "--input", "fe.cif", "--pattern", "pattern.csv", "--cutoff", "3", "--output", "match.csv"],
                 ["--analyze", "sculpt", "--input", "cu.cif", "--slabs", "1 0 0 0 1 1;0 1 0 0 1 1;0 0 1 0 2 1", "--output", "sc.vasp"],
                 ["--analyze", "--help"],
                 ["--analyze", "cna", "--input", "cu.cif", "--output", "x.csv", "--bins", "5"],
                 ["--analyze", "rdf", "--input", "cu.cif", "--output", "x.csv", "--bins", "0"],
                 ["--analyze", "rdf", "--input", "cu.cif", "--output", "x.csv", "--cutoff", "8", "--cutoff", "9"],
                 ["--analyze", "cna", "--input", "cu.cif"],
                 ["--analyze", "magic", "--input", "cu.cif", "--output", "x.csv"],
                 ["--analyze", "cna", "--input", "missing.cif", "--output", "x.csv"]], {}, None),
    ("defect analyses", [B("bulk", "--a", "3.615", "--atom", "Cu 0 0 0", "--output", "cu.vasp"),
                         ["--pipe", "replicate 3 3 3", "--input", "cu.vasp", "--output", "ref.vasp", "--quiet"],
                         ["--pipe", "random-displacement 0.05", "--input", "ref.vasp", "--output", "def.vasp", "--quiet"],
                         ["--analyze", "nye", "--input", "def.vasp", "--reference", "ref.vasp", "--output", "nye.csv"],
                         ["--analyze", "vitek", "--input", "def.vasp", "--reference", "ref.vasp", "--line", "0 0 1", "--burgers", "2.55",
                          "--output", "vitek.csv"],
                         ["--analyze", "fit-dislocation", "--input", "def.vasp", "--reference", "ref.vasp", "--line", "0 0 1",
                          "--output", "fit.csv"],
                         ["--analyze", "drag", "--input", "ref.vasp", "--final", "def.vasp", "--zeta", "0.3", "--output", "mid.vasp"]],
     {}, None),
    ("convert", [["--convert", "--input", "cu.cif", "--output", "cu.xyz"],
                 ["--convert", "--input", "cu.cif", "--output", "cu.extxyz"],
                 ["--convert", "--input", "aspirin.mol", "--output", "a.pdb"],
                 ["--convert", "--input", "benzene.pdb", "--output", "b.mol2"],
                 ["--convert", "--input", "nacl.vasp", "--output", "n.pwi"],
                 ["--convert", "--input", "hcp.cif", "--output", "POSCAR"],
                 ["--convert", "--input", "cu.cif", "--output", "x.lmp", "--format", "lammps"],
                 ["--convert", "--input", "cu.cif", "--output", "x.cif", "--cutoff", "3"]], {}, None),
    ("science", [["--science", "--help"], ["--science"],
                 ["--science", "harmonic-thermodynamics", "--input", "thermo.json", "--output", "t.json", "--report", "t.txt"],
                 ["--science", "harmonic-thermodynamics", "--input", "thermo.json", "--output", "t.json"],
                 ["--science", "symmetry", "--input", "sym.json", "--output", "s.json", "--structures", "cells.extxyz"],
                 ["--science", "symmetry", "--input", "sym.json", "--output", "sym.json"],
                 ["--science", "nonsense", "--input", "sym.json", "--output", "n.json"],
                 ["--science", "symmetry", "--input", "sym.json"],
                 ["--science", "symmetry", "--bogus"],
                 ["--science-batch", "batch.json", "--output", "batch_result.json", "--csv", "batch.csv"]],
     {"thermo.json": THERMO, "sym.json": SYMMETRY, "batch.json": BATCH}, None),
    ("pipe", [["--pipe", "--help"], ["--pipe", "--list"], ["--pipe", "--list-json"],
              ["--pipe", "replicate 2 2 2 | select-random 0.1 | delete-selected", "--input", "cu.cif", "--output", "v.xyz",
               "--report", "r.json", "--save-pipeline", "p.json"],
              ["--pipeline", "p.json", "--input", "cu.cif", "--output", "v2.cif"],
              ["--pipe", "replicate 2 2 2 | select-expression fz > 0.5 | assign-element Ni", "--input", "cu.cif", "--output", "-"],
              ["--pipe", "build-bulk | build-vacancy --count 1 --seed 3", "--input", "cu.cif", "--output", "b.xyz"],
              ["--pipe", "build-sss", "--input", "cu.cif", "--output", "b2.cif"],
              ["--pipe", "replicate 2 2 2 | warp-drive", "--input", "cu.cif", "--output", "x.xyz"],
              ["--pipe", "delete-selected", "--input", "cu.cif", "--output", "x.xyz", "--quite"],
              ["--pipe", "select-element Cu | delete-selected", "--input", "cu.cif", "--output", "x.xyz"]], {}, None),
    ("pipe stdin", [["--pipe", "add-vacuum c 15", "--input", "-", "--output", "out.cif"]], {},
     "4\nLattice=\"3.615 0 0 0 3.615 0 0 0 3.615\" Properties=species:S:1:pos:R:3\nCu 0 0 0\nCu 0 1.8075 1.8075\nCu 1.8075 0 1.8075\n"
     "Cu 1.8075 1.8075 0\n"),
    ("cloud", [["--cloud", "generate", "--lattice", "fcc", "--a", "3.615", "--cells", "4", "4", "4", "--element", "Cu", "--output",
                "c.afcloud", "--quiet", "--threads", "1"],
               ["--cloud", "info", "c.afcloud"],
               ["--cloud", "build", "--input", "cu.cif", "--output", "d.afcloud", "--quiet", "--threads", "1"],
               ["--cloud", "--help"], ["--cloud"], ["--cloud", "warp", "--output", "x.afcloud"]], {}, None),
]

MASKS = [re.compile(rb"\(\d+\.\d ms\)"), re.compile(rb"done in \d+\.\d s \(\d+\.\d M atoms/s\)"),
         re.compile(rb'"milliseconds": [0-9.e+-]+'), re.compile(rb'"elapsed_s": [0-9.e+-]+')]


# Files compared numerically (docs/rust-port-plan.md 6a: Eigen least-squares
# fits are nalgebra in Rust, 1e-9 relative plus 1e-10 absolute).
TOLERANT = {"nye.csv", "fit.csv"}
# Files the C++ writes through Open Babel and the Rust with its native writers
# (header, Gasteiger charges and SYBYL types differ); the switch-over notes list them.
KNOWN = {"a.pdb": "Open Babel PDB writer", "b.mol2": "Open Babel MOL2 writer"}


def close(a, b):
    ta, tb = re.split(rb"[,\s]+", a.strip()), re.split(rb"[,\s]+", b.strip())
    if len(ta) != len(tb):
        return False
    for x, y in zip(ta, tb):
        if x == y:
            continue
        try:
            fx, fy = float(x), float(y)
        except ValueError:
            return False
        if abs(fx - fy) > 1e-10 + 1e-9 * max(abs(fx), abs(fy)):
            return False
    return True


def mask(data):
    for m in MASKS:
        data = m.sub(b"<masked>", data)
    return data


def run_case(exe, folder, commands, files, stdin):
    for name, source in INPUTS.items():
        shutil.copy(source, folder / name)
    for name, text in files.items():
        (folder / name).write_text(text, encoding="utf-8", newline="\n")
    results = []
    env = dict(os.environ)
    env["PATH"] = r"C:\msys64\ucrt64\bin;" + env.get("PATH", "")
    for args in commands:
        done = subprocess.run([str(exe), *args], cwd=folder, capture_output=True, timeout=600, env=env,
                              input=stdin.encode() if stdin is not None else None)
        results.append((done.returncode, mask(done.stdout), mask(done.stderr)))
    written = {}
    for path in sorted(folder.rglob("*")):
        if path.is_file():
            rel = path.relative_to(folder).as_posix()
            if rel in INPUTS or rel in files:
                continue
            written[rel] = mask(path.read_bytes())
    return results, written


def first_difference(a, b):
    n = next((i for i in range(min(len(a), len(b))) if a[i] != b[i]), min(len(a), len(b)))
    return n, a[max(0, n - 60):n + 80], b[max(0, n - 60):n + 80]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-k")]
    filters = [a[2:] or "" for a in sys.argv[1:] if a.startswith("-k")]
    cpp = Path(args[0]) if len(args) > 0 else ROOT / "build-global-safety/AtomForge.exe"
    rust = Path(args[1]) if len(args) > 1 else ROOT / "rust/target/release/atomforge.exe"
    failures = 0
    for name, commands, files, stdin in CASES:
        if filters and not any(f in name for f in filters):
            continue
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            ra, fa = run_case(cpp, Path(a), commands, files, stdin)
            rb, fb = run_case(rust, Path(b), commands, files, stdin)
        problems, known = [], set()
        for k, (x, y) in enumerate(zip(ra, rb)):
            for what, i in (("exit", 0), ("stdout", 1), ("stderr", 2)):
                if x[i] != y[i]:
                    detail = "%s vs %s" % (x[i], y[i]) if what == "exit" else "at %d:\n      C++:  %r\n      Rust: %r" % first_difference(x[i], y[i])
                    problems.append("command %d %s %s differs %s" % (k + 1, " ".join(commands[k])[:70], what, detail))
        for f in sorted(set(fa) | set(fb)):
            if f not in fa or f not in fb:
                problems.append("file %s only from %s" % (f, "C++" if f in fa else "Rust"))
            elif f in KNOWN:
                known.add("%s (%s)" % (f, KNOWN[f]))
            elif fa[f] != fb[f] and not (f in TOLERANT and close(fa[f], fb[f])):
                problems.append("file %s differs at %d:\n      C++:  %r\n      Rust: %r" % ((f,) + first_difference(fa[f], fb[f])))
        status = "ok" if not problems else "DIFF"
        print("%-4s %s (%d commands, %d files)" % (status, name, len(commands), len(fa)))
        for p in problems[:8]:
            print("     " + p)
        for k in sorted(known):
            print("     known difference: " + k)
        failures += bool(problems)
    print("%d cases differ" % failures)
    sys.exit(1 if failures else 0)


main()
