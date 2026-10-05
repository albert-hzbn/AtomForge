"""Headless integration tests: python tests/cli_smoke.py path/to/AtomForge."""
from pathlib import Path
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())


def run(*args, success=True):
    result = subprocess.run([exe, *map(str, args)], capture_output=True, text=True, timeout=60)
    if (result.returncode == 0) != success:
        raise AssertionError(f"{args}: exit {result.returncode}\n{result.stdout}\n{result.stderr}")
    return result


with tempfile.TemporaryDirectory(prefix="atomforge_cli_") as folder:
    root = Path(folder)
    run("--version")
    for mode in ("bulk", "gb", "poly", "nano", "amorphous", "sss", "dislocation", "custom"):
        run("--help", mode)
    source = root / "cu.cif"
    run("--build", "bulk", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", source)
    assert source.stat().st_size > 0
    for mode, options in (
        ("nano", ["--radius", "5"]),
        ("poly", ["--grains", "2", "--sizex", "10", "--sizey", "10", "--sizez", "10"]),
        ("sss", ["--frac", "Cu=0.5,Ni=0.5"]),
        ("gb", ["--sigma", "5", "--axis", "0 0 1"]),
    ):
        output = root / f"{mode}.cif"
        run("--build", mode, "--input", source, *options, "--output", output)
        assert output.stat().st_size > 0
    # Defect, lattice and nanostructure builders on a conventional fcc cell.
    conventional = root / "cu_conv.vasp"
    run("--build", "bulk", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", conventional)
    for mode in ("vacancy", "strain", "primitive", "surface", "sqs", "nanowire", "core-shell"):
        assert mode.upper().replace("-", "") in run("--help", mode).stdout.upper().replace("-", "")
    for mode, options in (
        ("vacancy", ["--count", "1"]),
        ("strain", ["--exx", "0.01"]),
        ("primitive", []),
        ("surface", ["--h", "1", "--k", "1", "--l", "1", "--layers", "3"]),
        ("sqs", ["--element", "Cu 0.5", "--element", "Ni 0.5", "--steps", "200"]),
        ("nanowire", ["--radius", "6"]),
        ("core-shell", ["--core-radius", "2", "--core-element", "Au", "--shell-element", "Ag"]),
    ):
        output = root / f"{mode}.xyz"
        run("--build", mode, "--input", conventional, *options, "--output", output)
        lines = output.read_text().splitlines()
        assert int(lines[0]) > 0 and len(lines) >= int(lines[0]) + 2, mode
    run("--build", "vacancy", "--input", conventional, "--count", "99", "--output", root / "too_many.xyz", success=False)
    for fraction in ("Cu=nan", "Cu=inf", "Cu=0.5garbage"):
        output = root / "invalid.cif"
        run("--build", "sss", "--input", source, "--frac", fraction, "--output", output, success=False)
        assert not output.exists()
    run("--build", "bulk", "--a", "nan", "--output", root / "invalid.cif", success=False)
    run("--build", "unknown", success=False)
print("CLI smoke tests passed (version, 15 help modes, 12 builders, 6 invalid inputs)")
