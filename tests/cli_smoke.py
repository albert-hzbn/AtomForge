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
    for fraction in ("Cu=nan", "Cu=inf", "Cu=0.5garbage"):
        output = root / "invalid.cif"
        run("--build", "sss", "--input", source, "--frac", fraction, "--output", output, success=False)
        assert not output.exists()
    run("--build", "bulk", "--a", "nan", "--output", root / "invalid.cif", success=False)
    run("--build", "unknown", success=False)
print("CLI smoke tests passed (version, 8 help modes, 5 builders, 5 invalid inputs)")
