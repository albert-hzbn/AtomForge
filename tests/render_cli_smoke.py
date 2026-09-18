"""Headless --render CLI test: python tests/render_cli_smoke.py path/to/AtomForge.

Exits 77 (ctest SKIP_RETURN_CODE) if no GL context/display is available,
mirroring electronic_viewport's GPU-less skip convention.
"""
from pathlib import Path
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())

with tempfile.TemporaryDirectory(prefix="atomforge_render_") as folder:
    root = Path(folder)
    source = root / "cu.cif"
    subprocess.run(
        [exe, "--build", "bulk", "--a", "3.61", "--atom", "Cu 0 0 0", "--output", str(source)],
        check=True, capture_output=True, text=True, timeout=60)

    help_result = subprocess.run([exe, "--render", "--help"], capture_output=True, text=True, timeout=30)
    if help_result.returncode != 0 or "--render" not in help_result.stdout:
        raise AssertionError(f"--render --help failed:\n{help_result.stdout}\n{help_result.stderr}")

    output = root / "cu.png"
    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(output),
         "--width", "64", "--height", "64", "--color", "Cu 0.9 0.5 0.1"],
        capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        if "OpenGL context" in result.stderr:
            print("Skipping: no GPU/display available for --render:", result.stderr.strip())
            sys.exit(77)
        raise AssertionError(f"--render failed: exit {result.returncode}\n{result.stdout}\n{result.stderr}")

    data = output.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError("Rendered file is not a valid PNG")
    if len(data) < 100:
        raise AssertionError("Rendered PNG is suspiciously small")

    # Turntable mode writes numbered frames instead of a single file.
    turntable_base = root / "turn.png"
    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(turntable_base),
         "--width", "64", "--height", "64", "--frames", "2", "--yaw-step", "45"],
        capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        raise AssertionError(f"--render turntable failed: exit {result.returncode}\n{result.stdout}\n{result.stderr}")
    for index in range(2):
        frame_path = root / f"turn-{index:03d}.png"
        if not frame_path.is_file() or frame_path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
            raise AssertionError(f"Missing or invalid turntable frame {frame_path}")

    # Invalid element symbols must be rejected, not silently ignored.
    result = subprocess.run(
        [exe, "--render", "--input", str(source), "--output", str(root / "bad.png"),
         "--color", "Xx 0.5 0.5 0.5"],
        capture_output=True, text=True, timeout=30)
    if result.returncode == 0:
        raise AssertionError("--render accepted an unknown element symbol in --color")

print("Render CLI smoke tests passed (help, single PNG, turntable frames, invalid input)")
