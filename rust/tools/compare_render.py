"""Renders scenes with the C++ and the Rust command lines and compares the
images: the GPUs' rasterisation differs, so the comparison is statistical
(mean absolute difference per channel, the share of pixels that differ by
more than 24 levels, and a coarse 8x8-block comparison that ignores
one-pixel edge shifts). Run from the repository root:

    python rust/tools/compare_render.py [CPP_EXE] [RUST_EXE] [-k FILTER] [--keep DIR]
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
EX = ROOT / "docs/manual/examples"
RF = ROOT / "rust/tests/fixtures/render"
INPUTS = {"cu.cif": EX / "cu_fcc.cif", "si.xyz": RF / "si_super.xyz", "mol.xyz": RF / "molecule.xyz", "nacl.xyz": RF / "nacl.xyz",
          "mg.xyz": RF / "mg_hcp.xyz", "fe.cif": EX / "fe_bcc.cif"}
DISPLAY = ("ambient=0.3\nshadow_strength=0.6\ncontrast=1.1\nsaturation=1.2\nlight_azimuth=120\nlight_elevation=40\n"
           "depth_cueing=0.6\noutline_width=0.25\noutline_colour=0.1 0.1 0.1\nambient_occlusion=0.5\natom_scale=0.8\n"
           "bond_scale=1.5\ncustom_background=1\nbackground=0.1 0.12 0.18\ncell_line_width=2\n")

CASES = [
    ("cu default", ["--input", "cu.cif", "--width", "400", "--height", "300"]),
    ("cu distance", ["--input", "cu.cif", "--width", "320", "--height", "320", "--distance", "14", "--yaw", "30", "--pitch", "20"]),
    ("cu ortho", ["--input", "cu.cif", "--width", "300", "--height", "400", "--orthographic", "--distance", "12", "--roll", "15"]),
    ("si bonds", ["--input", "si.xyz", "--width", "480", "--height", "360"]),
    ("si no bonds no box", ["--input", "si.xyz", "--width", "300", "--height", "300", "--no-bonds", "--no-box"]),
    ("molecule", ["--input", "mol.xyz", "--width", "400", "--height", "300", "--background", "0", "0", "0"]),
    ("nacl colours", ["--input", "nacl.xyz", "--width", "360", "--height", "360", "--color", "Na", "0.2", "0.4", "0.9",
                      "--radius", "Cl", "1.4", "--radius-scale", "0.8"]),
    ("hcp display", ["--input", "mg.xyz", "--width", "400", "--height", "300", "--display", "look.txt"]),
    ("fe display", ["--input", "fe.cif", "--width", "400", "--height", "300", "--display", "look.txt", "--distance", "9"]),
]


def run(exe, folder, args, output):
    env = dict(os.environ)
    env["PATH"] = r"C:\msys64\ucrt64\bin;" + env.get("PATH", "")
    done = subprocess.run([str(exe), "--render", *args, "--output", output], cwd=folder, capture_output=True, text=True, timeout=300, env=env)
    return done


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-k") and not a.startswith("--keep")]
    filters = [a[2:] for a in sys.argv[1:] if a.startswith("-k")]
    keep = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--keep=")), None)
    cpp = Path(args[0]) if args else ROOT / "build-global-safety/AtomForge.exe"
    rust = Path(args[1]) if len(args) > 1 else ROOT / "rust/target/release/atomforge.exe"
    bad = 0
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        for name, source in INPUTS.items():
            shutil.copy(source, folder / name)
        (folder / "look.txt").write_text(DISPLAY)
        for name, case in CASES:
            if filters and not any(f in name for f in filters):
                continue
            stem = name.replace(" ", "_")
            a, b = run(cpp, folder, case, "cpp_%s.png" % stem), run(rust, folder, case, "rust_%s.png" % stem)
            if a.returncode or b.returncode:
                print("FAIL %s: C++ %d %s | Rust %d %s" % (name, a.returncode, a.stderr.strip(), b.returncode, b.stderr.strip()))
                bad += 1
                continue
            x = np.asarray(Image.open(folder / ("cpp_%s.png" % stem)).convert("RGB"), dtype=np.int32)
            y = np.asarray(Image.open(folder / ("rust_%s.png" % stem)).convert("RGB"), dtype=np.int32)
            d = np.abs(x - y)
            mean = d.mean()
            share = (d.max(axis=2) > 24).mean() * 100
            h, w = (x.shape[0] // 8) * 8, (x.shape[1] // 8) * 8
            block = lambda img: img[:h, :w].reshape(h // 8, 8, w // 8, 8, 3).mean(axis=(1, 3))
            blocks = np.abs(block(x) - block(y)).max()
            ok = mean < 3.0 and share < 3.0 and blocks < 40
            bad += not ok
            print("%-4s %-22s mean %.2f  >24: %.2f%%  worst 8x8 block %.1f" % ("ok" if ok else "DIFF", name, mean, share, blocks))
        if keep:
            os.makedirs(keep, exist_ok=True)
            for png in folder.glob("*.png"):
                shutil.copy(png, Path(keep) / png.name)
    print("%d renders differ" % bad)
    sys.exit(1 if bad else 0)


main()
