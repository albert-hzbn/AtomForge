"""Build and plot all seven native nanocrystal cuts for the tutorials."""
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
MANUAL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
if os.environ.get("MANUAL_PYTHON_DEPS"):
    sys.path.insert(0, os.environ["MANUAL_PYTHON_DEPS"])
import atomforge as af
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SHAPES = [
    ("sphere", ["--radius", "10"], "Radius 10 A"),
    ("ellipsoid", ["--rx", "12", "--ry", "8", "--rz", "5"], "Radii 12, 8, 5 A"),
    ("box", ["--hx", "10", "--hy", "7", "--hz", "5"], "Half-lengths 10, 7, 5 A"),
    ("cylinder", ["--cylradius", "7", "--cylheight", "20", "--cylaxis", "2"], "Radius 7 A; height 20 A; Z"),
    ("octahedron", ["--radius", "12"], "Octahedral bound 12 A"),
    ("truncated-octahedron", ["--radius", "15", "--trunc", "8"], "Octahedral bound 15 A; cap 8 A"),
    ("cuboctahedron", ["--radius", "12"], "Pairwise absolute-sum bound 12 A"),
]


def main():
    executable = os.environ.get("ATOMFORGE_PATH", str(ROOT / "build" /
                               ("AtomForge.exe" if os.name == "nt" else "AtomForge")))
    directory = MANUAL / "examples/shapes"
    directory.mkdir(exist_ok=True)
    figures = MANUAL / "figures/tutorials"
    figures.mkdir(exist_ok=True)
    records = []
    for shape, options, dimensions in SHAPES:
        path = directory / (shape + ".xyz")
        command = [executable, "--build", "nano", "--input", str(MANUAL / "examples/cu_fcc.cif"),
                   "--shape", shape, *options, "--vacuum", "5", "--output", str(path)]
        subprocess.run(command, check=True, capture_output=True, text=True)
        atoms = af.load(str(path))
        xyz = np.array([(a.x, a.y, a.z) for a in atoms.atoms])
        if not len(atoms) or not np.isfinite(xyz).all():
            raise RuntimeError(f"Invalid shape output: {shape}")
        fig = plt.figure(figsize=(6.8, 3.8), facecolor="white")
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(*xyz.T, s=16, color="#b36b3f", edgecolors="#724426", linewidths=.15)
        ax.set_box_aspect(np.maximum(np.ptp(xyz, axis=0), 1))
        ax.view_init(elev=22, azim=35)
        ax.set_axis_off()
        ax.set_title(f"{shape.replace('-', ' ').title()} | {len(atoms)} Cu atoms\n{dimensions}", fontsize=12)
        fig.tight_layout()
        fig.savefig(figures / ("shape-" + shape + ".png"), dpi=180, bbox_inches="tight")
        plt.close(fig)
        records.append({"shape": shape, "options": options, "atoms": len(atoms),
                        "file": "shapes/" + path.name, "dimensions": dimensions})
    (MANUAL / "examples/shape-results.json").write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
