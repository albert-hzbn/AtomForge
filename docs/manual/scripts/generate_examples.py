"""Generate manual examples using AtomForge's actual CLI and Python engines.

Run from any directory: python generate_examples.py --chgcar /path/to/CHGCAR
The original electronic file is read only and is not copied into this manual.
"""
import argparse
import hashlib
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
os.environ.setdefault("ATOMFORGE_PATH", str(ROOT / "build" / ("AtomForge.exe" if os.name == "nt" else "AtomForge")))
import atomforge as af
from atomforge.electronic import load_volume
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chgcar", type=Path)
    args = parser.parse_args()
    examples, figures = MANUAL / "examples", MANUAL / "figures"
    examples.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    records = []

    def build(name, *options):
        output = examples / name
        command = [os.environ["ATOMFORGE_PATH"], "--build", *options, "--output", str(output)]
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        structure = af.load(str(output))
        assert len(structure) > 0
        records.append({"file": name, "atoms": len(structure), "options": list(options)})
        (examples / (name + ".log")).write_text(result.stdout + result.stderr, encoding="utf-8")
        return structure

    cu = build("cu_fcc.cif", "bulk", "--system", "cubic", "--spacegroup", "225", "--a", "3.61", "--atom", "Cu 0 0 0")
    fe = build("fe_bcc.cif", "bulk", "--system", "cubic", "--spacegroup", "229", "--a", "2.87", "--atom", "Fe 0 0 0")
    host = cu.repeat(5, 5, 5)
    host.save(str(examples / "cu_host.cif"))
    build("cu_dislocation.cif", "dislocation", "--input", str(examples / "cu_host.cif"),
          "--character", "edge", "--shape", "cylinder", "--cyl-radius", "5",
          "--core", "1.2", "--cutoff", "8")
    alloy = build("cu_ni_alloy.cif", "sss", "--input", str(examples / "cu_host.cif"), "--frac", "Cu=0.7,Ni=0.3", "--seed", "42")
    nano = build("cu_sphere.xyz", "nano", "--input", str(examples / "cu_fcc.cif"), "--shape", "sphere", "--radius", "10", "--vacuum", "5")
    gb = build("cu_sigma5.cif", "gb", "--input", str(examples / "cu_fcc.cif"), "--axis", "0 0 1", "--sigma", "5", "--plane", "0", "--uca", "3", "--ucb", "3", "--overlap", "1.5")
    poly = build("cu_polycrystal.cif", "poly", "--input", str(examples / "cu_fcc.cif"), "--sizex", "25", "--sizey", "25", "--sizez", "25", "--grains", "4", "--seed", "7")
    glass = build("sio2_pack.xyz", "amorphous", "--element", "Si 40", "--element", "O 80", "--density", "1.5", "--seed", "7")
    mesh = examples / "octahedron.obj"
    mesh.write_text("# Unit octahedron; generated example, not an imported asset\n"
                    "v 1 0 0\nv -1 0 0\nv 0 1 0\nv 0 -1 0\nv 0 0 1\nv 0 0 -1\n"
                    "f 1 3 5\nf 3 2 5\nf 2 4 5\nf 4 1 5\nf 3 1 6\nf 2 3 6\nf 4 2 6\nf 1 4 6\n", encoding="ascii")
    custom = build("cu_mesh.xyz", "custom", "--input", str(examples / "cu_fcc.cif"), "--mesh", str(mesh), "--scale", "12", "--vacuum", "5")
    # Keep paths portable in the published provenance file.
    for row in records:
        row["options"] = [str(v).replace(str(examples), "examples") for v in row["options"]]
    (examples / "builder-results.json").write_text(json.dumps(records, indent=2), encoding="utf-8")

    if args.chgcar:
        grid = load_volume(str(args.chgcar)).fields[0]
        values = np.array(grid.values)
        summary = {"sha256": hashlib.sha256(args.chgcar.read_bytes()).hexdigest(), "shape": grid.shape,
                   "cell": grid.cell, "unit": grid.unit, "sites": len(grid.sites), "min": float(values.min()),
                   "max": float(values.max()), "integral": grid.integrate(), "periodic": grid.periodic}
        (examples / "chgcar-summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        cell = np.array(grid.cell)
        section = grid.section(np.array(grid.origin) + .5*cell[2], cell[0], cell[1], (100, 100))
        planar = np.array(grid.planar_average(2))
        cumulative = np.array(grid.cumulative_charge(2))
        np.savetxt(examples / "chgcar-planar.csv", planar, delimiter=",", header="distance_A,density_e_per_A3", comments="")
        np.savetxt(examples / "chgcar-cumulative.csv", cumulative, delimiter=",", header="distance_A,electrons", comments="")
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.3), layout="constrained")
        image = axes[0].imshow(section.section_values, origin="lower", extent=(0, 1, 0, 1), cmap="viridis", aspect="equal")
        axes[0].set(xlabel="fraction along a", ylabel="fraction along b", title="Section at c/2")
        fig.colorbar(image, ax=axes[0], label="e / A$^3$", shrink=.75)
        axes[1].plot(*planar.T, color="#176d83")
        axes[1].set(xlabel="Normal distance (A)", ylabel="Mean density (e / A$^3$)", title="Planar average")
        axes[2].plot(*cumulative.T, color="#a2602d")
        axes[2].set(xlabel="Normal distance (A)", ylabel="Cumulative electrons", title="Integrated profile")
        fig.savefig(figures / "chgcar-analysis.png", dpi=190, bbox_inches="tight")
        plt.close(fig)
        assert abs(cumulative[-1, 1] - summary["integral"]) < 1e-7 * max(1, abs(summary["integral"]))
        zero = grid.density_difference(grid)
        assert max(abs(v) for v in zero.values) == 0
        print(json.dumps(summary, indent=2))
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
