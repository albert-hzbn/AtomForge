"""Prepare structures for real AtomForge cover screenshots.

Run generate_examples.py first. Open these files in the desktop application;
this script does not render or alter screenshot pixels or atomic radii.
"""
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
sys.path.insert(0, str(ROOT / "python"))
import atomforge as af


def main():
    af.load(EXAMPLES / "cu_fcc.cif").repeat(2, 2, 2).save(EXAMPLES / "cover-bulk.cif")
    af.load(EXAMPLES / "cu_sigma5.cif").repeat(1, 8, 6).save(EXAMPLES / "cover-gb.cif")
    # An illustrative isolated octahedral coordination environment, not a crystal.
    cluster = af.Structure()
    cluster.add_atom("Ti", 0, 0, 0)
    for axis in range(3):
        for sign in (-1, 1):
            position = [0, 0, 0]
            position[axis] = sign * 2
            cluster.add_atom("O", *position)
    cluster.save(EXAMPLES / "cover-tio6.xyz")
    executable = os.environ.get("ATOMFORGE_PATH", str(ROOT / "build" /
                                ("AtomForge.exe" if os.name == "nt" else "AtomForge")))
    af.wulff(EXAMPLES / "cu_fcc.cif", [(1, 0, 0, 1), (1, 1, 1, 1)], radius=20,
             output=EXAMPLES / "cover-wulff.xyz", executable=executable)


if __name__ == "__main__":
    main()
