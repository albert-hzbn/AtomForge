"""Structures for the af-render golden comparison (rust/tests/fixtures/render):
XYZ files with a "cell a1 a2 a3 b1 b2 b3 c1 c2 c3 [origin x y z]" comment
line for periodic ones. Run from the repository root:
python rust/tools/golden/make_render_fixtures.py"""
import math
import os
import random

OUT = "rust/tests/fixtures/render/"
os.makedirs(OUT, exist_ok=True)


def write(name, atoms, cell=None, origin=None):
    comment = name
    if cell:
        comment = "cell " + " ".join("%.10g" % v for row in cell for v in row)
        if origin:
            comment += " origin " + " ".join("%.10g" % v for v in origin)
    lines = [str(len(atoms)), comment] + ["%s %.10f %.10f %.10f" % (s, *p) for s, p in atoms]
    with open(OUT + name + ".xyz", "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


a = 3.615
fcc = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
write("cu_fcc", [("Cu", [f[0] * a, f[1] * a, f[2] * a]) for f in fcc], [[a, 0, 0], [0, a, 0], [0, 0, a]])
s = 5.431
diamond = [p for f in fcc for p in (f, (f[0] + .25, f[1] + .25, f[2] + .25))]
write("si_diamond", [("Si", [f[0] * s, f[1] * s, f[2] * s]) for f in diamond], [[s, 0, 0], [0, s, 0], [0, 0, s]])
n = 5.64
nacl = [atom for f in fcc for atom in (("Na", [f[0] * n, f[1] * n, f[2] * n]), ("Cl", [(f[0] + .5) * n, f[1] * n, f[2] * n]))]
write("nacl", nacl, [[n, 0, 0], [0, n, 0], [0, 0, n]], origin=[0.3, -0.2, 0.1])
# Hexagonal (non-orthogonal) cell with atoms on the cell faces.
ah, ch = 3.21, 5.21
hexcell = [[ah, 0, 0], [-ah / 2, ah * math.sqrt(3) / 2, 0], [0, 0, ch]]
def frac(f):
    return [sum(f[k] * hexcell[k][i] for k in range(3)) for i in range(3)]
write("mg_hcp", [("Mg", frac([1 / 3, 2 / 3, 0.25])), ("Mg", frac([2 / 3, 1 / 3, 0.75])), ("O", frac([0, 0, 0])), ("O", frac([0.5, 0, 0.5]))], hexcell)
# A molecule without a cell (water, methane and a C-C chain).
random.seed(5)
molecule = [("C", [0, 0, 0]), ("C", [1.54, 0, 0]), ("C", [2.05, 1.45, 0]), ("O", [3.45, 1.5, 0.1]), ("H", [-0.36, -1.03, 0]),
            ("H", [-0.36, 0.51, 0.89]), ("H", [-0.36, 0.51, -0.89]), ("H", [1.9, -0.5, 0.89]), ("H", [1.9, -0.5, -0.89]),
            ("H", [1.7, 2.0, 0.89]), ("H", [3.8, 2.4, 0.1]), ("N", [0.2, 3.1, -0.4]), ("H", [0.9, 3.5, 0.2])]
write("molecule", molecule)
# A distorted Si supercell (2x2x2, 64 atoms) for many bonds.
big = []
for i in range(2):
    for j in range(2):
        for k in range(2):
            for f in diamond:
                big.append(("Si", [(f[0] + i) * s + random.uniform(-0.1, 0.1), (f[1] + j) * s + random.uniform(-0.1, 0.1),
                                   (f[2] + k) * s + random.uniform(-0.1, 0.1)]))
write("si_super", big, [[2 * s, 0, 0], [0, 2 * s, 0], [0, 0, 2 * s]])
print("render fixtures written")
