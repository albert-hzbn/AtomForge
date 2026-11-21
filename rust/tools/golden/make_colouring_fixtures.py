"""Structures for the grain-boundary colouring comparison
(rust/tests/fixtures/colouring): XYZ files with a "cell ... [origin x y z]"
comment line for periodic ones and an optional fifth column, the grain region.
Run from the repository root: python rust/tools/golden/make_colouring_fixtures.py"""
import os
import random

OUT = "rust/tests/fixtures/colouring/"
os.makedirs(OUT, exist_ok=True)


def write(name, atoms, cell=None, regions=None):
    comment = "cell " + " ".join("%.10g" % v for row in cell for v in row) if cell else name
    lines = [str(len(atoms)), comment]
    for i, (s, p) in enumerate(atoms):
        line = "%s %.10f %.10f %.10f" % (s, *p)
        if regions:
            line += " %d" % regions[i]
        lines.append(line)
    with open(OUT + name + ".xyz", "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


def lattice(basis, a, n, symbol):
    return [(symbol, [(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a]) for i in range(n) for j in range(n) for k in range(n) for b in basis]


fcc = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
bcc = [(0, 0, 0), (.5, .5, .5)]
random.seed(11)
a, n = 3.615, 4
cube = [[a * n, 0, 0], [0, a * n, 0], [0, 0, a * n]]
cu = lattice(fcc, a, n, "Cu")
# A vacancy in the middle of a periodic crystal (coordination seeds).
write("cu_vacancy", [x for k, x in enumerate(cu) if k != 2 * 4 * 16 + 2 * 4 * 4 + 2 * 4], cube)
# A cluster without a cell with a vacancy (bounding-box edges).
write("cu_cluster", [x for k, x in enumerate(lattice(fcc, a, 5, "Cu")) if k != 250], None)
# Perfect periodic crystals (common-neighbour signatures).
write("cu_perfect", cu, cube)
af = 2.87
write("fe_bcc", lattice(bcc, af, 4, "Fe"), [[af * 4, 0, 0], [0, af * 4, 0], [0, 0, af * 4]])
jitter = [(s, [c + random.uniform(-0.25, 0.25) for c in p]) for s, p in cu]
write("cu_jitter", jitter, cube)
# Two grains (regions) in one periodic box.
write("bicrystal", cu, cube, [0 if p[0] < a * n / 2 else 1 for _, p in cu])
# Random packing without a cell, two elements.
write("amorphous", [(random.choice(["Ni", "Al"]), [random.uniform(0, 12) for _ in range(3)]) for _ in range(150)], None)
# Too few atoms.
write("three", [("Cu", [0, 0, 0]), ("Cu", [2.5, 0, 0]), ("Cu", [0, 2.5, 0])], None)
