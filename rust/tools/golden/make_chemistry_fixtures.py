"""Molecules for the af-chemistry golden comparison (rust/tests/fixtures/chemistry):
plain XYZ files (one molecule each; a "cell" comment line gives a periodic
cell) from the io fixtures and hand-built geometries with known stereo,
charges and periodic bonds. Run from the repository root:
python rust/tools/golden/make_chemistry_fixtures.py"""
import math
import os

IO = "rust/tests/fixtures/io/"
OUT = "rust/tests/fixtures/chemistry/"
os.makedirs(OUT, exist_ok=True)


def write(name, atoms, cell=None):
    comment = "cell " + " ".join("%.6f" % v for row in cell for v in row) if cell else name
    lines = [str(len(atoms)), comment] + ["%s %.6f %.6f %.6f" % (s, *p) for s, p in atoms]
    with open(OUT + name + ".xyz", "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


def mol_v2000(path):
    lines = open(path).read().splitlines()
    n = int(lines[3][:3])
    return [(l[31:34].strip(), [float(l[0:10]), float(l[10:20]), float(l[20:30])]) for l in lines[4:4 + n]]


def mol2(path):
    lines = open(path).read().splitlines()
    start = lines.index("@<TRIPOS>ATOM") + 1
    out = []
    for l in lines[start:]:
        if l.startswith("@"):
            break
        w = l.split()
        out.append((w[5].split(".")[0], [float(w[2]), float(w[3]), float(w[4])]))
    return out


def pdb(path):
    return [(l[76:78].strip(), [float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def gjf(path):
    lines = open(path).read().splitlines()
    out = []
    for l in lines[6:]:
        w = l.split()
        if len(w) != 4:
            break
        out.append((w[0], [float(v) for v in w[1:]]))
    return out


write("aspirin", mol_v2000(IO + "aspirin.mol"))
write("caffeine", mol2(IO + "caffeine.mol2"))
write("benzene", pdb(IO + "benzene.pdb"))
write("ethanol", gjf(IO + "ethanol.gjf"))

# Tetrahedral CHFClBr (one enantiomer) and its mirror image.
t = 1.0 / math.sqrt(3)
dirs = [[t, t, t], [t, -t, -t], [-t, t, -t], [-t, -t, t]]
lengths = {"H": 1.09, "F": 1.35, "Cl": 1.77, "Br": 1.94}
chfclbr = [("C", [0, 0, 0])] + [(s, [d * lengths[s] for d in dirs[k]]) for k, s in enumerate(["H", "F", "Cl", "Br"])]
write("chfclbr", chfclbr)
write("chfclbr_mirror", [(s, [-p[0], p[1], p[2]]) for s, p in chfclbr])


# 2-Butene, E and Z (planar sp2 carbons, methyl hydrogens staggered).
def butene(cis):
    c1, c2 = [0, 0, 0], [1.34, 0, 0]
    a = math.radians(120)
    m1 = [c1[0] + 1.50 * math.cos(a), 1.50 * math.sin(a), 0]
    m2 = [c2[0] - 1.50 * math.cos(a), (1 if cis else -1) * 1.50 * math.sin(a), 0]
    h1 = [c1[0] + 1.09 * math.cos(-a), 1.09 * math.sin(-a), 0]
    h2 = [c2[0] - 1.09 * math.cos(-a), (1 if cis else -1) * 1.09 * math.sin(-a), 0]
    atoms = [("C", c1), ("C", c2), ("C", m1), ("C", m2), ("H", h1), ("H", h2)]
    for m, c in ((m1, c1), (m2, c2)):
        axis = [m[k] - c[k] for k in range(3)]
        n = math.sqrt(sum(v * v for v in axis))
        axis = [v / n for v in axis]
        perp = [0, 0, 1]
        other = [axis[1] * perp[2] - axis[2] * perp[1], axis[2] * perp[0] - axis[0] * perp[2], axis[0] * perp[1] - axis[1] * perp[0]]
        for k in range(3):
            phi = 2 * math.pi * k / 3
            d = [axis[i] / 3 + math.sqrt(8 / 9) * (math.cos(phi) * perp[i] + math.sin(phi) * other[i]) for i in range(3)]
            atoms.append(("H", [m[i] + 1.09 * d[i] for i in range(3)]))
    return atoms


write("z_butene", butene(True))
write("e_butene", butene(False))
# Ammonium acetate: N with four H (N+), acetate O- by valence.
write("ammonium_acetate", [
    ("N", [5, 0, 0]), ("H", [5.6, 0.6, 0.6]), ("H", [5.6, -0.6, -0.6]), ("H", [4.4, 0.6, -0.6]), ("H", [4.4, -0.6, 0.6]),
    ("C", [0, 0, 0]), ("C", [1.52, 0, 0]), ("O", [2.15, 1.07, 0]), ("O", [2.15, -1.07, 0]),
    ("H", [-0.36, 1.03, 0]), ("H", [-0.36, -0.51, 0.89]), ("H", [-0.36, -0.51, -0.89])])
# A periodic polyethylene chain along c (two CH2 per cell, bonds across the boundary).
cpc = 2.55
pe = [("C", [0.0, 0.0, 0.0]), ("C", [0.0, 0.85, 1.275]),
      ("H", [0.89, -0.63, 0.0]), ("H", [-0.89, -0.63, 0.0]), ("H", [0.89, 1.48, 1.275]), ("H", [-0.89, 1.48, 1.275])]
write("polyethylene", [(s, [p[0] + 3, p[1] + 3, p[2] + 0.4]) for s, p in pe], cell=[[8, 0, 0], [0, 8, 0], [0, 0, cpc]])

# Point groups: symmetric molecules (Td, Oh, C2v, C3v, D3h, D3d, C2h, D2d,
# Dinfh, Cinfv, D6h) and a slightly distorted benzene for symmetrisation.
write("methane", [("C", [0, 0, 0])] + [("H", [d * 1.09 for d in dirs[k]]) for k in range(4)])
write("sf6", [("S", [0, 0, 0])] + [("F", v) for v in ([1.56, 0, 0], [-1.56, 0, 0], [0, 1.56, 0], [0, -1.56, 0], [0, 0, 1.56], [0, 0, -1.56])])
w = math.radians(104.5 / 2)
write("water", [("O", [0, 0, 0]), ("H", [0.96 * math.sin(w), 0, 0.96 * math.cos(w)]), ("H", [-0.96 * math.sin(w), 0, 0.96 * math.cos(w)])])
nh3 = [("N", [0, 0, 0])]
for k in range(3):
    phi = 2 * math.pi * k / 3
    nh3.append(("H", [0.94 * math.cos(phi), 0.94 * math.sin(phi), -0.38]))
write("ammonia", nh3)
pcl5 = [("P", [0, 0, 0]), ("Cl", [0, 0, 2.12]), ("Cl", [0, 0, -2.12])]
for k in range(3):
    phi = 2 * math.pi * k / 3
    pcl5.append(("Cl", [2.02 * math.cos(phi), 2.02 * math.sin(phi), 0]))
write("pcl5", pcl5)
ethane = [("C", [0, 0, 0.765]), ("C", [0, 0, -0.765])]
for k in range(3):
    phi = 2 * math.pi * k / 3
    ethane.append(("H", [1.02 * math.cos(phi), 1.02 * math.sin(phi), 1.16]))
    ethane.append(("H", [1.02 * math.cos(phi + math.pi / 3), 1.02 * math.sin(phi + math.pi / 3), -1.16]))
write("ethane", ethane)
write("trans_dichloroethylene", [("C", [0.665, 0, 0]), ("C", [-0.665, 0, 0]), ("Cl", [1.55, 1.45, 0]), ("Cl", [-1.55, -1.45, 0]),
                                 ("H", [1.25, -0.92, 0]), ("H", [-1.25, 0.92, 0])])
write("allene", [("C", [0, 0, 0]), ("C", [0, 0, 1.31]), ("C", [0, 0, -1.31]), ("H", [0.93, 0, 1.87]), ("H", [-0.93, 0, 1.87]),
                 ("H", [0, 0.93, -1.87]), ("H", [0, -0.93, -1.87])])
write("co2", [("C", [0, 0, 0]), ("O", [0, 0, 1.16]), ("O", [0, 0, -1.16])])
write("hcn", [("H", [0, 0, -1.066]), ("C", [0, 0, 0]), ("N", [0, 0, 1.156])])
bz = []
for k in range(6):
    phi = 2 * math.pi * k / 6
    noise = [0.02 * math.sin(3.1 * k + 0.4), 0.015 * math.cos(1.7 * k), 0.01 * math.sin(k)]
    bz.append(("C", [1.39 * math.cos(phi) + noise[0], 1.39 * math.sin(phi) + noise[1], noise[2]]))
    bz.append(("H", [2.48 * math.cos(phi) - noise[1], 2.48 * math.sin(phi) + noise[0], -noise[2]]))
write("benzene_distorted", bz)
# Butane (anti) for scans and conformers, with explicit hydrogens.
butane = [("C", [0, 0, 0]), ("C", [1.53, 0, 0]), ("C", [2.04, 1.44, 0]), ("C", [3.57, 1.44, 0])]
for c, (sx, other) in ((0, (-1, 1)), (3, (1, 2))):
    p = butane[c][1]
    butane.append(("H", [p[0] + 0.36 * sx, p[1] - 1.03 if c == 0 else p[1] + 1.03, p[2]]))
    butane.append(("H", [p[0] + 0.36 * sx, p[1] + 0.51 if c == 0 else p[1] - 0.51, p[2] + 0.89]))
    butane.append(("H", [p[0] + 0.36 * sx, p[1] + 0.51 if c == 0 else p[1] - 0.51, p[2] - 0.89]))
for c in (1, 2):
    p = butane[c][1]
    shift = -0.51 if c == 1 else 0.51
    butane.append(("H", [p[0] + (0.36 if c == 1 else -0.36), p[1] + shift, p[2] + 0.89]))
    butane.append(("H", [p[0] + (0.36 if c == 1 else -0.36), p[1] + shift, p[2] - 0.89]))
write("butane", butane)
print("chemistry fixtures written")
