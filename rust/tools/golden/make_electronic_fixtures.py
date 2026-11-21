"""Fixtures for the af-electronic golden comparison (rust/tests/fixtures/electronic):
volumetric files in every supported format (and malformed ones), LOBSTER
COHPCAR/ICOHPLIST files and a Wannier90 Hamiltonian. Run from the repository
root: python rust/tools/golden/make_electronic_fixtures.py"""
import cmath
import math
import os

OUT = "rust/tests/fixtures/electronic/"
os.makedirs(OUT, exist_ok=True)
BOHR = 0.529177210903


def gaussian_density(points, centres, widths, heights):
    return [sum(h * math.exp(-sum((p[k] - c[k]) ** 2 for k in range(3)) / (2 * w * w)) for c, w, h in zip(centres, widths, heights)) + 1e-4 for p in points]


def write(name, text):
    with open(OUT + name, "w", newline="\n") as f:
        f.write(text)


# ---- Cube: a finite grid (Bohr), x outermost, z innermost.
def cube(name, shape, origin_bohr, steps_bohr, atoms, channels_values, header_channels=None, orbital_ids=None, negative_axis=False):
    lines = ["density test cube", "generated for af-electronic golden tests"]
    count = -len(atoms) if orbital_ids else len(atoms)
    head = "%5d %12.6f %12.6f %12.6f" % (count, *origin_bohr)
    if header_channels is not None:
        head += " %d" % header_channels
    lines.append(head)
    for i in range(3):
        n = -shape[i] if (negative_axis and i == 1) else shape[i]
        lines.append("%5d %12.6f %12.6f %12.6f" % (n, *steps_bohr[i]))
    for z, pos in atoms:
        lines.append("%5d %12.6f %12.6f %12.6f %12.6f" % (z, float(z), *pos))
    if orbital_ids:
        lines.append(" ".join(str(v) for v in [len(orbital_ids)] + orbital_ids))
    data = []
    for x in range(shape[0]):
        for y in range(shape[1]):
            for z in range(shape[2]):
                idx = (z * shape[1] + y) * shape[0] + x
                for ch in channels_values:
                    data.append("%.6E" % ch[idx])
    for i in range(0, len(data), 6):
        lines.append(" ".join(data[i:i + 6]))
    write(name, "\n".join(lines) + "\n")


shape = [12, 10, 8]
steps = [[0.5, 0, 0], [0.05, 0.45, 0], [0, 0.04, 0.55]]
origin = [-0.2, 0.1, 0.0]
points = []
for z in range(shape[2]):
    for y in range(shape[1]):
        for x in range(shape[0]):
            points.append([BOHR * (origin[k] + x * steps[0][k] + y * steps[1][k] + z * steps[2][k]) for k in range(3)])
atoms = [(8, [2.6, 2.4, 1.9]), (1, [4.0, 2.4, 2.0])]
centres = [[BOHR * c for c in a[1]] for a in atoms]
rho = gaussian_density(points, centres, [0.35, 0.3], [2.0, 0.6])
cube("density.cube", shape, origin, steps, atoms, [[v * BOHR ** 3 for v in rho]])
cube("negative_axis.cube", shape, origin, steps, atoms, [[v * BOHR ** 3 for v in rho]], negative_axis=True)
cube("two_channel.cube", shape, origin, steps, atoms, [[v * BOHR ** 3 for v in rho], [math.sin(0.3 * i) for i in range(len(rho))]], header_channels=2)
cube("orbitals.cube", shape, origin, steps, atoms, [[math.cos(0.2 * i) for i in range(len(rho))], [math.sin(0.1 * i) * 0.5 for i in range(len(rho))]], header_channels=1, orbital_ids=[5, 6])
with open(OUT + "density.cube") as f:
    full = f.read()
write("truncated.cube", full[: len(full) // 2])
write("extra.cube", full + "1.0 2.0\n")
write("D_exponent.cube", full.replace("E-", "D-", 3))

# ---- VASP: a periodic grid (density times volume), Fortran column order x fastest.
cell = [[4.0, 0.0, 0.0], [0.3, 3.8, 0.0], [0.0, 0.2, 3.5]]
vol = cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1]) - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0]) + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0])
vshape = [8, 8, 6]
frac_atoms = [("Cu", [0.1, 0.2, 0.3]), ("O", [0.6, 0.55, 0.7])]


def periodic_density(shape, cell, fracs, widths, heights, base=0.01):
    values = []
    for z in range(shape[2]):
        for y in range(shape[1]):
            for x in range(shape[0]):
                f = [x / shape[0], y / shape[1], z / shape[2]]
                total = base
                for (fc, w, h) in zip(fracs, widths, heights):
                    for i in (-1, 0, 1):
                        for j in (-1, 0, 1):
                            for k in (-1, 0, 1):
                                d = [f[0] - fc[0] + i, f[1] - fc[1] + j, f[2] - fc[2] + k]
                                r = [sum(d[a] * cell[a][c] for a in range(3)) for c in range(3)]
                                total += h * math.exp(-sum(c * c for c in r) / (2 * w * w))
                values.append(total)
    return values


prho = periodic_density(vshape, cell, [a[1] for a in frac_atoms], [0.5, 0.4], [3.0, 1.5])
spin = [0.2 * math.sin(i * 0.37) for i in range(len(prho))]


def vasp_text(title, scale_line, symbols, counts, mode_lines, positions, grids, augmentation=False):
    lines = [title, scale_line] + ["  %.10f %.10f %.10f" % tuple(r) for r in cell]
    if symbols:
        lines.append(" ".join(symbols))
    lines.append(" ".join(str(c) for c in counts))
    lines += mode_lines
    lines += ["  %.8f %.8f %.8f" % tuple(p) for p in positions]
    for g in grids:
        lines.append("")
        lines.append(" %d %d %d" % tuple(vshape))
        vals = ["%.11E" % v for v in g]
        for i in range(0, len(vals), 5):
            lines.append(" " + " ".join(vals[i:i + 5]))
        if augmentation:
            lines.append("augmentation occupancies   1  4")
            lines.append("  0.1 0.2 0.3")
            lines.append("  0.4")
            lines.append("augmentation occupancies   2  2")
            lines.append("  0.5 0.6")
    return "\n".join(lines) + "\n"


positions = [a[1] for a in frac_atoms]
write("CHGCAR", vasp_text("Cu O charge", "1.0", ["Cu", "O"], [1, 1], ["Direct"], positions, [[v * vol for v in prho], [v * vol for v in spin]], augmentation=True))
cart = [[sum(p[a] * cell[a][c] for a in range(3)) for c in range(3)] for p in positions]
write("LOCPOT", vasp_text("potential", "-%.10f" % vol, ["Cu", "O"], [1, 1], ["Selective dynamics", "Cartesian"], [[c / 1.0 for c in p] for p in cart], [[-5 + v for v in prho]]))
write("ELFCAR", vasp_text("elf", "1.0 1.0 1.0", [], [1, 1], ["direct"], positions, [[min(1.0, v / 3) for v in prho]]))
write("CHGCAR_four", vasp_text("noncollinear", "1.0", ["Cu", "O"], [1, 1], ["Direct"], positions, [[v * vol for v in prho]] + [[v * vol * 0.1 for v in spin]] * 3))
write("CHGCAR_three", vasp_text("bad", "1.0", ["Cu", "O"], [1, 1], ["Direct"], positions, [[v * vol for v in prho]] * 3))
write("CHGCAR_traj", vasp_text("traj", "1.0", ["Cu", "O"], [1, 1], ["Direct"], positions, [[v * vol for v in prho]]) + "Direct configuration= 2\n")
write("grid.raw", vasp_text("raw grid", "1.0", ["Cu", "O"], [1, 1], ["Direct"], positions, [prho]))

# ---- XSF: two grids in a block, ATOMS and PRIMCOORD sections.
xs = [5, 6, 4]
xcell = [[3.0, 0, 0], [0.4, 2.8, 0], [0, 0, 2.5]]
xvals = [math.exp(-((x - 2) ** 2 + (y - 2.5) ** 2 + (z - 1.5) ** 2) / 3.0) for z in range(xs[2]) for y in range(xs[1]) for x in range(xs[0])]
xsf = ["# comment", "CRYSTAL", "PRIMVEC"] + ["  %.6f %.6f %.6f" % tuple(r) for r in xcell] + ["PRIMCOORD", "2 1", "29 0.1 0.2 0.3", "O 1.5 1.4 1.2", "ATOMS", "6 0.5 0.5 0.5", "H 1.0 1.0 1.0", "BEGIN_BLOCK_DATAGRID_3D", "test"]
for name, scale in (("BEGIN_DATAGRID_3D_rho", 1.0), ("BEGIN_DATAGRID_3D", 2.0)):
    xsf += [name, "%d %d %d" % tuple(xs), "0.1 0.0 -0.1"] + ["%.6f %.6f %.6f" % tuple(r) for r in xcell]
    xsf += [" ".join("%.8f" % (v * scale) for v in xvals[i:i + 6]) for i in range(0, len(xvals), 6)]
    xsf.append("END_DATAGRID_3D")
xsf.append("END_BLOCK_DATAGRID_3D")
write("grid.xsf", "\n".join(xsf) + "\n")
write("nogrid.xsf", "CRYSTAL\nPRIMVEC\n1 0 0\n0 1 0\n0 0 1\n")
write("anim.xsf", "ANIMSTEPS 2\n")

# ---- LOBSTER
energies = [-5.0 + 0.5 * i for i in range(21)]


def cohp_values(b, s, e):
    return -0.3 * (b + 1) * math.exp(-((e - 0.5 * b) ** 2)) * (1 + 0.1 * s), 0.05 * (b + 1) * e * (1 - 0.2 * s)


def cohpcar(name, headers, spins, title="COHPCAR.lobster"):
    lines = [title, "%7d %7d %7d %12.6f %12.6f %12.6f" % (len(headers) + 1, spins, len(energies), energies[0], energies[-1], 5.4321), "Average"] + headers
    for e in energies:
        row = [e]
        for s in range(spins):
            row += [0.1 * e * (1 + s), 0.01 * e]
            for b in range(len(headers)):
                row += list(cohp_values(b, s, e))
        lines.append(" ".join("%.8f" % v for v in row))
    write(name, "\r\n".join(lines) + "\r\n")


cohpcar("COHPCAR.lobster", ["No.1:Fe1->Fe2(2.48000)", "No.2:Fe1[3d_xy]->Fe2[3d_xy](2.48000)", "No.3:Fe1->O13(1.95321)"], 1)
cohpcar("COHPCAR.spin", ["No.1:Ni1->Ni2(2.49)", "No.2:Ni2->O3(2.08)"], 2)
cohpcar("COHPCAR.noheaders", ["No.1:Fe1[3d]->Fe2[3d](2.48)"], 1)
cohpcar("COHPCAR.badheader", ["No.1:Fe1->Fe2 2.48"], 1)
cohpcar("COHPCAR.nolabel", ["No.1:Fe->Fe2(2.48)"], 1)
write("COHPCAR.short", "a\n1 1\n")
write("ICOHPLIST.lobster",
      "  COHP#  atomMU  atomNU  distance  translation  ICOHP(eF)  for spin  1\n"
      "     1     Fe1     Fe2   2.48000   0   0   0   -1.23456\n"
      "     1     Fe1_3d     Fe2_3d   2.48000   0   0   0   -0.5\n"
      "     2     Fe1     O13   1.95321   1   0  -1   -2.5\n")
write("ICOHPLIST.spin",
      "  COHP#  atomMU  atomNU  distance  translation  ICOHP(eF)  for spin  1\n"
      "     1     Ni1     Ni2   2.49   0   0   0   -1.1\n"
      "     2     Ni2     O3   2.08   0   0   0   -2.2\n"
      "  COHP#  atomMU  atomNU  distance  translation  ICOHP(eF)  for spin  2\n"
      "     1     Ni1     Ni2   2.49   0   0   0   -1.0\n"
      "     2     Ni2     O3   2.08   0   0   0   -2.0\n\n")
write("ICOHPLIST.old",
      "  COHP#  atomMU  atomNU  distance  ICOHP(eF)  #bonds\n"
      "     1     Na1     Cl2   2.82   -0.75   6\n"
      "     2     Na1     Na3   3.99   -0.05   12\n")
write("ICOHPLIST.bad", "header\n  1 Fe1 Fe2 2.48 -1.2\n")

# ---- Wannier90: a two-band Haldane model on a triangular Bravais lattice.
t1, t2, phi, m = -1.0, 0.15, math.pi / 2, 0.2
hop = {}


def add_hop(R, a, b, v):
    hop.setdefault(R, [[0j, 0j], [0j, 0j]])
    hop[R][a][b] += v
    mr = (-R[0], -R[1], -R[2])
    hop.setdefault(mr, [[0j, 0j], [0j, 0j]])
    if mr != R or a != b:
        hop[mr][b][a] += v.conjugate()


hop[(0, 0, 0)] = [[m + 0j, 0j], [0j, -m + 0j]]
for R in [(0, 0, 0), (-1, 0, 0), (0, -1, 0)]:
    if R == (0, 0, 0):
        hop[R][0][1] += t1
        hop[R][1][0] += t1
    else:
        add_hop(R, 0, 1, complex(t1))
for R, sign in [((1, 0, 0), 1), ((0, 1, 0), -1), ((-1, -1, 0), 1)]:
    add_hop(R, 0, 0, t2 * cmath.exp(1j * phi * sign))
    add_hop(R, 1, 1, t2 * cmath.exp(-1j * phi * sign))
vectors = sorted(hop)
deg = [1 if R != (0, 0, 0) else 1 for R in vectors]
lines = ["written for af-electronic golden tests", "2", str(len(vectors))]
for i in range(0, len(deg), 15):
    lines.append(" ".join("%5d" % d for d in deg[i:i + 15]))
for R in vectors:
    for n in range(2):
        for mm in range(2):
            v = hop[R][mm][n]
            lines.append("%5d%5d%5d%5d%5d%12.6f%12.6f" % (R[0], R[1], R[2], mm + 1, n + 1, v.real, v.imag))
write("haldane_hr.dat", "\n".join(lines) + "\n")
nonherm = list(lines)
nonherm[-1] = nonherm[-1][:-12] + "%12.6f" % 0.9
write("nonhermitian_hr.dat", "\n".join(nonherm) + "\n")
write("short_hr.dat", "x\n2\n")
print("electronic fixtures written")
