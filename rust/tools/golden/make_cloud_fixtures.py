"""Inputs for the af-cloud golden comparison (rust/tests/fixtures/cloud): XYZ
and LAMMPS dump files, valid and malformed. Run from the repository root:
python rust/tools/golden/make_cloud_fixtures.py"""
import os
import random

OUT = "rust/tests/fixtures/cloud/"
os.makedirs(OUT, exist_ok=True)
random.seed(4242)


def write(name, text):
    with open(OUT + name, "w", newline="\n") as f:
        f.write(text)


def atoms(n, box, symbols):
    return [(random.choice(symbols), [random.uniform(0, box[k]) for k in range(3)]) for _ in range(n)]


# Extended XYZ with a lattice and an extra column before the species.
box = [30.0, 25.0, 20.0]
frame = atoms(3000, box, ["Cu", "Ni", "Al"])
lines = [str(len(frame)), 'Lattice="30.0 0.0 0.0 0.0 25.0 0.0 0.0 0.0 20.0" Properties=id:I:1:species:S:1:pos:R:3 pbc="T T T"']
for i, (s, p) in enumerate(frame):
    lines.append("%d %s %.6f %.6f %.6f" % (i + 1, s, *p))
write("big.extxyz", "\n".join(lines) + "\n")
# A two-frame plain XYZ (CRLF line ends).
frames = []
for f in range(2):
    fr = atoms(1500, [12, 14, 9], ["O", "H"])
    frames.append(str(len(fr)))
    frames.append("frame %d" % f)
    frames += ["%s %.5f %.5f %.5f" % (s, *p) for s, p in fr]
write("two_frames.xyz", "\r\n".join(frames) + "\r\n")
write("bad_coordinate.xyz", "2\ncomment\nC 0 0 0\nC 1 x 0\n")
write("short.xyz", "5\ncomment\nC 0 0 0\n")


def dump(timestep, items, columns, rows, box_lines, triclinic=False):
    head = ["ITEM: TIMESTEP", str(timestep), "ITEM: NUMBER OF ATOMS", str(len(rows))]
    head.append("ITEM: BOX BOUNDS xy xz yz pp pp pp" if triclinic else "ITEM: BOX BOUNDS pp pp pp")
    head += box_lines
    head.append("ITEM: ATOMS " + " ".join(columns))
    return "\n".join(head + rows) + "\n"


rows = []
for i in range(2500):
    e = random.choice(["Cu", "Zr"])
    p = [random.uniform(-5, 35), random.uniform(0, 20), random.uniform(10, 30)]
    rows.append("%d %s %.6f %.6f %.6f %.3f" % (i + 1, e, *p, random.uniform(-1, 1)))
write("elements.dump", dump(0, None, ["id", "element", "x", "y", "z", "c_pe"], rows, ["-5 35", "0 20", "10 30"]))
rows2 = []
for i in range(2200):
    t = random.choice([1, 2, 3])
    rows2.append("%d %d %.6f %.6f %.6f" % (i + 1, t, random.uniform(0, 25), random.uniform(0, 18), random.uniform(0, 15)))
tri = dump(100, None, ["id", "type", "xu", "yu", "zu"], rows2, ["-2.0 27.0 1.5", "0.0 18.0 -1.0", "0.0 15.0 0.5"], triclinic=True)
write("triclinic.dump", tri)
rows3 = []
for i in range(1800):
    rows3.append("%d %d %.6f %.6f %.6f" % (i + 1, random.choice([1, 2]), random.random(), random.random(), random.random()))
first = dump(0, None, ["id", "type", "xs", "ys", "zs"], rows3[:900], ["0 20", "0 21", "0 22"])
second = dump(50, None, ["id", "type", "xs", "ys", "zs"], rows3[900:], ["0 20.5", "0 21.5", "0 22.5"])
write("scaled_two_frames.dump", first + second)
write("no_positions.dump", dump(0, None, ["id", "type", "vx", "vy", "vz"], ["1 1 0 0 0"], ["0 1", "0 1", "0 1"]))
write("short_line.dump", dump(0, None, ["id", "type", "x", "y", "z"], ["1 1 0 0 0", "2 1 0.5"], ["0 1", "0 1", "0 1"]))
write("not_a_dump.dump", "ITEM: TIMESTEP\n0\nITEM: BOX\n")
write("text.txt", "hello\n")
print("cloud fixtures written")
