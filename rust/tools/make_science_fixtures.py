"""Writes the science-data fixtures (rust/tests/fixtures/science) read by both
the C++ golden dumper and the Rust tests."""
import os
import struct
import sys

out = sys.argv[1] if len(sys.argv) > 1 else "rust/tests/fixtures/science"
os.makedirs(out, exist_ok=True)


def put(name, text, binary=False):
    mode = "wb" if binary else "w"
    with open(os.path.join(out, name), mode, **({} if binary else {"newline": "\n"})) as f:
        f.write(text)


put("momenta.extxyz", """2
Lattice="4 0 0 0 4 0 0 0 4" Properties=species:S:1:pos:R:3:momenta:R:3 time_fs=0 pbc="T T F"
Cu 0 0 0 0.1 0.2 0.3
Au 1 1.5 2 -0.4 0 0.25
2
Lattice="4 0 0 0 4 0 0 0 4" Properties=species:S:1:pos:R:3:momenta:R:3 time_fs=2.5 pbc="T T F"
Cu 0.1 0 0 0.1 0.2 0.3
Au 1 1.6 2 -0.4 0 0.25
""")
put("velocities.extxyz", """3
Properties=species:S:1:pos:R:3:vel:R:3:masses:R:1:property:R:1 comment="a b=c" flag
O 0 0 0 1 2 3 15.5 0.25
H 0.96 0 0 0 0 0.5 1.1 -1
H -0.24 0.93 0 0 0.1 0 1.0 7
""")
put("znumbers.extxyz", """2
Properties=Z:I:1:pos:R:3
29 0 0 0
79 1 1 1
""")
put("crlf.xyz", "2\r\nplain\r\nC 0 0 0\r\nO 1.128 0 0\r\n\r\n")
put("frames.xyz", "1\nf0\nNa 0 0 0\n1\nf1\nNa 0.5 0 0\n1\nf2\nNa 1 0 0\n")
put("bad_props.extxyz", "1\nProperties=species:S:1:pos:R\nCu 0 0 0\n")
put("noncubic_POSCAR", """Si
 -40.0
 2.7 2.7 0
 0 2.7 2.7
 2.7 0 2.7
 Si_pv C/abc
 1 1
Selective dynamics
Cartesian
 0 0 0 T T T
 1.35 1.35 1.35 F F F
""")
put("XDATCAR", """Ar
1.0
5 0 0
0 5 0
0 0 5
Ar
2
Direct configuration=     1
 0 0 0
 0.5 0.5 0.5
Direct configuration=     2
 0.1 0 0
 0.5 0.6 0.5
""")
put("XDATCAR_variable", """Ar
1.0
5 0 0
0 5 0
0 0 5
Ar
2
Direct configuration=     1
 0 0 0
 0.5 0.5 0.5
Ar
1.0
5.5 0 0
0 5 0
0 0 5
Ar
2
Cartesian configuration=     2
 0.1 0 0
 2.5 3.0 2.5
""")
put("structure.json", """{"symbols": ["Fe", "Fe"], "positions": [[0,0,0],[1.435,1.435,1.435]],
 "cell": [[2.87,0,0],[0,2.87,0],[0,0,2.87]], "pbc": [true, true, false],
 "magmoms": [2.2, -2.2], "fixed": [[true,false,true],[false,false,false]],
 "velocities": [[0.01,0,0],[0,0,-0.01]], "charges": [0.1, -0.1]}""")
put("molecule.json", """{"symbols": ["N", "N"], "positions": [[0,0,0],[0,0,1.1]]}""")
put("bad_key.json", """{"symbols": ["N"], "positions": [[0,0,0]], "colour": 1}""")
put("bad_cell.json", """{"symbols": ["N"], "positions": [[0,0,0]], "cell": [[1,0,0],[2,0,0],[0,0,1]]}""")
put("table.csv", "# x, y\n1, 2.5\n3,4e-1 # tail\n\n5,6\n")
put("column.dat", "1\n2\n3.5\n")
put("ragged.txt", "1 2\n3\n")
put("EIGENVAL", """    4    4    1    2
  0.1 0.1 0.1 0.1
  1e-16
  CAR
  unknown system
    8    2    3

  0.0 0.0 0.0 0.5
    1   -5.0  -5.1
    2    1.0   1.1
    3    4.0   4.2

  0.5 0.0 0.0 0.5
    1   -4.0  -4.1
    2    2.0   2.1
    3    5.0   5.2
""")


def npy(name, descr, shape, data, fortran=False, version=1):
    header = "{'descr': '%s', 'fortran_order': %s, 'shape': (%s), }" % (
        descr, "True" if fortran else "False", "".join(f"{s}, " for s in shape) if len(shape) != 1 else f"{shape[0]},")
    pad = 64 - (len(header) + (10 if version == 1 else 12) + 1) % 64
    header = header + " " * pad + "\n"
    head = b"\x93NUMPY" + bytes([version, 0])
    head += struct.pack("<H", len(header)) if version == 1 else struct.pack("<I", len(header))
    put(name, head + header.encode() + data, binary=True)


npy("f8.npy", "<f8", (2, 3), struct.pack("<6d", 1, 2, 3, 4, 5.5, -6))
npy("f4_fortran.npy", "<f4", (2, 3), struct.pack("<6f", 1, 4, 2, 5, 3, 6), fortran=True)
npy("i4.npy", "<i4", (4,), struct.pack("<4i", 1, -2, 3, 40000))
npy("i2_v2.npy", "<i2", (2, 2), struct.pack("<4h", 1, -2, 3, 4), version=2)
npy("u1.npy", "|u1", (3,), bytes([0, 128, 255]))
npy("b1.npy", "|b1", (3,), bytes([0, 1, 2]))
npy("big.npy", ">f8", (1,), struct.pack(">d", 1))
npy("short.npy", "<f8", (5,), struct.pack("<2d", 1, 2))
print("wrote", out)

import bz2, gzip, lzma
raw = open(os.path.join(out, "EIGENVAL"), "rb").read()
open(os.path.join(out, "EIGENVAL.gz"), "wb").write(gzip.compress(raw, mtime=0))
open(os.path.join(out, "EIGENVAL.bz2"), "wb").write(bz2.compress(raw))
open(os.path.join(out, "EIGENVAL.xz"), "wb").write(lzma.compress(raw))
open(os.path.join(out, "EIGENVAL_bad.gz"), "wb").write(b"not gzip at all")
print("compressed")
