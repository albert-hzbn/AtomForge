"""Trajectory fixtures for the af-io tests (tests/fixtures/traj). Run from rust/."""
import os

out = 'tests/fixtures/traj'
os.makedirs(out, exist_ok=True)


def put(name, text):
    with open(os.path.join(out, name), 'w', newline='\n') as f:
        f.write(text)


put('frames.xyz', '''2
frame 1
Cu 0 0 0
cu 1.8 0 0

2
Lattice="4 0 0 0 4 0 0 0 4" Properties=species:S:1:pos:R:3
CU 0.1 0 0
Cu 1.9 0 0
''')
put('props.extxyz', '''3
Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=id:I:1:species:S:1:pos:R:3:forces:R:3 pbc="T T T"
1 Ni 0.0 0.0 0.0 0.1 0.2 0.3
2 Al 2.5 2.5 0.0 0.0 0.0 0.0
3 Ni 2.5 0.0 2.5 0.0 0.0 0.0
''')
put('element.dump', '''ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0.0 4.0
0.0 4.0
0.0 4.0
ITEM: ATOMS id element x y z
3 Cu 2.0 2.0 0.0
1 cu 0.0 0.0 0.0
2 Ni 2.0 0.0 2.0
ITEM: TIMESTEP
100
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0.0 4.1
0.0 4.1
0.0 4.1
ITEM: ATOMS id element xu yu zu
1 Cu 0.05 0.0 0.0
2 Ni 2.05 0.0 2.0
3 Cu 2.05 2.0 0.0
''')
put('scaled_triclinic.lammpstrj', '''ITEM: TIMESTEP
5
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz pp pp pp
-0.5 5.5 1.0
0.0 4.33 0.5
0.0 6.0 0.25
ITEM: ATOMS id type xs ys zs
2 2 0.5 0.5 0.5
1 1 0.0 0.0 0.0
''')
put('type_only.dump', '''ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 3
0 3
0 3
ITEM: ATOMS id type x y z
1 1 0 0 0
2 2 1.5 1.5 1.5
''')
put('bad_count.xyz', '''abc
x
''')
put('truncated.dump', '''ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 3
0 3
0 3
ITEM: ATOMS id element x y z
1 Cu 0 0 0
''')
print('trajectory fixtures in', out)
