"""Writes the structure-file fixtures for the af-io tests into
tests/fixtures/io. Run from rust/. Open Babel-converted fixtures (MOL, SDF,
MOL2, PDB, GJF) are produced separately by tools/make_io_fixtures.sh."""
import os

out = 'tests/fixtures/io'
os.makedirs(out, exist_ok=True)


def put(name, text):
    with open(os.path.join(out, name), 'w', newline='\n', encoding='utf-8') as f:
        f.write(text)


# ---------------------------------------------------------------- XYZ
put('labels.xyz', '''5
plain xyz with labels
Cu1 0 0 0
CA\t1.5 0.0 0.0
29 0 1.8 0
o 0 0 1.2
Fe_up 2.0 2.0 2.0
''')
put('blank_first.xyz', '''

  3
comment
H 0 0 0
H 0 0 0.74
He 5 5 5
''')
put('cell.extxyz', '''4
Lattice="4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0" Properties=species:S:1:pos:R:3:initial_magmoms:R:1:initial_charges:R:1:velo:R:3 pbc="T T F"
Fe 0.0 0.0 0.0 2.2 0.1 0.01 0.0 0.0
Fe 2.0 2.0 2.0 -2.2 -0.1 0.0 0.02 0.0
O 4.5 0.0 0.0 0.0 0.0 0.0 0.0 0.03
O -0.5 3.9999999 1.0 0.0 0.0 0.0 0.0 0.0
''')
put('molecule_box.extxyz', '''3
Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3 pbc="F F F"
O 0 0 0
H 0.76 0.59 0
H -0.76 0.59 0
''')
put('z_column.extxyz', '''2
Lattice="3 0 0 0 3 0 0 0 3" Properties=id:I:1:Z:I:1:pos:R:3
1 29 0 0 0
2 28 1.5 1.5 0
''')
put('bad_count.xyz', '''5
too few
H 0 0 0
''')
put('bad_element.xyz', '''1
x
Qq 0 0 0
''')
put('triclinic_wrap.extxyz', '''3
Lattice="3.0 0.0 0.0 1.5 2.6 0.0 0.5 0.5 3.0" Properties=species:S:1:pos:R:3
Si -0.3 0.1 0.2
Si 4.6 2.7 3.1
Si 1.0 1.0 1.0
''')

# ---------------------------------------------------------------- POSCAR
put('nacl.vasp', '''Na Cl rocksalt
5.64
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Na Cl
4 4
Direct
0.0 0.0 0.0
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5
0.5 0.0 0.0
0.0 0.5 0.0
0.0 0.0 0.5
0.5 0.5 0.5
''')
put('vasp4.vasp', '''Na Cl
1.0
5.64 0 0
0 5.64 0
0 0 5.64
1 1
Cartesian
0 0 0
2.82 2.82 2.82
''')
put('selective.vasp', '''slab
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 12.0
Fe_pv O
2 1
Selective dynamics
Direct
0.0 0.0 0.0 F F F
0.5 0.5 0.1 T T T
0.0 0.0 0.2 T F T
''')
put('negative_scale.vasp', '''volume given
-64.0
2.0 0.0 0.0
0.0 2.0 0.0
0.0 0.0 2.0
Cu
1
Cartesian
1.0 1.0 1.0
''')
put('per_axis.vasp', '''scales
2.0 3.0 4.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Al
2
Direct
0 0 0
0.5 0.5 1.25
''')
put('POSCAR-001', '''phonopy displacement
1.0
3.6 0 0
0 3.6 0
0 0 3.6
Cu
4
Direct
0.01 0 0
0.5 0.5 0
0.5 0 0.5
0 0.5 0.5
''')
put('bad_poscar.vasp', '''title
1.0
1 0 0
0 1 0
''')

# ---------------------------------------------------------------- CIF
put('p1.cif', '''data_p1
_cell_length_a 4.0
_cell_length_b 5.0
_cell_length_c 6.0
_cell_angle_alpha 90
_cell_angle_beta 100.5(3)
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1 Fe3+ 0.1 0.2 0.3
O1 O2- 0.6 0.7 1.05
''')
put('fm3m_symops.cif', '''data_cu
_cell_length_a 3.615(2)
_cell_length_b 3.615
_cell_length_c 3.615
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F m -3 m'
loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'
'x, y+1/2, z+1/2'
'x+1/2, y, z+1/2'
'x+1/2, y+1/2, z'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cu 0.0 0.0 0.0
''')
put('fd3m_name.cif', '''data_si
_cell_length_a 5.431
_cell_length_b 5.431
_cell_length_c 5.431
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F d -3 m'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 0.125 0.125 0.125
''')
put('fd3m_origin1.cif', '''data_si
_cell_length_a 5.431
_cell_length_b 5.431
_cell_length_c 5.431
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F d -3 m'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 0 0 0
''')
put('hcp_number.cif', '''data_mg
_cell_length_a 3.21
_cell_length_b 3.21
_cell_length_c 5.21
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_space_group_IT_number 194
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Mg1 0.3333 0.6667 0.25
''')
put('r3m_hex.cif', '''data_bi
_cell_length_a 4.546
_cell_length_b 4.546
_cell_length_c 11.862
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_symmetry_space_group_name_H-M 'R -3 m H'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Bi1 0 0 0.2339
''')
put('hall.cif', '''data_rutile
_cell_length_a 4.594
_cell_length_b 4.594
_cell_length_c 2.959
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_Hall '-P 4n 2n'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Ti1 0 0 0
O1 0.3048 0.3048 0
''')
put('disorder.cif', '''data_alloy
_cell_length_a 3.6
_cell_length_b 3.6
_cell_length_c 3.6
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Cu1 Cu 0 0 0 0.4
Ni1 Ni 0 0 0 0.6
Au1 Au 0.5 0.5 0 0.5
Au2 Au 0.52 0.5 0 0.5
O1 O 0.5 0 0.5 0.7
X1 C 0 0.5 0.5 0
H1 H ? ? ? 1
''')
put('blocks.cif', '''# two structures and a text field
data_first
_publ_section_comment
;
A text field with _fake_tag and loop_ inside.
;
_cell_length_a 3
_cell_length_b 3
_cell_length_c 3
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_calc_flag
Na1 0 0 0 .
Cg1 0.5 0.5 0.5 dum
data_second
_cell_length_a 4
_cell_length_b 4
_cell_length_c 4
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
K1 0 0 0
''')
put('magnetic.cif', '''data_mag
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_space_group_symop_magn_operation.xyz
x,y,z,+1
-x,-y,-z,+1
loop_
_space_group_symop_magn_centering.xyz
x,y,z,+1
x+1/2,y+1/2,z+1/2,-1
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Mn1 0.1 0.1 0.1
''')
put('special_rounded.cif', '''data_rounded
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_space_group_symop_operation_xyz
x,y,z
-x,-y,z
-x,y,-z
x,-y,-z
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Zn1 0.0001 0.25 0.25
''')
put('cif2_triple.cif', '''#\\#CIF_2.0
data_c2
_cell.length_a 2.5
_cell.length_b 2.5
_cell.length_c 4.0
_cell.angle_gamma 120
_note """a triple-quoted
value over two lines"""
loop_
_atom_site.label
_atom_site.fract_x
_atom_site.fract_y
_atom_site.fract_z
C1 0 0 0
C2 0.3333333 0.6666667 0
''')
put('bad_cell.cif', '''data_bad
_cell_length_a 3
_cell_length_c 3
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1 0 0 0
''')
put('unknown_element.cif', '''data_bad
_cell_length_a 3
_cell_length_b 3
_cell_length_c 3
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
1Xx 0 0 0
''')

# ---------------------------------------------------------------- LAMMPS data
put('atomic.data', '''LAMMPS data file written by someone

4 atoms
2 atom types

0.0 4.0 xlo xhi
0.0 4.0 ylo yhi
-1.0 3.0 zlo zhi

Masses

1 63.546  # Cu
2 58.6934

Atoms  # atomic

3 2 2.0 2.0 0.0
1 1 0.0 0.0 0.0 0 0 1
2 1 2.0 0.0 2.0
4 2 0.0 2.0 2.0
''')
put('charge_triclinic.data', '''tilted box

3 atoms
2 atom types
0 5 xlo xhi
0 4.33 ylo yhi
0 6 zlo zhi
2.5 0.5 0.25 xy xz yz

Masses

1 22.99
2 35.45

Atoms # charge

1 1 1.0 0.0 0.0 0.0
2 2 -1.0 2.5 2.165 3.0
3 1 0.5 1.0 1.0 1.0

Velocities

1 1.0 2.0 3.0
3 -1.0 0.0 0.0
''')
put('full_labels.data', '''full style with type labels

2 atoms
2 atom types
0 10 xlo xhi
0 10 ylo yhi
0 10 zlo zhi

Atom Type Labels

1 O
2 H

Atoms # full

1 1 1 -0.8 1.0 1.0 1.0 0 0 0
2 1 2 0.4 1.96 1.0 1.0 0 0 0
''')
put('detect_charge.data', '''detect charge style from columns

2 atoms
1 atom types
0 3 xlo xhi
0 3 ylo yhi
0 3 zlo zhi

Masses

1 26.98

Atoms

1 1 0.5 0 0 0
2 1 -0.5 1.5 1.5 1.5
''')
put('ambiguous.data', '''molecular or charge?

1 atoms
1 atom types
0 3 xlo xhi
0 3 ylo yhi
0 3 zlo zhi

Masses

1 12.011

Atoms

1 1 1 0 0 0
''')
put('forcefield_types.data', '''CA is calcium only if its mass is right

2 atoms
2 atom types
0 10 xlo xhi
0 10 ylo yhi
0 10 zlo zhi

Masses

1 12.011 # CA
2 1.008 # HO

Atoms # atomic

1 1 1 1 1
2 2 2 2 2
''')

# ---------------------------------------------------------------- simulation inputs
put('crystal.xsf', '''# comment
CRYSTAL
PRIMVEC
 3.0 0.0 0.0
 0.0 3.0 0.0
 0.0 0.0 3.0
PRIMCOORD
 2 1
26 0.0 0.0 0.0
Fe 1.5 1.5 1.5
''')
put('slab.xsf', '''SLAB
PRIMVEC
 3.0 0.0 0.0
 0.0 3.0 0.0
 0.0 0.0 20.0
PRIMCOORD
 1 1
Pt 0.0 0.0 5.0
''')
put('molecule.xsf', '''ATOMS
 8 0.0 0.0 0.0
 1 0.76 0.59 0.0
 1 -0.76 0.59 0.0
''')
put('fractional.gen', '''2 F
Ga As
1 1 0.0 0.0 0.0
2 2 0.25 0.25 0.25
0.0 0.0 0.0
0.0 2.83 2.83
2.83 0.0 2.83
2.83 2.83 0.0
''')
put('cluster.gen', '''# a cluster
3 C
O H
1 1 0.0 0.0 0.0
2 2 0.76 0.59 0.0
3 2 -0.76 0.59 0.0
''')
put('periodic.gen', '''1 S
Cu
1 1 0.5 0.5 0.5
1.0 0.0 0.0
2.0 0.0 0.0
0.0 2.0 0.0
0.0 0.0 2.0
''')
put('abs.cell', '''%BLOCK LATTICE_CART
ANG
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
%ENDBLOCK LATTICE_CART

%block positions_abs
Fe 0.0 0.0 0.0 SPIN=2.5
Fe 2.0 2.0 2.0 spin = -2.5
O 2.0 0.0 0.0 SPIN : 0.1
%endblock positions_abs

%BLOCK IONIC_CONSTRAINTS
1 Fe 1 1 0 0
2 Fe 2 0 0 1
3 O 1 0.5 0.5 0
%ENDBLOCK IONIC_CONSTRAINTS
''')
put('abc.cell', '''%BLOCK LATTICE_ABC
bohr
6.0 6.0 8.0
90 90 120
%ENDBLOCK LATTICE_ABC
%BLOCK POSITIONS_FRAC
! comment
Zn 0.333333 0.666667 0.0
O 0.333333 0.666667 0.38
%ENDBLOCK POSITIONS_FRAC
''')
put('ibrav0.pwi', '''&CONTROL
  calculation = 'relax'
/
&SYSTEM
  ibrav = 0, celldm(1) = 10.2, nat = 2, ntyp = 1
/
&ELECTRONS
/
ATOMIC_SPECIES
  Si 28.086 Si.UPF
CELL_PARAMETERS alat
  -0.5 0.0 0.5
   0.0 0.5 0.5
  -0.5 0.5 0.0
ATOMIC_POSITIONS crystal
  Si 0.00 0.00 0.00 0 0 0
  Si 0.25 0.25 0.25
K_POINTS automatic
  4 4 4 0 0 0
''')
put('angstrom.pwi', '''&SYSTEM
  ibrav=0
  nat=2
/
CELL_PARAMETERS {angstrom}
  3.0 0.0 0.0
  0.0 3.0 0.0
  0.0 0.0 3.0
ATOMIC_POSITIONS (bohr)
  Cu 0.0 0.0 0.0
  Cu 2.8345 2.8345 2.8345
''')
put('a_param.pwi', '''&SYSTEM
  ibrav = 0, A = 4.0, nat = 1
/
CELL_PARAMETERS
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0d0
ATOMIC_POSITIONS angstrom
  Ag 1.0 1.0 1.0
''')

# ---------------------------------------------------------------- misc
put('unsupported.foo', 'nothing\n')
put('empty.xyz', '')
put('missing_ext', 'x\n')
print('fixtures written to', out)
