#!/bin/sh
# Open Babel-written fixtures (MOL, SDF, MOL2, PDB, GJF) for the af-io tests.
# Run from rust/ with the MSYS2 Open Babel on the PATH.
set -e
OUT=tests/fixtures/io
B=/c/msys64/ucrt64
export BABEL_LIBDIR=$(ls -d $B/lib/openbabel/*/ | head -1)
export BABEL_DATADIR=$(ls -d $B/share/openbabel/*/ | head -1)
OB=$B/bin/obabel
$OB -:"CC(=O)Oc1ccccc1C(=O)O aspirin" --gen3d -omol -O $OUT/aspirin.mol 2>/dev/null
$OB -:"Cn1cnc2c1c(=O)n(C)c(=O)n2C caffeine" --gen3d -omol2 -O $OUT/caffeine.mol2 2>/dev/null
$OB -:"CCO ethanol" --gen3d -ogjf -O $OUT/ethanol.gjf 2>/dev/null
$OB -:"c1ccccc1 benzene" --gen3d -opdb -O $OUT/benzene.pdb 2>/dev/null
# Two records in one SDF (the first is read).
$OB -:"O water" --gen3d -osdf -O $OUT/first.sdf 2>/dev/null
$OB -:"N ammonia" --gen3d -osdf -O $OUT/second.sdf 2>/dev/null
cat $OUT/first.sdf $OUT/second.sdf > $OUT/two.sdf
rm $OUT/first.sdf $OUT/second.sdf
# A periodic PDB from a CIF (CRYST1 record).
$OB $OUT/p1.cif -opdb -O $OUT/cell.pdb 2>/dev/null
cat > $OUT/hetatm.pdb <<'EOF'
HEADER    HAND WRITTEN
CRYST1   10.000   11.000   12.000  90.00  95.00  90.00 P 1           1
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.400   2.000   3.000  1.00  0.00           C
HETATM    3 FE   HEM A   2       5.000   5.000   5.000  1.00  0.00          FE
HETATM    4  O   HOH A   3       7.000   7.000   7.000  1.00  0.00
ATOM      5 CL   UNK A   4       8.000   1.000   1.000  1.00  0.00
END
EOF
ls $OUT/*.mol $OUT/*.mol2 $OUT/*.gjf $OUT/*.pdb $OUT/*.sdf
