"""Native grain/IPF sidecars, matched by species and Cartesian position."""

import math
from pathlib import Path
import warnings

_SYMBOLS=("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
          "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
          "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th "
          "Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()


def save_metadata(structure, path):
    """Write a coordinate-matched native sidecar, including optional grain IDs."""
    import os
    import tempfile
    sidecar=Path(path).with_suffix(".atomforge-ipf")
    if not structure.atoms or not any(a.grain_color is not None for a in structure.atoms):
        for stale in {sidecar,Path(str(path)+".atomforge-ipf")}:
            if stale.exists(): stale.unlink()
        return
    if any(a.grain_color is None for a in structure.atoms):
        raise ValueError("Grain colors must be supplied for every atom")
    regions=all(a.grain_region is not None for a in structure.atoms)
    rows=["ATOMFORGE_IPF_V2" if regions else "ATOMFORGE_IPF_V1",str(len(structure))]
    for atom in structure.atoms:
        values=(atom.x,atom.y,atom.z,*atom.grain_color)
        if len(values)!=6 or not all(math.isfinite(v) for v in values):
            raise ValueError("Invalid grain metadata coordinates or colors")
        row=str(_SYMBOLS.index(atom.symbol))+" "+" ".join(format(v,".17g") for v in values)
        if regions:
            if not isinstance(atom.grain_region,int): raise ValueError("Grain IDs must be integers")
            row+=" "+str(atom.grain_region)
        rows.append(row)
    with tempfile.NamedTemporaryFile(mode="w",encoding="utf-8",dir=sidecar.parent,delete=False) as handle:
        temporary=Path(handle.name)
        handle.write("\n".join(rows)+"\n")
    try: os.replace(temporary,sidecar)
    finally: temporary.unlink(missing_ok=True)


def load_metadata(structure, path):
    candidates = [Path(str(path) + ".atomforge-ipf"), Path(path).with_suffix(".atomforge-ipf")]
    sidecar = next((item for item in candidates if item.is_file()), None)
    if sidecar is None:
        return
    try:
        lines = sidecar.read_text(encoding="utf-8").splitlines()
        if lines[0] not in ("ATOMFORGE_IPF_V1", "ATOMFORGE_IPF_V2"):
            raise ValueError("Unknown sidecar version")
        count = int(lines[1])
        if count != len(structure) or len(lines) != count + 2:
            raise ValueError("Sidecar atom count mismatch")
        # Atomic-number order is explicit; never infer it from dictionary order.
        symbols = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
                   "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
                   "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th "
                   "Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()
        pending = []
        remaining = set(range(len(structure)))
        from collections import defaultdict
        from itertools import product
        buckets = defaultdict(list)
        for index, atom in enumerate(structure.atoms):
            buckets[(atom.symbol, *(math.floor(v / .002) for v in (atom.x, atom.y, atom.z)))].append(index)
        for line in lines[2:]:
            values = line.split()
            if len(values) != (8 if lines[0].endswith("V2") else 7):
                raise ValueError("Invalid sidecar row")
            number = int(values[0]); xyz = [float(v) for v in values[1:4]]
            color = tuple(float(v) for v in values[4:7])
            if not 0 < number < len(symbols) or not all(math.isfinite(v) for v in xyz + list(color)):
                raise ValueError("Invalid sidecar values")
            bucket = tuple(math.floor(v / .002) for v in xyz)
            candidates = (i for shift in product((-1, 0, 1), repeat=3)
                          for i in buckets.get((symbols[number], *(bucket[a]+shift[a] for a in range(3))), ()))
            matches = [i for i in candidates if i in remaining
                       and max(abs(a-b) for a,b in zip(xyz, (structure.atoms[i].x,structure.atoms[i].y,structure.atoms[i].z))) < 2e-3]
            if not matches:
                raise ValueError("Sidecar coordinates do not match this structure")
            index = min(matches)
            remaining.remove(index)
            pending.append((index, color, int(values[7]) if len(values)==8 else None))
        for index, color, region in pending:
            structure.atoms[index].grain_color = color
            structure.atoms[index].grain_region = region
    except (ValueError, IndexError, OSError) as error:
        warnings.warn("Ignoring grain metadata: " + str(error), RuntimeWarning, stacklevel=2)
