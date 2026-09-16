"""Population analysis and WAVECAR orbital density with explicit conventions."""

import math
from pathlib import Path


def hirshfeld(density, proatoms, *, reference_electrons=None):
    """Stockholder partition using supplied neutral proatom density grids.

    All grids must share geometry, units and periodicity. The caller supplies
    physically appropriate neutral-atom densities; this routine does not invent
    Gaussian substitutes. Returns electron populations and optional net atomic
    charges, reference_electrons - population. A valence density requires valence
    proatoms/reference electron counts; do not mix all-electron and valence data.
    """
    from ..electronic import Grid
    proatoms = list(proatoms)
    if not proatoms or density.unit != "e/A^3":
        raise ValueError("Supply an electron density and at least one proatom")
    rho = density.values
    if any(value < 0 for value in rho):
        raise ValueError("Population analysis requires nonnegative electron density")
    profiles = []
    for reference in proatoms:
        # Native arithmetic verifies complete lattice, origin, shape and units.
        density + reference
        values = reference.values
        if any(value < 0 for value in values):
            raise ValueError("Proatom densities must be nonnegative")
        profiles.append(values)
    denominator = [sum(values) for values in zip(*profiles)]
    if any(total==0 and value!=0 for total,value in zip(denominator,rho)):
        raise ValueError("Proatoms do not cover all nonzero density samples")
    populations = []
    for profile in profiles:
        values = [value * atom / total if total else 0 for value,atom,total in zip(rho,profile,denominator)]
        field = Grid(density.shape,density.cell,values,origin=density.origin,
                     periodic=density.periodic,unit=density.unit)
        populations.append(field.integrate())
    result = {"electron_populations":populations, "total_electrons":density.integrate()}
    if reference_electrons is not None:
        reference_electrons = list(reference_electrons)
        if len(reference_electrons)!=len(populations) or any(not math.isfinite(v) or v<0 for v in reference_electrons):
            raise ValueError("One finite nonnegative reference electron count per atom is required")
        result["net_atomic_charges"] = [reference-population for reference,population in zip(reference_electrons,populations)]
    return result


def _bader_chgcar(source, destination):
    """Canonicalize the total field for legacy fixed-column Henkelman readers.

    Solver input uses six decimal places for lattice and fractional sites, as
    required by the official Fortran reader. Scalar samples retain 16 digits.
    This temporary file never replaces the original density or PAW data.
    """
    import numpy as np
    from pymatgen.io.vasp.outputs import Chgcar
    volume=Chgcar.from_file(source)
    structure=volume.structure
    if not len(structure) or len(volume.poscar.natoms)>110:
        raise ValueError("Bader requires sites and at most 110 contiguous species groups")
    lattice=np.asarray(structure.lattice.matrix)
    if not np.isfinite(lattice).all() or np.abs(lattice).max()>=1e6:
        raise ValueError("Lattice is outside Bader's input range")
    values=np.asarray(volume.data["total"])
    if not np.isfinite(values).all():
        raise ValueError("Density contains non-finite samples")
    with open(destination,"w",encoding="ascii") as handle:
        handle.write("AtomForge Bader input\n1.0000000000000000\n")
        for vector in lattice: handle.write("".join(format(v,"13.6f") for v in vector)+"\n")
        handle.write(" ".join(volume.poscar.site_symbols)+"\n")
        handle.write(" ".join(map(str,volume.poscar.natoms))+"\nDirect\n")
        for point in structure.frac_coords % 1:
            handle.write("".join(format(v,"10.6f") for v in point)+"\n")
        handle.write("\n"+" ".join(map(str,values.shape))+"\n")
        flattened=values.ravel(order="F")
        for start in range(0,len(flattened),5):
            handle.write(" ".join(format(v,".16e") for v in flattened[start:start+5])+"\n")
    return str(destination.resolve())


def bader(path, *, reference=None, potcar=None, executable=None, format="vasp"):
    """Run the Henkelman Bader executable through pymatgen.

    reference can be AECCAR0+AECCAR2 for basin definition while integrating the
    supplied valence density. POTCAR supplies reference valence counts for charge
    transfer; raw electron populations are always returned. No Voronoi fallback.
    Temporary solver inputs use six-decimal lattice/site coordinates for legacy
    fixed-column readers; scalar samples retain 16 digits. Calls use pymatgen's
    working-directory-based runner and should run in separate processes, not threads.
    """
    from pymatgen.command_line.bader_caller import BaderAnalysis
    if format not in ("vasp","cube"):
        raise ValueError("Bader input must be vasp or cube")
    import tempfile
    with tempfile.TemporaryDirectory(prefix="af_bader_") as directory:
        primary=str(Path(path).resolve())
        ref=str(Path(reference).resolve()) if reference else ""
        if format=="vasp":
            if ref:
                import numpy as np
                from pymatgen.io.vasp.outputs import Chgcar
                first,second=Chgcar.from_file(primary),Chgcar.from_file(ref)
                if first.data["total"].shape!=second.data["total"].shape or not np.allclose(
                        first.structure.lattice.matrix,second.structure.lattice.matrix,rtol=1e-8,atol=1e-8):
                    raise ValueError("Bader reference must share the primary grid and lattice")
            primary=_bader_chgcar(primary,Path(directory)/"CHGCAR")
            if ref:
                ref=_bader_chgcar(ref,Path(directory)/"CHGREF")
        analysis = BaderAnalysis(chgcar_filename=primary if format=="vasp" else "",
            cube_filename=primary if format=="cube" else "", chgref_filename=ref,
            potcar_filename=str(Path(potcar).resolve()) if potcar else "", bader_path=executable)
    result = {"atoms":analysis.data, "total_electrons":analysis.nelectrons,
              "vacuum_charge":analysis.vacuum_charge, "vacuum_volume":analysis.vacuum_volume}
    if potcar:
        # Pymatgen's charge_transfer is electrons gained; net atomic charge has the opposite sign.
        result["net_atomic_charges"] = [-analysis.get_charge_transfer(i) for i in range(len(analysis.data))]
    return result


def ddec(directory, *, atomic_densities=None, run=True):
    """Run/read Chargemol DDEC6 and CM5 results through pymatgen.

    Requires the Chargemol executable and its atomic-density database. Input
    directory contains the calculation's CHGCAR/POTCAR/AECCAR files. External
    solver failures propagate; an approximate partition is never substituted.
    """
    from pymatgen.command_line.chargemol_caller import ChargemolAnalysis
    return ChargemolAnalysis(path=str(Path(directory).resolve()),
        atomic_densities_path=str(Path(atomic_densities).resolve()) if atomic_densities else None,
        run_chargemol=run)


def orbital_density(wavecar, poscar, output, *, kpoint=0, band=0, spin=None, spinor=None, scale=2):
    """Export a selected WAVECAR orbital as PARCHG, loadable in the desktop.

    Indices are zero-based. Uses pymatgen's plane-wave reconstruction; PAW
    augmentation is not included, so this is not an all-electron density.
    Occupation and k-point weights are not applied. A collinear spin selection
    returns one component; omitted spin in ISPIN=2 returns total and difference.
    """
    if any(not isinstance(v,int) or v<0 for v in (kpoint,band)) or not isinstance(scale,int) or scale<1:
        raise ValueError("Invalid orbital indices or FFT scale")
    if spin not in (None,0,1) or spinor not in (None,0,1):
        raise ValueError("Spin and spinor must be None, 0 or 1")
    from pymatgen.io.vasp.outputs import Wavecar
    from pymatgen.io.vasp.inputs import Poscar
    import numpy as np
    wave = Wavecar(str(wavecar))
    structure = Poscar.from_file(str(poscar))
    if not np.allclose(wave.a,structure.structure.lattice.matrix,rtol=1e-7,atol=1e-7):
        raise ValueError("POSCAR and WAVECAR lattice vectors differ")
    if kpoint>=wave.nk or band>=wave.nb:
        raise ValueError("Orbital index is outside WAVECAR")
    noncollinear=str(wave.vasp_type).lower().startswith("n")
    if not noncollinear and spinor is not None:
        raise ValueError("Spinor selection requires a noncollinear WAVECAR")
    if (noncollinear and spin is not None) or (wave.spin==1 and spin==1):
        raise ValueError("Spin selection is incompatible with WAVECAR")
    if int(np.prod(wave.ng))*scale**3>100000000:
        raise ValueError("Orbital FFT exceeds 100 million samples")
    # For scalar wavefunctions fft_mesh ignores spinor; explicitly selecting
    # component zero avoids summing the same scalar orbital twice upstream.
    component=spinor if noncollinear else 0
    result = wave.get_parchg(structure,kpoint,band,spin=spin,spinor=component,scale=scale)
    result.write_file(str(output))
    return result
