"""COHP/COOP/COBI bonding curves and their Fermi-level integral from LOBSTER.

Thin ctypes wrapper over the native engine's LOBSTER readers (see
src/electronic/Lobster.h): Local Orbital Basis Suite Towards
Electronic-Structure Reconstruction (Nelson et al., J. Comput. Chem. 2020),
implementing the Crystal Orbital Hamilton Population method (Dronskowski and
Bloechl, J. Phys. Chem. 1993). COHPCAR/COOPCAR/COBICAR hold the full pCOHP
(etc.) vs energy curve and its running integral per bond; ICOHPLIST/
ICOOPLIST/ICOBILIST hold only that integral evaluated at the Fermi level,
one row per bond -- the compact per-bond bond-strength summary most analyses
actually want. Orbital-resolved and multi-center-COBI files are not
supported and raise a clear error rather than being silently mis-parsed.
"""

import ctypes as ct

from ._native import library


def _pull(pointer, lib, free, extract):
    if not pointer:
        raise ValueError(lib.af_lobster_error().decode("utf-8", errors="replace"))
    try:
        return extract(pointer)
    finally:
        free(pointer)


def read_cohpcar(path):
    """Return a dict with the full COHP/COOP/COBI curve data from a COHPCAR-style file.

    Keys: spin_polarized, fermi_energy_eV, energies_eV (Fermi level at 0 eV),
    bonds (list of (atom1, atom2, length_A) with 0-based site indices),
    average (one {"cohp":[...],"icohp":[...]} dict per spin channel) and
    bond_curves (one list of per-bond {"cohp":[...],"icohp":[...]} dicts,
    in bond order, per spin channel).
    """
    lib = library()
    pointer = lib.af_cohpcar_load(str(path).encode("utf-8"))

    def extract(pointer):
        spin_polarized, num_energies, num_bonds = ct.c_int(), ct.c_size_t(), ct.c_int()
        if not lib.af_cohpcar_info(pointer, ct.byref(spin_polarized), ct.byref(num_energies), ct.byref(num_bonds)):
            raise ValueError(lib.af_lobster_error().decode("utf-8", errors="replace"))
        n, spins = num_energies.value, 2 if spin_polarized.value else 1
        energies = (ct.c_double * n)()
        lib.af_cohpcar_energies(pointer, energies)
        bonds = []
        for b in range(num_bonds.value):
            atom1, atom2, length = ct.c_int(), ct.c_int(), ct.c_double()
            lib.af_cohpcar_bond(pointer, b, ct.byref(atom1), ct.byref(atom2), ct.byref(length))
            bonds.append((atom1.value, atom2.value, length.value))
        average, curves = [], []
        for spin in range(spins):
            cohp, icohp = (ct.c_double * n)(), (ct.c_double * n)()
            lib.af_cohpcar_average(pointer, spin, cohp, icohp)
            average.append({"cohp": list(cohp), "icohp": list(icohp)})
            per_bond = []
            for b in range(num_bonds.value):
                bond_cohp, bond_icohp = (ct.c_double * n)(), (ct.c_double * n)()
                lib.af_cohpcar_bond_curve(pointer, b, spin, bond_cohp, bond_icohp)
                per_bond.append({"cohp": list(bond_cohp), "icohp": list(bond_icohp)})
            curves.append(per_bond)
        return {"spin_polarized": bool(spin_polarized.value), "fermi_energy_eV": lib.af_cohpcar_fermi_energy(pointer),
                "energies_eV": list(energies), "bonds": bonds, "average": average, "bond_curves": curves}

    return _pull(pointer, lib, lib.af_cohpcar_free, extract)


def read_icohplist(path):
    """Return a dict with the per-bond Fermi-level ICOHP/ICOOP/ICOBI summary.

    Keys: spin_polarized, entries (list of dicts with atom1, atom2 [0-based
    site indices], length_A, num_bonds, icohp_up, icohp_down [icohp_down is
    0.0 when not spin_polarized]).
    """
    lib = library()
    pointer = lib.af_icohplist_load(str(path).encode("utf-8"))

    def extract(pointer):
        entries = []
        for i in range(lib.af_icohplist_count(pointer)):
            atom1, atom2, length, num_bonds = ct.c_int(), ct.c_int(), ct.c_double(), ct.c_int()
            icohp_up, icohp_down = ct.c_double(), ct.c_double()
            if not lib.af_icohplist_entry(pointer, i, ct.byref(atom1), ct.byref(atom2), ct.byref(length),
                                           ct.byref(num_bonds), ct.byref(icohp_up), ct.byref(icohp_down)):
                raise ValueError(lib.af_lobster_error().decode("utf-8", errors="replace"))
            entries.append({"atom1": atom1.value, "atom2": atom2.value, "length_A": length.value,
                             "num_bonds": num_bonds.value, "icohp_up": icohp_up.value, "icohp_down": icohp_down.value})
        return {"spin_polarized": bool(lib.af_icohplist_spin_polarized(pointer)), "entries": entries}

    return _pull(pointer, lib, lib.af_icohplist_free, extract)
