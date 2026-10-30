"""Builder and editing steps checked against ASE and geometric invariants.

Each new structure tool is compared with an independent construction
(ase.build, ase.cluster, ASE's Niggli reduction and adsorption sites) or with
exact properties: bond counts and lengths, symmetry, conservation of atoms,
distances preserved by rigid motions, densities.
"""
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest

if sys.platform == "win32" and os.path.isdir(r"C:\msys64\ucrt64\bin"):
    os.add_dll_directory(r"C:\msys64\ucrt64\bin")  # runtime of a development build
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, os.path.dirname(__file__))
from _optional import missing, requires
import atomforge as af
from atomforge.pipeline import Pipeline, PipelineError

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None


def run(text, structure=None, **kwargs):
    return Pipeline.parse(text).run(structure or af.Structure(), **kwargs)


def positions(s):
    return np.array([[a.x, a.y, a.z] for a in s.atoms])


def symbols(s):
    return [a.symbol for a in s.atoms]


def from_ase(atoms):
    s = af.Structure()
    for symbol, p in zip(atoms.get_chemical_symbols(), atoms.positions):
        s.add_atom(symbol, *map(float, p))
    if atoms.cell.rank == 3:
        s.cell = [list(map(float, row)) for row in atoms.cell]
        s.pbc = [bool(b) for b in atoms.pbc]
    return s


def to_ase(s):
    from ase import Atoms
    return Atoms(symbols(s), positions=positions(s), cell=s.cell if s.cell else None, pbc=s.pbc if s.cell else False)


def neighbour_counts(s, cutoff):
    """Neighbours of each atom within the cutoff (periodic along the periodic axes)."""
    from ase.neighborlist import neighbor_list
    i = neighbor_list("i", to_ase(s), cutoff)
    return np.bincount(i, minlength=len(s.atoms))


def bond_lengths(s, cutoff):
    from ase.neighborlist import neighbor_list
    return neighbor_list("d", to_ase(s), cutoff)


def min_intermolecular(s, labels):
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list("ijd", to_ase(s), 3.0)
    mask = labels[i] != labels[j]
    return d[mask].min() if mask.any() else math.inf


@unittest.skipIf(np is None or missing("ase"), "numpy and ase are required")
class LowDimensional(unittest.TestCase):
    def test_nanotubes_match_ase_and_bond_every_atom_three_times(self):
        from ase.build import nanotube
        for n, m in ((6, 6), (10, 0), (6, 5), (8, 4), (5, 3)):
            ours = run(f"build-nanotube {n} {m} length=2").structure
            reference = nanotube(n, m, length=2, bond=1.421)
            self.assertEqual(len(ours.atoms), len(reference), (n, m))
            # Same radius and axial period.
            p = positions(ours)
            centre = np.array(ours.cell[0][0] / 2), np.array(ours.cell[1][1] / 2)
            radius = np.hypot(p[:, 0] - centre[0], p[:, 1] - centre[1])
            ref_radius = np.hypot(*(reference.positions[:, :2] - reference.positions[:, :2].mean(axis=0)).T)
            self.assertAlmostEqual(radius.mean(), ref_radius.mean(), delta=1e-3, msg=(n, m))
            self.assertLess(np.ptp(radius), 1e-9)
            self.assertAlmostEqual(ours.cell[2][2], reference.cell[2][2], delta=1e-6)
            self.assertEqual(ours.pbc, [False, False, True])
            # Every carbon has three neighbours at 1.38-1.43 A (rolling shortens chords).
            self.assertTrue((neighbour_counts(ours, 1.6) == 3).all(), (n, m))
            d = bond_lengths(ours, 1.6)
            self.assertGreater(d.min(), 1.37)
            self.assertLess(d.max(), 1.422)

    def test_boron_nitride_and_multiwall_tubes(self):
        bn = run("build-nanotube 5 5 elements=B,N".replace(",", " ")).structure
        self.assertEqual(symbols(bn).count("B"), symbols(bn).count("N"))
        from ase.neighborlist import neighbor_list
        i, j = neighbor_list("ij", to_ase(bn), 1.6)
        self.assertTrue(all(symbols(bn)[a] != symbols(bn)[b] for a, b in zip(i, j)), "B bonds only to N")
        double = run("build-nanotube 10 0 walls=2").structure
        p = positions(double)
        radius = np.hypot(p[:, 0] - double.cell[0][0] / 2, p[:, 1] - double.cell[1][1] / 2)
        walls = sorted(set(np.round(radius, 3)))
        self.assertEqual(len(walls), 2)
        self.assertAlmostEqual(walls[1] - walls[0], 3.4, delta=0.2)
        with self.assertRaises(PipelineError):
            run("build-nanotube 6 5 walls=2")

    def test_monolayers_match_ase(self):
        from ase.build import graphene, mx2
        g = run("build-monolayer graphene na=3 nb=3").structure
        self.assertEqual(len(g.atoms), 18)
        self.assertTrue((neighbour_counts(g, 1.6) == 3).all())
        self.assertAlmostEqual(bond_lengths(g, 1.6).mean(), 1.421, delta=1e-3)
        ref = graphene("C2", a=2.4612, vacuum=7.5)
        np.testing.assert_allclose(np.array(run("build-monolayer graphene").structure.cell)[:2, :2], ref.cell[:2, :2], atol=1e-9)
        for kind, ase_kind in (("mx2-2h", "2H"), ("mx2-1t", "1T")):
            ours = run(f"build-monolayer {kind} na=2 nb=2").structure
            reference = mx2("MoS2", kind=ase_kind, a=3.18, thickness=3.19, size=(2, 2, 1), vacuum=7.5)
            self.assertEqual(sorted(symbols(ours)), sorted(reference.get_chemical_symbols()))
            # Mo has six S neighbours at sqrt(a^2/3 + t^2/4) in both phases.
            from ase.neighborlist import neighbor_list
            for structure in (to_ase(ours), reference):
                i, d = neighbor_list("id", structure, 2.6)
                mo = [k for k, s in enumerate(structure.get_chemical_symbols()) if s == "Mo"]
                for k in mo:
                    self.assertEqual((i == k).sum(), 6)
                np.testing.assert_allclose(d, math.sqrt(3.18 ** 2 / 3 + 3.19 ** 2 / 4), atol=1e-6)
        hbn = run("build-monolayer hbn").structure
        self.assertEqual(sorted(symbols(hbn)), ["B", "N"])

    def test_nanoribbons_match_ase_counts_and_are_saturated(self):
        from ase.build import graphene_nanoribbon
        for edge, width in (("armchair", 7), ("armchair", 9), ("zigzag", 4), ("zigzag", 6)):
            ours = run(f"build-nanoribbon {edge} {width} length=3").structure
            # N-AGNR counts dimer lines; ASE counts armchair widths in pairs of them.
            reference = graphene_nanoribbon(width / 2 if edge == "armchair" else width, 3, type=edge, saturated=True, C_C=1.421, C_H=1.09)
            self.assertEqual(sorted(symbols(ours)), sorted(reference.get_chemical_symbols()), (edge, width))
            # Carbon: three bonds (C or H); hydrogen: one, at 1.09 A, in the ribbon plane.
            counts = neighbour_counts(ours, 1.6)
            for k, symbol in enumerate(symbols(ours)):
                self.assertEqual(counts[k], 3 if symbol == "C" else 1, (edge, width, k))
            self.assertLess(np.ptp(positions(ours)[:, 2]), 1e-9)
            from ase.neighborlist import neighbor_list
            i, j, d = neighbor_list("ijd", to_ase(ours), 1.3)
            ch = [dist for a, b, dist in zip(i, j, d) if symbols(ours)[a] == "H"]
            np.testing.assert_allclose(ch, 1.09, atol=1e-9)
        bare = run("build-nanoribbon armchair 7 passivate=false").structure
        self.assertEqual(symbols(bare).count("H"), 0)
        resaturated = run("build-nanoribbon armchair 7 passivate=false | saturate-hydrogen").structure
        self.assertEqual(symbols(resaturated).count("H"), 4)

    def test_twisted_bilayers(self):
        for index in (1, 2, 3, 5):
            result = run(f"build-twisted-bilayer index={index}")
            s = result.structure
            theta = math.degrees(math.acos((3 * index ** 2 + 3 * index + 0.5) / (3 * index ** 2 + 3 * index + 1)))
            reported = float(result.stages[0]["notes"][0].split()[1])
            self.assertAlmostEqual(reported, theta, places=3)
            self.assertEqual(len(s.atoms), 2 * 2 * (3 * index ** 2 + 3 * index + 1))
            z = positions(s)[:, 2]
            layers = sorted(set(np.round(z, 6)))
            self.assertEqual(len(layers), 2)
            self.assertAlmostEqual(layers[1] - layers[0], 3.35, delta=1e-9)
            # Each layer on its own is perfect graphene in the moire cell: three bonds of 1.421 A.
            for level in layers:
                layer = af.Structure()
                for a in s.atoms:
                    if abs(a.z - level) < 1e-6:
                        layer.add_atom(a.symbol, a.x, a.y, a.z)
                layer.cell = s.cell
                layer.pbc = [True, True, False]
                self.assertTrue((neighbour_counts(layer, 1.6) == 3).all(), index)
                np.testing.assert_allclose(bond_lengths(layer, 1.6), 2.4612 / math.sqrt(3), atol=1e-9)
            a, b = np.array(s.cell[0]), np.array(s.cell[1])
            self.assertAlmostEqual(np.linalg.norm(a), np.linalg.norm(b), delta=1e-9)
            self.assertAlmostEqual(math.degrees(math.acos(a @ b / np.linalg.norm(a) ** 2)), 60, delta=1e-9)
        nearest = run("build-twisted-bilayer angle=1.1").stages[0]["notes"][0]
        self.assertIn("index 30", nearest)
        mos2 = run("build-twisted-bilayer mx2-2h index=1").structure
        self.assertEqual(len(mos2.atoms), 2 * 3 * 7)


@unittest.skipIf(np is None or missing("ase"), "numpy and ase are required")
class ClustersAndMolecules(unittest.TestCase):
    def radial(self, s):
        p = positions(s)
        return np.sort(np.linalg.norm(p - p.mean(axis=0), axis=1))

    def test_icosahedra_and_cuboctahedra_match_ase(self):
        from ase.cluster import Icosahedron, Octahedron
        a = 3.615
        for shells in (1, 2, 3, 4, 5):
            ours = run(f"build-cluster icosahedron Cu {shells} distance={a / math.sqrt(2)}").structure
            reference = Icosahedron("Cu", shells, latticeconstant=a)
            self.assertEqual(len(ours.atoms), len(reference))
            r_ref = np.sort(np.linalg.norm(reference.positions - reference.positions.mean(axis=0), axis=1))
            np.testing.assert_allclose(self.radial(ours), r_ref, atol=1e-6)
            cubo = run(f"build-cluster cuboctahedron Cu {shells} distance={a / math.sqrt(2)}").structure
            self.assertEqual(len(cubo.atoms), [1, 13, 55, 147, 309][shells - 1])
            if shells > 1:
                ref_cubo = Octahedron("Cu", 2 * shells - 1, cutoff=shells - 1, latticeconstant=a)
                np.testing.assert_allclose(self.radial(cubo), np.sort(np.linalg.norm(ref_cubo.positions - ref_cubo.positions.mean(axis=0), axis=1)), atol=1e-6)

    def test_decahedra(self):
        from ase.cluster import Decahedron
        counts = []
        for shells in (2, 3, 4, 5):
            ours = run(f"build-cluster decahedron Cu {shells}").structure
            counts.append(len(ours.atoms))
            # The five-fold axis and pentagon edges: every atom at least the
            # nearest-neighbour distance from the others.
            from ase.neighborlist import neighbor_list
            d = neighbor_list("d", to_ase(ours), 2.4)
            self.assertEqual(len(d), 0)
        self.assertEqual(counts, [7, 23, 54, 105])
        self.assertEqual(len(Decahedron("Cu", 2, 1, 0)), 7)
        self.assertEqual(len(Decahedron("Cu", 3, 1, 0)), 23)
        with self.assertRaises(PipelineError):
            run("build-cluster decahedron Cu 1")

    def test_molecule_geometries(self):
        def distance(s, i, j):
            return np.linalg.norm(positions(s)[i] - positions(s)[j])

        def angle(s, i, j, k):
            p = positions(s)
            u, v = p[i] - p[j], p[k] - p[j]
            return math.degrees(math.acos(u @ v / np.linalg.norm(u) / np.linalg.norm(v)))
        water = run("build-molecule water").structure
        self.assertEqual(symbols(water), ["O", "H", "H"])
        self.assertAlmostEqual(distance(water, 0, 1), 0.9572, places=6)
        self.assertAlmostEqual(angle(water, 1, 0, 2), 104.52, places=6)
        ammonia = run("build-molecule NH3").structure
        for a, b in ((1, 2), (1, 3), (2, 3)):
            self.assertAlmostEqual(angle(ammonia, a, 0, b), 106.67, places=6)
        methane = run("build-molecule CH4").structure
        for a in range(1, 5):
            self.assertAlmostEqual(distance(methane, 0, a), 1.087, places=6)
        self.assertAlmostEqual(angle(methane, 1, 0, 2), 109.4712, places=3)
        ethane = run("build-molecule C2H6").structure
        self.assertAlmostEqual(distance(ethane, 0, 1), 1.535, places=6)
        for h in range(2, 8):
            self.assertAlmostEqual(min(distance(ethane, 0, h), distance(ethane, 1, h)), 1.094, places=6)
        benzene = run("build-molecule benzene").structure
        self.assertTrue((neighbour_counts(benzene, 1.5) == [3] * 6 + [1] * 6).all())
        c60 = run("build-molecule C60").structure
        self.assertEqual(len(c60.atoms), 60)
        self.assertTrue((neighbour_counts(c60, 1.6) == 3).all())
        np.testing.assert_allclose(self.radial(c60), self.radial(c60)[0], atol=1e-9)  # all on one sphere
        # Centre of mass at the box centre.
        boxed = run("build-molecule CO2 vacuum=10").structure
        np.testing.assert_allclose(positions(boxed)[0], np.array(boxed.cell).diagonal() / 2, atol=1e-9)
        with self.assertRaises(PipelineError):
            run("build-molecule unobtainium")
        self.assertEqual(symbols(run("build-molecule Fe").structure), ["Fe"])

    def test_packing_respects_distance_density_and_geometry(self):
        result = run("build-packing H2O 120 CO2 10 density=0.95 min_distance=2.0 seed=3")
        s = result.structure
        self.assertEqual(symbols(s).count("O"), 120 + 20)
        self.assertEqual(symbols(s).count("C"), 10)
        mass = 120 * 18.015 + 10 * 44.009
        volume = abs(np.linalg.det(np.array(s.cell)))
        self.assertAlmostEqual(mass * 1.66053906660 / volume, 0.95, delta=0.002)
        # Molecules: CO2 first (larger), then water; labels by position in the output.
        labels = np.array([k // 3 for k in range(len(s.atoms))])
        self.assertGreaterEqual(min_intermolecular(s, labels), 2.0 - 1e-9)
        # Water (listed first, as large as CO2) comes first.
        self.assertEqual(symbols(s)[:3], ["O", "H", "H"])
        for k in range(120):
            o = positions(s)[3 * k]
            for h in (1, 2):
                self.assertAlmostEqual(np.linalg.norm(positions(s)[3 * k + h] - o), 0.9572, places=6)
        again = run("build-packing H2O 120 CO2 10 density=0.95 min_distance=2.0 seed=3").structure
        np.testing.assert_allclose(positions(again), positions(s))
        with self.assertRaises(PipelineError):
            run("build-packing H2O 200 box=10,10,10 min_distance=2.5 attempts=200")


@unittest.skipIf(np is None or missing("ase"), "numpy and ase are required")
class SurfacesAndDefects(unittest.TestCase):
    def test_adsorption_sites_match_ase(self):
        from ase.build import fcc111
        slab = fcc111("Pt", size=(3, 3, 4), vacuum=10, a=3.92)
        info = slab.info["adsorbate_info"]
        cell2 = info["cell"]
        top = max(slab.positions[:, 2])
        for site, ase_site in (("top", "ontop"), ("bridge", "bridge"), ("fcc", "fcc"), ("hcp", "hcp")):
            result = run(f"add-adsorbate O site={site} height=1.5 position=0.5,0.5,0", from_ase(slab))
            o = positions(result.structure)[-1]
            self.assertAlmostEqual(o[2], top + 1.5, places=9)
            # The site is one of ASE's sites of that kind (modulo the surface lattice).
            lattice = np.array(slab.cell[:2, :2]) / 3  # the 1 x 1 surface cell
            frac = np.linalg.solve(lattice.T, o[:2] - slab.positions[np.argmax(slab.positions[:, 2]), :2])
            if site == "bridge":
                # A midpoint of two top atoms (three equivalent orientations).
                np.testing.assert_allclose(2 * frac, np.round(2 * frac), atol=1e-6)
                self.assertGreater(np.abs(frac - np.round(frac)).max(), 0.4)
            else:
                reference = np.array(info["sites"][ase_site]) @ cell2
                d = np.linalg.solve(lattice.T, o[:2] - reference)
                np.testing.assert_allclose(d, np.round(d), atol=1e-6, err_msg=site)
            self.assertEqual(result.selected_indices, [len(slab)])
        co = run("add-adsorbate CO site=top height=1.85", from_ase(slab)).structure
        c, o = positions(co)[-2], positions(co)[-1]
        self.assertEqual(symbols(co)[-2:], ["C", "O"])
        self.assertAlmostEqual(c[2], top + 1.85, places=9)
        np.testing.assert_allclose(o - c, [0, 0, 1.128], atol=1e-9)

    def test_antisite_and_frenkel_pairs(self):
        from ase.build import bulk
        b2 = bulk("NiAl", "cesiumchloride", a=2.88, cubic=True).repeat((3, 3, 3))
        result = run("create-antisite Ni Al seed=4", from_ase(b2))
        s = result.structure
        self.assertEqual(sorted(symbols(s)), sorted(b2.get_chemical_symbols()))
        changed = [k for k in range(len(b2)) if symbols(s)[k] != b2.get_chemical_symbols()[k]]
        self.assertEqual(len(changed), 2)
        self.assertEqual(result.selected_indices, sorted(changed))
        d = to_ase(s).get_distance(changed[0], changed[1], mic=True)
        self.assertAlmostEqual(d, 2.88 * math.sqrt(3) / 2, places=9)  # nearest Ni-Al
        cu = bulk("Cu", "fcc", a=3.615, cubic=True).repeat((3, 3, 3))
        result = run("create-frenkel-pair Cu min_distance=5 seed=2", from_ase(cu))
        s = result.structure
        self.assertEqual(len(s.atoms), len(cu))
        new = positions(s)[-1]
        # The interstitial sits in a void: no atom closer than 1.5 A.
        atoms = to_ase(s)
        self.assertGreater(min(atoms.get_distances(len(s.atoms) - 1, range(len(s.atoms) - 1), mic=True)), 1.5)
        # The vacancy: the lattice site no atom occupies, at least 5 A from the interstitial.
        occupied = to_ase(s)
        vacancies = [p for p in cu.positions
                     if min(np.linalg.norm(occupied.positions[:-1] - p, axis=1)) > 0.1]
        self.assertEqual(len(vacancies), 1)
        delta = vacancies[0] - new
        frac = np.linalg.solve(cu.cell.T, delta)
        frac -= np.round(frac)
        self.assertGreaterEqual(np.linalg.norm(frac @ cu.cell), 5 - 1e-9)

    def test_saturating_a_diamond_surface(self):
        from ase.build import diamond100
        slab = diamond100("C", size=(2, 2, 6), vacuum=8, a=3.567)
        result = run("saturate-hydrogen", from_ase(slab))
        s = result.structure
        # Each surface carbon (two bonds) gets two hydrogens at 1.09 A, H-C-H tetrahedral.
        hydrogens = [k for k, x in enumerate(symbols(s)) if x == "H"]
        self.assertEqual(len(hydrogens), 2 * 2 * 4)  # two surfaces of 4 atoms, two H each
        self.assertEqual(result.selected_indices, hydrogens)
        # Bonds: C-C up to 1.7 A, C-H up to 1.2 A (neighbouring dihydrides crowd their H atoms).
        from ase.neighborlist import neighbor_list
        i, j, d = neighbor_list("ijd", to_ase(s), 1.7)
        bonded = [(a, b) for a, b, dist in zip(i, j, d)
                  if (symbols(s)[a] == "C" and symbols(s)[b] == "C") or ("H" in (symbols(s)[a], symbols(s)[b]) and "C" in (symbols(s)[a], symbols(s)[b]) and dist < 1.2)]
        counts = np.bincount([a for a, b in bonded], minlength=len(s.atoms))
        for k, x in enumerate(symbols(s)):
            self.assertEqual(counts[k], 4 if x == "C" else 1, k)
        atoms = to_ase(s)
        for h in hydrogens:
            carbon = int(np.argmin([atoms.get_distance(h, k, mic=True) if symbols(s)[k] == "C" else 9 for k in range(len(s.atoms))]))
            self.assertAlmostEqual(atoms.get_distance(h, carbon, mic=True), 1.09, places=9)
        # Bulk silicon needs nothing; a lone carbon becomes methane, a lone oxygen water.
        from ase.build import bulk
        si = run("saturate-hydrogen", from_ase(bulk("Si", "diamond", a=5.43, cubic=True)))
        self.assertIn("0 hydrogen", si.stages[0]["notes"][0])
        methane = run("build-molecule C vacuum=0 | saturate-hydrogen").structure
        self.assertEqual(sorted(symbols(methane)), ["C", "H", "H", "H", "H"])
        p = positions(methane)
        angles = [math.degrees(math.acos((p[a] - p[0]) @ (p[b] - p[0]) / 1.09 ** 2)) for a in range(1, 5) for b in range(a + 1, 5)]
        np.testing.assert_allclose(angles, 109.4712, atol=1e-3)
        water = run("build-molecule O vacuum=0 | saturate-hydrogen").structure
        p = positions(water)
        self.assertAlmostEqual(math.degrees(math.acos((p[1] - p[0]) @ (p[2] - p[0]) / 0.97 ** 2)), 104.5, places=6)


@unittest.skipIf(np is None or missing("ase"), "numpy and ase are required")
class GeometryEdits(unittest.TestCase):
    def copper(self, repeat=2):
        from ase.build import bulk
        return from_ase(bulk("Cu", "fcc", a=3.615, cubic=True).repeat((repeat,) * 3))

    def test_rigid_motions_preserve_distances(self):
        s = self.copper()
        reference = to_ase(s).get_all_distances(mic=True)
        for step in ("mirror normal=1,2,3", "invert", "rotate axis=1,1,0 angle=33", "rotate axis=0,0,1 angle=90 center=cell",
                     "align-direction 1 1 1 axis=z", "align-direction 1 1 1 axis=z second_direction=1,-1,0 second_axis=x"):
            out = run(step, s).structure
            np.testing.assert_allclose(to_ase(out).get_all_distances(mic=True), reference, atol=1e-9, err_msg=step)
        twice = to_ase(run("mirror normal=1,2,3 center=origin | mirror normal=1,2,3 center=origin", s).structure)
        # Back on the same sites (up to lattice vectors) in a cell spanning the same lattice.
        original = to_ase(s)
        for k in range(len(original)):
            d = twice.positions[k] - original.positions[k]
            frac = np.linalg.solve(original.cell.T, d)
            np.testing.assert_allclose(frac, np.round(frac), atol=1e-9)
        self.assertAlmostEqual(abs(np.linalg.det(twice.cell)), abs(np.linalg.det(original.cell)), places=9)
        self.assertGreater(np.linalg.det(to_ase(run("mirror normal=1,2,3", s).structure).cell), 0)
        # Rotating with the cell keeps the fractional coordinates.
        turned = run("rotate axis=1,2,3 angle=40", s).structure
        np.testing.assert_allclose(to_ase(turned).get_scaled_positions(), to_ase(s).get_scaled_positions(), atol=1e-9)

    def test_align_direction(self):
        s = self.copper(1)
        aligned = run("align-direction 1 1 1 axis=z second_direction=1,-1,0 second_axis=x", s).structure
        c = np.array(aligned.cell)
        body = c[0] + c[1] + c[2]
        np.testing.assert_allclose(body / np.linalg.norm(body), [0, 0, 1], atol=1e-12)
        second = c[0] - c[1]
        np.testing.assert_allclose(second / np.linalg.norm(second), [1, 0, 0], atol=1e-12)
        cartesian = run("align-direction 0 1 0 lattice=false axis=x", s).structure
        np.testing.assert_allclose(np.array(cartesian.cell)[1], [3.615, 0, 0], atol=1e-12)

    def test_set_and_measure_geometry(self):
        ethane = "build-molecule C2H6 vacuum=0"
        for target in (60.0, -75.0, 180.0, 10.0):
            result = run(f"{ethane} | set-geometry 2 0 1 5 value={target} | measure-geometry 2 0 1 5")
            value = float(result.stages[-1]["notes"][0].split("=")[1].split()[0])
            self.assertAlmostEqual(math.remainder(value - target, 360), 0, places=4)
        before = run(ethane).structure
        after = run(f"{ethane} | set-geometry 2 0 1 5 value=60").structure
        # The methyl group on carbon 1 turned rigidly: distances within it unchanged; carbon 0's side fixed.
        right = [1, 5, 6, 7]
        pa, pb = positions(before), positions(after)
        for i in right:
            for j in right:
                self.assertAlmostEqual(np.linalg.norm(pa[i] - pa[j]), np.linalg.norm(pb[i] - pb[j]), places=9)
        np.testing.assert_allclose(pa[[0, 2, 3, 4]], pb[[0, 2, 3, 4]], atol=1e-12)
        water = run("build-molecule H2O vacuum=0 | set-geometry 1 0 2 value=120 | measure-geometry 1 0 2")
        self.assertIn("= 120", water.stages[-1]["notes"][0])
        co = run("build-molecule CO vacuum=0 | set-geometry 0 1 value=1.5 move=atom | measure-geometry 0 1")
        self.assertIn("= 1.5 ", co.stages[-1]["notes"][0] + " ")
        with self.assertRaises(PipelineError):
            run("build-molecule benzene vacuum=0 | set-geometry 0 1 value=1.6")  # a ring
        # Across a periodic boundary: minimum-image distances.
        cu = self.copper(1)
        self.assertIn("= 2.556", run("measure-geometry 0 1", cu).stages[0]["notes"][0])

    def test_sort_duplicate_swap_scale(self):
        s = run("build-molecule C2H6 vacuum=0").structure
        by_z = run("sort-atoms key=z", s).structure
        self.assertEqual(list(positions(by_z)[:, 2]), sorted(positions(s)[:, 2]))
        by_element = run("sort-atoms key=element descending=true", s).structure
        self.assertEqual(symbols(by_element), ["C", "C"] + ["H"] * 6)
        copies = run("select-element C | duplicate-selected offset=0,0,5 copies=2", s)
        self.assertEqual(len(copies.structure.atoms), 12)
        np.testing.assert_allclose(positions(copies.structure)[8:10], positions(s)[:2] + [0, 0, 5])
        np.testing.assert_allclose(positions(copies.structure)[10:12], positions(s)[:2] + [0, 0, 10])
        self.assertEqual(copies.selected_indices, list(range(8, 12)))
        from ase.build import bulk
        nial = from_ase(bulk("NiAl", "cesiumchloride", a=2.88).repeat((2, 2, 2)))
        swapped = run("swap-elements Ni Al", nial).structure
        self.assertEqual(symbols(swapped), ["Al" if x == "Ni" else "Ni" for x in symbols(nial)])
        cu = self.copper(1)
        dense = run("scale-cell mode=density value=9.0", cu).structure
        mass = 4 * 63.546
        self.assertAlmostEqual(mass * 1.66053906660 / abs(np.linalg.det(np.array(dense.cell))), 9.0, places=9)
        np.testing.assert_allclose(to_ase(dense).get_scaled_positions(), to_ase(cu).get_scaled_positions(), atol=1e-12)
        volume = run("scale-cell mode=volume value=50", cu).structure
        self.assertAlmostEqual(abs(np.linalg.det(np.array(volume.cell))), 50, places=9)

    def test_displace_fractional_and_wrap(self):
        cu = self.copper(1)
        moved = run("displace 0.25 0.5 0.75 fractional=true wrap=true", cu).structure
        frac = to_ase(moved).get_scaled_positions(wrap=False)
        expected = (to_ase(cu).get_scaled_positions() + [0.25, 0.5, 0.75]) % 1.0
        np.testing.assert_allclose(frac, expected, atol=1e-12)
        self.assertTrue(((frac >= 0) & (frac < 1)).all())

    def test_niggli_reduction_matches_ase(self):
        from ase.build import bulk, niggli_reduce
        rng = np.random.default_rng(7)
        for name, atoms in (("fcc", bulk("Cu", "fcc", a=3.615)), ("bcc", bulk("Fe", "bcc", a=2.87)),
                            ("hcp", bulk("Mg", "hcp", a=3.21, c=5.21)), ("tri", bulk("Cu", "fcc", a=3.615))):
            if name == "tri":
                atoms.set_cell(np.array([[3.1, 0.2, 0.1], [0.7, 4.3, 0.4], [1.1, 0.9, 5.7]]), scale_atoms=True)
            for trial in range(3):
                # A random unimodular basis of the same lattice.
                while True:
                    m = rng.integers(-2, 3, size=(3, 3))
                    if round(abs(np.linalg.det(m))) == 1:
                        break
                skewed = atoms.copy()
                skewed.set_cell(m @ atoms.cell, scale_atoms=False)
                skewed.wrap()
                ours = run("niggli-reduce", from_ase(skewed)).structure
                reference = skewed.copy()
                niggli_reduce(reference)
                np.testing.assert_allclose(to_ase(ours).cell.cellpar(), reference.cell.cellpar(), atol=1e-6, err_msg=f"{name} {trial}")
                self.assertAlmostEqual(abs(np.linalg.det(np.array(ours.cell))), atoms.get_volume(), places=6)
                frac = to_ase(ours).get_scaled_positions(wrap=False)
                self.assertTrue(((frac > -1e-9) & (frac < 1 + 1e-9)).all())

    def test_near_cubic_supercells(self):
        from ase.build import bulk
        for atoms, target, edge in ((bulk("Cu", "fcc", a=3.615), 4, 3.615), (bulk("Cu", "fcc", a=3.615), 32, 7.23),
                                    (bulk("Fe", "bcc", a=2.87), 2, 2.87), (bulk("Fe", "bcc", a=2.87), 16, 5.74)):
            result = run(f"near-cubic-supercell target_atoms={target}", from_ase(atoms))
            cell = to_ase(result.structure).cell.cellpar()
            np.testing.assert_allclose(cell, [edge] * 3 + [90] * 3, atol=1e-9)
            self.assertEqual(len(result.structure.atoms), target)
            self.assertIn("deviation from a cube", result.stages[0]["notes"][0])
            # Every lattice site once: no two atoms coincide.
            self.assertGreater(to_ase(result.structure).get_all_distances(mic=True)[np.triu_indices(target, 1)].min(), 2.0)

    def test_interpolation_takes_the_short_way(self):
        s = af.Structure()
        s.add_atom("Cu", 0.5, 5, 5)
        s.add_atom("Ni", 5, 5, 5)
        s.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        end = s.copy()
        end.atoms[0].x = 9.5   # 1 A to the left across the boundary
        end.atoms[1].x = 7.0
        with tempfile.TemporaryDirectory() as folder:
            path = str(Path(folder) / "end.xyz")
            af.save(end, path)
            for t, cu_x, ni_x in ((0.0, 0.5, 5.0), (0.5, 0.0, 6.0), (1.0, -0.5, 7.0)):
                mid = run(f"interpolate {path} fraction={t}", s).structure
                self.assertAlmostEqual(mid.atoms[0].x, cu_x, places=9)
                self.assertAlmostEqual(mid.atoms[1].x, ni_x, places=9)
            wrong = s.copy()
            wrong.add_atom("H", 1, 1, 1)
            with self.assertRaises(PipelineError):
                run(f"interpolate {path}", wrong)

    def test_unwrap_and_select_fragments(self):
        s = af.Structure()
        s.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        # A water molecule split across x = 0, and a separate CO molecule.
        s.add_atom("O", 0.1, 5, 5)
        s.add_atom("H", 9.4, 5.6, 5)
        s.add_atom("H", 9.4, 4.4, 5)
        s.add_atom("C", 5, 5, 2)
        s.add_atom("O", 5, 5, 3.128)
        whole = run("unwrap-molecules", s).structure
        p = positions(whole)
        self.assertAlmostEqual(np.linalg.norm(p[1] - p[0]), math.hypot(0.7, 0.6), places=9)
        self.assertLess(np.linalg.norm(p[1] - p[0]), 1.0)
        touching = run("select-expression index == 1 | select-fragment", s)
        self.assertEqual(touching.selected_indices, [0, 1, 2])
        small = run("select-fragment rule=size min_atoms=1 max_atoms=2", s)
        self.assertEqual(small.selected_indices, [3, 4])

    def test_measure_and_set_geometry_need_valid_indices(self):
        with self.assertRaises(PipelineError):
            run("build-molecule H2O | measure-geometry 0 7")
        with self.assertRaises(PipelineError):
            run("build-molecule H2O | set-geometry 0 0 value=1")



def notes_of(result):
    return [n for stage in result.stages for n in stage.get("notes", [])]


@unittest.skipIf(np is None or missing("ase"), "numpy and ase are required")
class AdvancedBuilderOptions(unittest.TestCase):
    def test_octahedra_match_ase_and_nanoalloys(self):
        from ase.cluster import Octahedron
        for length in (2, 3, 4, 5):
            ours = run(f"build-cluster kind=octahedron element=Au shells={length} distance=2.884").structure
            theirs = Octahedron("Au", length, 0, latticeconstant=2.884 * math.sqrt(2))
            self.assertEqual(len(ours.atoms), len(theirs))
            a = np.sort(np.linalg.norm(positions(ours) - positions(ours).mean(axis=0), axis=1))
            b = np.sort(np.linalg.norm(theirs.positions - theirs.positions.mean(axis=0), axis=1))
            np.testing.assert_allclose(a, b, atol=1e-9)
        radius = lambda s: np.linalg.norm(positions(s) - positions(s).mean(axis=0), axis=1)
        for distribution in ("random", "core", "shell", "janus"):
            s = run(f"build-cluster element=Pt shells=4 second_element=Au fraction=0.3 distribution={distribution} seed=2").structure
            gold = np.array(symbols(s)) == "Au"
            self.assertEqual(gold.sum(), round(0.3 * 147))
            r = radius(s)
            if distribution == "core":
                self.assertLessEqual(r[gold].max(), r[~gold].min() + 1e-9)
            if distribution == "shell":
                self.assertGreaterEqual(r[gold].min(), r[~gold].max() - 1e-9)
            if distribution == "janus":
                z = positions(s)[:, 2]
                self.assertGreaterEqual(z[gold].min(), z[~gold].max() - 1e-9)

    def test_molecule_orientation_and_files(self):
        co2 = positions(run("build-molecule CO2 orientation=principal vacuum=0").structure)
        np.testing.assert_allclose(co2[:, :2], 0, atol=1e-9)  # the long axis on z
        benzene = positions(run("build-molecule C6H6 orientation=principal vacuum=0").structure)
        np.testing.assert_allclose(benzene[:, 1], 0, atol=1e-9)  # planar: in xz
        shuffled = positions(run("build-molecule C6H6 orientation=random seed=5 vacuum=0").structure)
        # A random orientation leaves the benzene plane (it lay in xz above).
        self.assertGreater(np.ptp(shuffled[:, 1]), 0.1)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "ammonia.xyz"
            path.write_text("4\nNH3\nN 0 0 0\nH 0 0.94 0.38\nH 0.81 -0.47 0.38\nH -0.81 -0.47 0.38\n", encoding="utf-8")
            mol = run(f"build-molecule {path} vacuum=6 periodic=true").structure
            self.assertEqual(symbols(mol), ["N", "H", "H", "H"])
            self.assertTrue(all(mol.pbc))
            packed = run(f"build-packing {path} 6 H2O 6 box=12,12,12 seed=3").structure
            self.assertEqual(symbols(packed).count("N"), 6)

    def test_packing_regions_and_solvation(self):
        sphere = run("build-packing H2O 30 box=20,20,20 region=sphere region_center=10,10,10 region_radius=7").structure
        self.assertLessEqual(np.linalg.norm(positions(sphere) - 10, axis=1).max(), 7 + 1e-9)
        outside = run("build-packing H2O 20 box=16,16,16 region=outside-sphere region_center=8,8,8 region_radius=5").structure
        self.assertGreaterEqual(np.linalg.norm(positions(outside) - 8, axis=1).min(), 5 - 1e-9)
        slab = run("build-packing H2O 0 box=15,15,15 density=1.0 region=slab region_low=0,0,5 region_high=0,0,10").structure
        z = positions(slab)[:, 2]
        self.assertTrue(z.min() >= 5 - 1e-9 and z.max() <= 10 + 1e-9)
        self.assertAlmostEqual(len(slab.atoms) // 3, 15 * 15 * 5 * 1.0 / (18.015 * 1.66053906660), delta=1)
        tube = run("build-packing CH4 10 box=20,20,20 region=cylinder region_center=10,10,0 region_radius=4 region_axis=z").structure
        self.assertLessEqual(np.linalg.norm(positions(tube)[:, :2] - 10, axis=1).max(), 4 + 1e-9)
        solvated = run("build-molecule C60 vacuum=10 periodic=true | build-packing H2O 0 mode=solvate density=1.0 solute_distance=3")
        s = solvated.structure
        self.assertEqual(symbols(s)[:60], ["C"] * 60)
        waters = (len(s.atoms) - 60) // 3
        self.assertGreater(waters, 120)
        carbon, rest = positions(s)[:60], positions(s)[60:]
        cell = np.array(s.cell)
        d = rest[:, None, :] - carbon[None, :, :]
        frac = d @ np.linalg.inv(cell)
        d = (frac - np.round(frac)) @ cell
        self.assertGreaterEqual(np.linalg.norm(d, axis=2).min(), 3 - 1e-9)
        self.assertEqual(len(solvated.selected_indices), len(s.atoms) - 60)
        with self.assertRaises(PipelineError):
            run("build-packing H2O 10 mode=solvate")  # no solute cell

    def test_adsorbate_coverage_orientation_and_both_surfaces(self):
        from ase.build import fcc111
        slab = fcc111("Pt", size=(4, 4, 4), vacuum=10, a=3.92)
        top, bottom = slab.positions[:, 2].max(), slab.positions[:, 2].min()
        full = run("add-adsorbate O site=fcc coverage=all height=1.2", from_ase(slab)).structure
        self.assertEqual(len(full.atoms) - len(slab), 16)  # one fcc hollow per surface atom
        bridges = run("add-adsorbate O site=bridge coverage=all", from_ase(slab)).structure
        self.assertEqual(len(bridges.atoms) - len(slab), 48)
        both = run("add-adsorbate O site=top coverage=all surface=both height=2", from_ase(slab)).structure
        z = positions(both)[len(slab):, 2]
        self.assertEqual(sorted(set(np.round(z, 9))), [round(bottom - 2, 9), round(top + 2, 9)])
        sparse = run("add-adsorbate O site=top coverage=fraction fraction=0.25 min_separation=5 seed=3", from_ase(slab)).structure
        o = positions(sparse)[len(slab):]
        self.assertEqual(len(o), 4)
        cell = np.array(slab.cell)
        for i in range(len(o)):
            for j in range(i + 1, len(o)):
                d = o[j] - o[i]
                f = d @ np.linalg.inv(cell)
                f[:2] -= np.round(f[:2])
                self.assertGreaterEqual(np.linalg.norm((f @ cell)[:2]), 5 - 1e-9)
        tilted = positions(run("add-adsorbate CO site=top height=1.85 tilt=30 azimuth=90", from_ase(slab)).structure)
        axis = tilted[-1] - tilted[-2]
        self.assertAlmostEqual(math.degrees(math.acos(axis[2] / np.linalg.norm(axis))), 30, places=6)
        self.assertAlmostEqual(axis[0], 0, places=9)  # azimuth 90: tilted towards y
        anchored = positions(run("add-adsorbate CO site=top height=2 anchor=1", from_ase(slab)).structure)
        self.assertAlmostEqual(anchored[-1][2], top + 2, places=9)  # O (atom 1) at the height


@unittest.skipIf(np is None or missing("ase"), "numpy and ase are required")
class AdvancedEditOptions(unittest.TestCase):
    def test_euler_rotation_and_sorting(self):
        water = run("build-molecule H2O vacuum=0").structure
        turned = run("rotate mode=euler euler=30,40,50 convention=zyz center=origin", water).structure
        rz = lambda a: np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]])
        ry = lambda a: np.array([[math.cos(a), 0, math.sin(a)], [0, 1, 0], [-math.sin(a), 0, math.cos(a)]])
        r = rz(math.radians(30)) @ ry(math.radians(40)) @ rz(math.radians(50))
        np.testing.assert_allclose(positions(turned), positions(water) @ r.T, atol=1e-12)
        mixed = run("build-packing H2O 4 CO2 2 box=10,10,10 | sort-atoms key=element element_order=O,C").structure
        self.assertEqual("".join(symbols(mixed)), "O" * 8 + "C" * 2 + "H" * 8)
        layered = run("build-packing H2O 6 box=10,10,10 | sort-atoms key=element then=z").structure
        for block in ("H", "O"):
            z = [p[2] for p, e in zip(positions(layered), symbols(layered)) if e == block]
            self.assertEqual(z, sorted(z))

    def test_saturation_options(self):
        from ase.build import diamond100
        slab = from_ase(diamond100("C", size=(2, 2, 4), vacuum=8))
        fluorinated = run("saturate-hydrogen terminator=F", slab).structure
        self.assertEqual(symbols(fluorinated).count("F"), 16)
        p = positions(fluorinated)
        f = p[np.array(symbols(fluorinated)) == "F"]
        c = p[np.array(symbols(fluorinated)) == "C"]
        self.assertAlmostEqual(np.linalg.norm(f[:, None] - c[None], axis=2).min(axis=1).max(), 0.76 + 0.57, delta=0.02)
        one = run("saturate-hydrogen max_per_atom=1", slab).structure
        self.assertEqual(symbols(one).count("H"), 8)
        long = run("saturate-hydrogen bond_length=1.2", slab).structure
        h = positions(long)[np.array(symbols(long)) == "H"]
        np.testing.assert_allclose(np.linalg.norm(h[:, None] - c[None], axis=2).min(axis=1), 1.2, atol=1e-9)
        top = max(a.z for a in slab.atoms)
        only_top = run(f"select-expression z>{top - 0.1} | saturate-hydrogen target=selected", slab).structure
        self.assertEqual(symbols(only_top).count("H"), 8)

    def test_set_geometry_moving_both_sides(self):
        from ase.build import molecule
        ethane = from_ase(molecule("C2H6"))
        stretched = run("set-geometry 0 1 value=1.7 move=both", ethane).structure
        p, q = positions(stretched), positions(ethane)
        self.assertAlmostEqual(np.linalg.norm(p[0] - p[1]), 1.7, places=9)
        np.testing.assert_allclose((p[0] + p[1]) / 2, (q[0] + q[1]) / 2, atol=1e-12)
        twisted = run("set-geometry 2 0 1 5 value=60 move=both | measure-geometry 2 0 1 5", ethane)
        self.assertIn("60", notes_of(twisted)[-1])

    def test_per_axis_scaling(self):
        from ase.build import bulk
        cu = from_ase(bulk("Cu", "fcc", a=3.6, cubic=True))
        s = run("scale-cell mode=axes factors=1.02,1.02,0.97", cu).structure
        np.testing.assert_allclose(np.diag(s.cell), [3.672, 3.672, 3.492], atol=1e-12)
        frac = positions(s) @ np.linalg.inv(np.array(s.cell))
        np.testing.assert_allclose(frac, positions(cu) / 3.6, atol=1e-12)
        fixed = run("scale-cell mode=factor value=1.1 scale_atoms=false", cu).structure
        np.testing.assert_allclose(positions(fixed), positions(cu))


@unittest.skipIf(np is None or missing("ase"), "numpy and ase are required")
class AdvancedDefects(unittest.TestCase):
    def test_near_cubic_shapes_and_tolerance(self):
        from ase.build import bulk
        cu = from_ase(bulk("Cu", "fcc", a=3.6))
        sc = run("near-cubic-supercell target_atoms=32 shape=sc", cu).structure
        metric = np.array(sc.cell) @ np.array(sc.cell).T
        np.testing.assert_allclose(metric, 7.2 ** 2 * np.eye(3), atol=1e-9)
        image = run("near-cubic-supercell target_atoms=32 shape=image", cu)
        self.assertIn("shortest image 7.2 ", notes_of(image)[0])
        loose = run("near-cubic-supercell target_atoms=30 tolerance=20", cu).structure
        self.assertEqual(len(loose.atoms), 32)  # 32 makes a perfect cube; 30 cannot

    def test_idpp_keeps_atoms_apart(self):
        from ase.build import bulk
        from ase.io import write
        start = bulk("Cu", "fcc", a=3.6, cubic=True) * (2, 2, 2)
        del start[0]
        end = start.copy()
        end.positions[0] = [0, 0, 0]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "end.xyz"
            write(str(path), end, format="extxyz")
            closest = {}
            for method in ("linear", "idpp"):
                p = positions(run(f"interpolate file={path} method={method} fraction=0.5", from_ase(start)).structure)
                closest[method] = np.linalg.norm(p[0] - p[1:], axis=1).min()
            self.assertGreater(closest["idpp"], closest["linear"] + 0.05)

    def test_multiple_antisites_and_dumbbells(self):
        from ase.build import bulk
        nial = from_ase(bulk("NiAl", "cesiumchloride", a=2.88) * (4, 4, 4))
        result = run("create-antisite Ni Al count=3 min_separation=6 seed=2", nial)
        swapped = result.selected_indices
        self.assertEqual(len(swapped), 6)
        self.assertEqual(symbols(result.structure).count("Ni"), 64)
        cu = from_ase(bulk("Cu", "fcc", a=3.6, cubic=True) * (3, 3, 3))
        result = run("create-frenkel-pair element=Cu interstitial=dumbbell direction=1,0,0 min_distance=5 max_distance=8", cu)
        s = result.structure
        self.assertEqual(len(s.atoms), 108)
        a, b = positions(s)[result.selected_indices]
        np.testing.assert_allclose(np.abs(b - a), [0.9 * 3.6 / math.sqrt(2), 0, 0], atol=1e-9)
        with self.assertRaises(PipelineError):
            run("create-frenkel-pair element=Cu interstitial=dumbbell min_distance=0 max_distance=1", cu)



def shape_step(script, structure, **parameters):
    from atomforge.pipeline import Step
    return Pipeline([Step("build-shape", {"script": script, **parameters})]).run(structure)


@unittest.skipIf(np is None or missing("ase"), "numpy and ase are required")
class ShapeModeller(unittest.TestCase):
    """build-shape against the same carving done independently with numpy."""

    @classmethod
    def setUpClass(cls):
        from ase.build import bulk
        cls.cell = bulk("Cu", "fcc", a=3.615, cubic=True)
        cls.crystal = from_ase(cls.cell)
        # The lattice around the origin, far enough for every shape below.
        big = cls.cell * (24, 24, 24)
        big.translate(-12 * 3.615 * np.ones(3))
        cls.points = big.positions

    def carved(self, script, **parameters):
        result = shape_step(script, self.crystal, **parameters)
        s = result.structure
        self.assertEqual(s.pbc, [False, False, False])
        return result, positions(s)

    def assertSameAtoms(self, script, inside, **parameters):
        result, p = self.carved(script, vacuum=0, **parameters)
        expected = self.points[inside(self.points)]
        self.assertEqual(len(p), len(expected), script)
        # Up to the box shift, the atoms are the same.
        shift = expected.min(axis=0) - p.min(axis=0)
        mine = np.round(p + shift, 6)
        theirs = np.round(expected, 6)
        self.assertEqual(sorted(map(tuple, mine)), sorted(map(tuple, theirs)), script)
        return result, p

    def test_primitives_match_independent_carving(self):
        r = np.linalg.norm
        self.assertSameAtoms("sphere(r = 14.9);", lambda q: r(q, axis=1) <= 14.9)
        self.assertSameAtoms("sphere(d = 20.1);", lambda q: r(q, axis=1) <= 10.05)
        self.assertSameAtoms("box(size = [20.1, 12.3, 8.7]);", lambda q: np.all(np.abs(q) <= [10.05, 6.15, 4.35], axis=1))
        self.assertSameAtoms("box(size = [10.1, 10.1, 10.1], center = false);", lambda q: np.all((q >= 0) & (q <= 10.1), axis=1))
        self.assertSameAtoms("cylinder(r = 9.1, h = 20.3);", lambda q: (r(q[:, :2], axis=1) <= 9.1) & (np.abs(q[:, 2]) <= 10.15))
        self.assertSameAtoms("cone(r1 = 10.1, r2 = 2.1, h = 16.1);",
                             lambda q: (np.abs(q[:, 2]) <= 8.05) & (r(q[:, :2], axis=1) <= 10.1 + (2.1 - 10.1) * (q[:, 2] + 8.05) / 16.1 + 1e-9))
        self.assertSameAtoms("torus(R = 12.1, r = 4.1);", lambda q: np.hypot(r(q[:, :2], axis=1) - 12.1, q[:, 2]) <= 4.1)
        self.assertSameAtoms("capsule(r = 5.1, h = 12.1);",
                             lambda q: r(q - np.stack([0 * q[:, 0], 0 * q[:, 1], np.clip(q[:, 2], -6.05, 6.05)], axis=1), axis=1) <= 5.1)
        self.assertSameAtoms("ellipsoid(r = [12.1, 8.1, 6.1]);", lambda q: r(q / [12.1, 8.1, 6.1], axis=1) <= 1)
        # A regular hexagon of circumradius 10.1 (apothem 10.1 cos 30).
        apothem = 10.1 * math.cos(math.pi / 6)
        normals = np.array([[math.cos(math.pi / 6 + k * math.pi / 3), math.sin(math.pi / 6 + k * math.pi / 3)] for k in range(6)])
        self.assertSameAtoms("prism(n = 6, r = 10.1, h = 12.1);",
                             lambda q: np.all(q[:, :2] @ normals.T <= apothem + 1e-9, axis=1) & (np.abs(q[:, 2]) <= 6.05))
        self.assertSameAtoms("polygon(points = [[0, 0], [12.1, 0], [0, 9.1]], h = 6.1);",
                             lambda q: (q[:, 0] >= 0) & (q[:, 1] >= 0) & (q[:, 0] / 12.1 + q[:, 1] / 9.1 <= 1 + 1e-9) & (np.abs(q[:, 2]) <= 3.05))

    def test_booleans_transforms_and_patterns(self):
        r = np.linalg.norm
        self.assertSameAtoms("difference() { cylinder(r = 12.1, h = 20.1); cylinder(r = 7.1, h = 22); }",
                             lambda q: (r(q[:, :2], axis=1) <= 12.1) & (r(q[:, :2], axis=1) > 7.1) & (np.abs(q[:, 2]) <= 10.05))
        self.assertSameAtoms("intersection() { sphere(r = 12.1); box(size = [30, 30, 10.1]); }",
                             lambda q: (r(q, axis=1) <= 12.1) & (np.abs(q[:, 2]) <= 5.05))
        self.assertSameAtoms("union() { translate([-8, 0, 0]) sphere(r = 6.1); translate([8, 0, 0]) sphere(r = 6.1); }",
                             lambda q: (r(q - [-8, 0, 0], axis=1) <= 6.1) | (r(q - [8, 0, 0], axis=1) <= 6.1))
        c, s_ = math.cos(math.radians(30)), math.sin(math.radians(30))
        turn = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])
        self.assertSameAtoms("rotate([0, 0, 30]) box(size = [30.1, 8.1, 8.1]);",
                             lambda q: np.all(np.abs(q @ turn) <= [15.05, 4.05, 4.05], axis=1))
        self.assertSameAtoms("rotate(a = 30, v = [0, 0, 1]) box(size = [30.1, 8.1, 8.1]);",
                             lambda q: np.all(np.abs(q @ turn) <= [15.05, 4.05, 4.05], axis=1))
        self.assertSameAtoms("scale([2, 1, 1]) sphere(r = 6.1);", lambda q: r(q / [2, 1, 1], axis=1) <= 6.1)
        self.assertSameAtoms("mirror([1, 0, 0]) translate([5, 0, 0]) sphere(r = 4.1);", lambda q: r(q - [-5, 0, 0], axis=1) <= 4.1)
        self.assertSameAtoms("shell(t = 3) sphere(r = 12.1);", lambda q: (r(q, axis=1) <= 12.1) & (r(q, axis=1) >= 9.1))
        self.assertSameAtoms("offset(r = 2) sphere(r = 8.1);", lambda q: r(q, axis=1) <= 10.1)
        self.assertSameAtoms("linear_array(n = 3, step = [14, 0, 0]) sphere(r = 4.1);",
                             lambda q: np.any([r(q - [14 * k, 0, 0], axis=1) <= 4.1 for k in range(3)], axis=0))
        self.assertSameAtoms("grid_array(n = [3, 2, 1], step = [12, 12, 0]) sphere(r = 4.1);",
                             lambda q: np.any([r(q - [12 * i, 12 * j, 0], axis=1) <= 4.1 for i in range(3) for j in range(2)], axis=0))
        centres = [[16 * math.cos(2 * math.pi * k / 4), 16 * math.sin(2 * math.pi * k / 4), 0] for k in range(4)]
        self.assertSameAtoms("polar_array(n = 4) translate([16, 0, 0]) sphere(r = 4.1);",
                             lambda q: np.any([r(q - c_, axis=1) <= 4.1 for c_ in centres], axis=0))
        # Smooth unions fill the neck between two spheres.
        _, plain = self.carved("union() { translate([-7, 0, 0]) sphere(r = 6); translate([7, 0, 0]) sphere(r = 6); }")
        _, smooth = self.carved("smooth_union(k = 4) { translate([-7, 0, 0]) sphere(r = 6); translate([7, 0, 0]) sphere(r = 6); }")
        self.assertGreater(len(smooth), len(plain))
        # Variables and arithmetic.
        self.assertSameAtoms("R = 10; t = R / 5; shell(t = t) sphere(r = R + 0.1);",
                             lambda q: (r(q, axis=1) <= 10.1) & (r(q, axis=1) >= 8.1))

    def test_crystal_facets(self):
        family = [np.array(v) / np.linalg.norm(v) for v in
                  {tuple(s * np.array(p)) for p in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)] for s in (1, -1)}]
        cube = [np.array(v) for v in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]]
        normals = np.array(family + cube)
        offsets = np.array([12.1] * 8 + [14.1] * 6)
        self.assertSameAtoms("facets(hkl = [[1, 1, 1], [1, 0, 0]], d = [12.1, 14.1]);",
                             lambda q: np.all(q @ normals.T <= offsets + 1e-9, axis=1))
        self.assertSameAtoms("octahedron(d = 10.1);", lambda q: np.abs(q).sum(axis=1) / math.sqrt(3) <= 10.1 + 1e-9)
        # Miller indices use the reciprocal lattice: (100) of a hexagonal cell is perpendicular to b and c.
        from ase.build import bulk
        mg = from_ase(bulk("Mg", "hcp", a=3.21, c=5.21) * (4, 4, 3))
        result = shape_step("halfspace(hkl = [1, 0, 0], d = 0);", mg, mode="select", anchor="origin")
        cell = np.array(mg.cell)
        astar = np.linalg.inv(cell)[:, 0]
        expected = [i for i, p in enumerate(positions(mg)) if p @ astar / np.linalg.norm(astar) <= 1e-9]
        self.assertEqual(result.selected_indices, expected)

    def test_regions_assign_elements_without_volume(self):
        r = np.linalg.norm
        result, p = self.assertSameAtoms('sphere(r = 14.1); region("Au") sphere(r = 8.1);', lambda q: r(q, axis=1) <= 14.1)
        s = result.structure
        gold = np.array(symbols(s)) == "Au"
        centre = p.mean(axis=0)
        self.assertEqual(gold.sum(), int((r(self.points, axis=1) <= 8.1).sum()))
        self.assertLessEqual(r(p[gold] - centre, axis=1).max(), 8.1 + 1e-6)
        self.assertIn(f"{gold.sum()} Au in regions", notes_of(result)[0])
        # A region adds no volume: a channel stays open.
        _, open_channel = self.carved('difference() { sphere(r = 12); cylinder(r = 3, h = 30); } region("Au") sphere(r = 8);')
        _, plain = self.carved('difference() { sphere(r = 12); cylinder(r = 3, h = 30); }')
        self.assertEqual(len(open_channel), len(plain))
        # Janus: an unbounded region; "none" regions remove atoms; nested regions win.
        janus = shape_step('sphere(r = 12.1); region("Ag") halfspace(normal = [0, 0, -1], d = 0);', self.crystal, vacuum=0).structure
        z = positions(janus)[:, 2] - positions(janus)[:, 2].mean()
        silver = np.array(symbols(janus)) == "Ag"
        self.assertTrue(np.all(z[silver] >= -1e-6) and np.all(z[~silver] <= 1e-6))
        holed = shape_step('sphere(r = 12.1); region("none") sphere(r = 5.1);', self.crystal).structure
        self.assertEqual(len(holed.atoms), int(((r(self.points, axis=1) <= 12.1) & (r(self.points, axis=1) > 5.1)).sum()))
        nested = shape_step('sphere(r = 12.1); region("Ni") sphere(r = 9.1); region("Au") sphere(r = 5.1);', self.crystal).structure
        counts = {e: symbols(nested).count(e) for e in ("Cu", "Ni", "Au")}
        self.assertEqual(counts["Au"], int((r(self.points, axis=1) <= 5.1).sum()))
        self.assertEqual(counts["Ni"], int(((r(self.points, axis=1) <= 9.1) & (r(self.points, axis=1) > 5.1)).sum()))

    def test_modes_on_the_structure_itself(self):
        from ase.build import bulk
        block = from_ase(bulk("Cu", "fcc", a=3.615, cubic=True) * (6, 6, 6))
        centre = np.array(block.cell).sum(axis=0) / 2
        inside = np.linalg.norm(positions(block) - centre, axis=1) <= 5.1
        void = shape_step("sphere(r = 5.1);", block, mode="delete").structure
        self.assertEqual(len(void.atoms), len(block.atoms) - inside.sum())
        self.assertEqual(void.cell, block.cell)  # still the periodic crystal
        kept = shape_step('sphere(r = 5.1); region("Au") sphere(r = 2.6);', block, mode="keep").structure
        self.assertEqual(len(kept.atoms), inside.sum())
        self.assertEqual(symbols(kept).count("Au"), int((np.linalg.norm(positions(block) - centre, axis=1) <= 2.6).sum()))
        selected = shape_step("sphere(r = 5.1);", block, mode="select").selected_indices
        self.assertEqual(selected, list(np.flatnonzero(inside)))
        moved = shape_step("sphere(r = 5.1);", block, mode="select", anchor="origin", position=list(centre)).selected_indices
        self.assertEqual(moved, selected)

    def test_errors_name_the_line(self):
        cases = {
            "sphere(r = 10;": "line 1",
            "R = 5;\nfoo(1);": "line 2: unknown command 'foo'",
            "sphere(radius = 3);": "no argument 'radius'",
            "translate([1, 2, 3]);": "needs shapes",
            "sphere(r = x);": "unknown variable 'x'",
            'region("Au") sphere(r = 3);': "only regions",
            "box(size = [1, -1, 1]);": "size must be positive",
        }
        for script, message in cases.items():
            with self.assertRaises(PipelineError) as raised:
                shape_step(script, self.crystal)
            self.assertIn(message, str(raised.exception), script)
        with self.assertRaises(PipelineError) as raised:
            shape_step("halfspace(normal = [0, 0, 1], d = 0);", self.crystal)
        self.assertIn("unbounded", str(raised.exception))
        with self.assertRaises(PipelineError) as raised:
            shape_step("sphere(r = 200);", self.crystal, max_atoms=1000)
        self.assertIn("max_atoms", str(raised.exception))

    def test_text_pipelines_and_surface_mesh(self):
        result = run("build-shape script='R = 9; sphere(r = R); region(\"Au\") sphere(r = 5);'", self.crystal)
        self.assertIn("Au", symbols(result.structure))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "ball.obj"
            shape_step("sphere(r = 10);", self.crystal, mesh_file=str(path), mesh_spacing=0.5)
            lines = path.read_text().splitlines()
            vertices = np.array([[float(x) for x in line.split()[1:4]] for line in lines if line.startswith("v ")])
            faces = [line for line in lines if line.startswith("f ")]
            self.assertGreater(len(faces), 1000)
            radii = np.linalg.norm(vertices - vertices.mean(axis=0), axis=1)
            self.assertLess(np.abs(radii - 10).max(), 0.3)
            # A closed surface: every edge is shared by two triangles.
            edges = {}
            for f in faces:
                a, b, c = (int(x) for x in f.split()[1:4])
                for e in ((a, b), (b, c), (c, a)):
                    key = tuple(sorted(e))
                    edges[key] = edges.get(key, 0) + 1
            self.assertTrue(all(n == 2 for n in edges.values()))


if __name__ == "__main__":
    unittest.main()
