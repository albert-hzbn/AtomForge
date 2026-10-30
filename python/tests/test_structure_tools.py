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


if __name__ == "__main__":
    unittest.main()
