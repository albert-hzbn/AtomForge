"""Checks the native analyses against independent implementations.

Each test computes the same quantity another way (scipy's Voronoi
tessellation, matscipy's ring statistics, ASE neighbour lists, brute-force
numpy) on crystals in rotated and triclinic cells, with thermal noise, open
boundaries and edge cases. The reference libraries are optional; tests that
need a missing one are skipped. Run with the native library available
(ATOMFORGE_SCIENCE_LIBRARY), e.g. python -m pytest python/tests/test_reference_validation.py
"""

import math
import os
from pathlib import Path
import sys
import unittest

if sys.platform == "win32" and os.path.isdir(r"C:\msys64\ucrt64\bin"):
    os.add_dll_directory(r"C:\msys64\ucrt64\bin")  # runtime of a development build
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

from atomforge import native


def available(module):
    try:
        __import__(module)
        return True
    except ImportError:
        return False


def lattice(cell, basis, reps):
    cell = np.asarray(cell, float)
    frac = np.array([np.add(b, (i, j, k)) for i in range(reps[0]) for j in range(reps[1]) for k in range(reps[2]) for b in basis])
    return frac @ cell, cell * np.asarray(reps)[:, None]


def rotation(axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * k @ k


FCC = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
BCC = [(0, 0, 0), (.5, .5, .5)]
HCP_BASIS = [(0, 0, 0), (1 / 3, 2 / 3, .5)]
DIAMOND = FCC + [tuple(np.add(b, .25)) for b in FCC] if np is not None else []


def hexagonal(a, ca):
    return [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, a * ca]]


def run(tool, **request):
    return native.run(tool, request).result


def periodic(positions, cell):
    return {"positions": np.asarray(positions).tolist(), "cell": np.asarray(cell).tolist(), "pbc": [True, True, True]}


def minimum_image(d, cell):
    f = d @ np.linalg.inv(cell)
    f -= np.round(f)
    return f @ cell


@unittest.skipIf(np is None, "numpy is required")
class StructureIdentification(unittest.TestCase):
    """Crystals in arbitrary orientation, in triclinic (primitive) cells and with noise."""

    def classify(self, tool, positions, cell, **extra):
        return np.array(run(tool, **periodic(positions, cell), **extra)["structure_type"])

    def test_fcc_primitive_triclinic_cell(self):
        a = 3.6
        primitive = np.array([[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]])
        pos, cell = lattice(primitive, [(0, 0, 0)], (4, 4, 4))
        r = rotation([1, 1, 3], 1.1)
        for tool in ("ptm", "cna"):
            types = self.classify(tool, pos @ r.T, cell @ r.T)
            self.assertTrue((types == 1).all(), tool)
        result = run("ptm", **periodic(pos @ r.T, cell @ r.T))
        self.assertTrue(np.allclose(result["lattice_constant_A"], a))

    def test_orientation_is_the_applied_rotation_up_to_symmetry(self):
        pos, cell = lattice(np.eye(3) * 3.6, FCC, (3, 3, 3))
        rng = np.random.default_rng(1)
        for _ in range(5):
            axis, angle = rng.normal(size=3), rng.uniform(0, math.pi)
            r = rotation(axis, angle)
            q = np.array(run("ptm", **periodic(pos @ r.T, cell @ r.T))["orientation"][0])
            # Rotation of the quaternion, compared with r up to the 24 cubic symmetries.
            w, x, y, z = q
            rq = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
                           [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
                           [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])
            delta = rq.T @ r
            # delta must be a cubic symmetry: a signed permutation matrix with det +1.
            self.assertTrue(np.allclose(np.abs(delta), np.round(np.abs(delta)), atol=1e-6), delta)
            self.assertAlmostEqual(np.linalg.det(delta), 1, places=6)
            self.assertLessEqual(2 * math.degrees(math.acos(min(1, q[0]))), 62.8 + 1e-6)  # fundamental zone of 432

    def test_hcp_and_bcc_with_noise(self):
        rng = np.random.default_rng(2)
        hcp, hcp_cell = lattice(hexagonal(2.95, 1.633), HCP_BASIS, (5, 5, 4))
        bcc, bcc_cell = lattice(np.eye(3) * 2.87, BCC, (5, 5, 5))
        for pos, cell, expected, tool_types in ((hcp, hcp_cell, 2, ("ptm", "cna")), (bcc, bcc_cell, 3, ("ptm", "cna", "structure-type"))):
            noisy = pos + rng.normal(0, 0.04, pos.shape)
            for tool in tool_types:
                types = self.classify(tool, noisy, cell)
                self.assertGreater((types == expected).mean(), 0.95, tool)

    def test_vacancy_neighbours_are_not_fcc_for_cna_and_ptm(self):
        pos, cell = lattice(np.eye(3) * 3.6, FCC, (4, 4, 4))
        pos = pos[1:]  # one vacancy at the origin
        for tool in ("cna", "ptm"):
            types = self.classify(tool, pos, cell)
            self.assertEqual((types != 1).sum(), 12, tool)  # its twelve neighbours

    def test_open_cluster_surface_is_other(self):
        pos, _ = lattice(np.eye(3) * 3.6, FCC, (5, 5, 5))
        pos -= pos.mean(0)
        ball = pos[np.linalg.norm(pos, axis=1) < 8]
        for tool in ("ptm", "cna"):
            types = np.array(run(tool, positions=ball.tolist())["structure_type"])
            inner = np.linalg.norm(ball, axis=1) < 4
            self.assertTrue((types[inner] == 1).all(), tool)
            self.assertTrue((types[np.linalg.norm(ball, axis=1) > 7] == 0).all(), tool)

    def test_edge_cases(self):
        for tool in ("ptm", "cna", "voronoi", "vorotop", "chill-plus", "local-entropy", "surface-mesh", "bond-orders", "fcc-planar-faults"):
            extra = {"types": [29]} if tool == "bond-orders" else {}
            extra.update({"local_density": True} if tool == "local-entropy" else {})
            result = run(tool, positions=[[0.0, 0.0, 0.0]], **extra)
            self.assertIsInstance(result, dict, tool)
        with self.assertRaises(ValueError):
            run("ptm", positions=[[0, 0, 0], [0, 0, 0]], cell=np.eye(3).tolist(), pbc=[True] * 3)
        with self.assertRaises(ValueError):
            run("voronoi", positions=[[0, 0, 0]], pbc=[True, True, True])  # periodic without a cell


@unittest.skipIf(np is None or not available("scipy"), "scipy is required")
class VoronoiAgainstScipy(unittest.TestCase):
    def test_volumes_and_neighbours_in_a_triclinic_cell(self):
        from scipy.spatial import Voronoi
        rng = np.random.default_rng(3)
        cell = np.array([[9.0, 0, 0], [2.5, 8.0, 0], [1.5, -2.0, 10.0]])
        frac = rng.uniform(0, 1, (60, 3))
        pos = frac @ cell
        result = run("voronoi", **periodic(pos, cell))
        # Reference: the central copy of a 3x3x3 tiling.
        shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
        tiled = np.vstack([pos + s @ cell for s in shifts])
        centre = np.flatnonzero((shifts == 0).all(1))[0] * len(pos)
        vor = Voronoi(tiled)
        from scipy.spatial import ConvexHull
        volumes, faces = [], []
        for i in range(len(pos)):
            region = vor.regions[vor.point_region[centre + i]]
            volumes.append(ConvexHull(vor.vertices[region]).volume)
            faces.append(sum(1 for a, b in vor.ridge_points if a == centre + i or b == centre + i))
        self.assertTrue(np.allclose(result["atomic_volume_A3"], volumes, rtol=1e-7), np.max(np.abs(np.subtract(result["atomic_volume_A3"], volumes))))
        self.assertTrue(np.array_equal(result["coordination"], faces))
        self.assertAlmostEqual(result["total_volume_A3"], abs(np.linalg.det(cell)), places=6)

    def test_open_boundaries_mark_surface_cells_open(self):
        rng = np.random.default_rng(4)
        pos = rng.uniform(0, 10, (40, 3))
        result = run("voronoi", positions=pos.tolist())
        from scipy.spatial import Voronoi
        vor = Voronoi(pos)
        unbounded = [(-1 in vor.regions[vor.point_region[i]]) for i in range(len(pos))]
        self.assertEqual(result["open_cell"], unbounded)
        for v, is_open in zip(result["atomic_volume_A3"], unbounded):
            self.assertEqual(v is None, is_open)


@unittest.skipIf(np is None or not available("matscipy") or not available("ase"), "matscipy and ase are required")
class RingsAgainstMatscipy(unittest.TestCase):
    def check(self, atoms, cutoff, max_size):
        """Ring counts up to max_size (matscipy's maxlength excludes rings of that length)."""
        from matscipy.rings import ring_statistics
        reference = np.zeros(max_size + 1, int)
        counts = np.asarray(ring_statistics(atoms, cutoff, maxlength=max_size + 2))
        reference[:min(len(counts), max_size + 1)] = counts[:max_size + 1]
        result = run("find-rings", positions=atoms.positions.tolist(), cell=atoms.cell.array.tolist(),
                     pbc=[bool(b) for b in atoms.pbc], cutoff_A=cutoff, max_size=max_size)
        ours = np.zeros(max_size + 1, int)
        for size, count in zip(result["ring_sizes"], result["ring_counts"]):
            ours[size] = count
        self.assertTrue(np.array_equal(ours, reference), (ours, reference))

    def test_diamond_and_disordered_networks(self):
        from ase.build import bulk
        si = bulk("Si", "diamond", a=5.43, cubic=True).repeat((3, 3, 3))
        self.check(si, 2.6, 8)
        rng = np.random.default_rng(5)
        for noise, cutoff in ((0.15, 2.75), (0.25, 2.8)):
            disordered = si.copy()
            disordered.positions += rng.normal(0, noise, si.positions.shape)
            self.check(disordered, cutoff, 10)

    def test_close_packed_metal(self):
        from ase.build import bulk
        cu = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((4, 4, 4))
        cu.positions += np.random.default_rng(1).normal(0, 0.05, cu.positions.shape)
        self.check(cu, 2.6, 8)


@unittest.skipIf(np is None or not available("ase"), "ase is required")
class BondsAgainstAse(unittest.TestCase):
    def test_bond_counts_and_angles_of_a_disordered_alloy(self):
        from ase.build import bulk
        from ase.neighborlist import neighbor_list
        rng = np.random.default_rng(6)
        atoms = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((3, 3, 3))
        atoms.numbers[rng.random(len(atoms)) < 0.4] = 28
        atoms.positions += rng.normal(0, 0.08, atoms.positions.shape)
        cutoff = 3.0
        i, j, d, vectors = neighbor_list("ijdD", atoms, cutoff)
        result = run("bond-distributions", positions=atoms.positions.tolist(), cell=atoms.cell.array.tolist(), pbc=[True] * 3,
                     cutoff_A=cutoff, types=atoms.numbers.tolist(), bins=60)
        self.assertEqual(result["bonds"], len(i) // 2)
        pairs = {}
        for a, b in zip(i, j):
            if a < b:
                key = f"{min(atoms.numbers[a], atoms.numbers[b])}-{max(atoms.numbers[a], atoms.numbers[b])}"
                pairs[key] = pairs.get(key, 0) + 1
        ours = {k: sum(v) for k, v in result["bond_length_counts_by_pair"].items()}
        self.assertEqual(ours, pairs)
        angles = sum(len(np.flatnonzero(i == a)) * (len(np.flatnonzero(i == a)) - 1) // 2 for a in range(len(atoms)))
        self.assertEqual(result["angles"], angles)
        self.assertAlmostEqual(result["mean_bond_length_A"], d.mean(), places=9)


@unittest.skipIf(np is None, "numpy is required")
class PropertyAndTrajectoryTools(unittest.TestCase):
    def test_local_entropy_matches_a_direct_integration(self):
        rng = np.random.default_rng(7)
        pos, cell = lattice(np.eye(3) * 3.6, FCC, (3, 3, 3))
        pos = pos + rng.normal(0, 0.1, pos.shape)
        rc, sigma = 5.0, 0.2
        result = run("local-entropy", **periodic(pos, cell), cutoff_A=rc, sigma_A=sigma, average=False)
        rho = len(pos) / abs(np.linalg.det(cell))
        r = np.linspace(0, rc, 4001)[1:]
        dr = r[1] - r[0]
        reference = []
        for i in (0, 5, 17):
            # Every image within the cutoff plus the Gaussian tail.
            dist = []
            for s in np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]):
                d = np.linalg.norm(pos + s @ cell - pos[i], axis=1)
                dist.extend(d[(d > 1e-9) & (d < rc + 3 * sigma)])
            dist = np.array(dist)
            g = np.exp(-(r[:, None] - dist[None, :]) ** 2 / (2 * sigma * sigma)).sum(1) / math.sqrt(2 * math.pi * sigma * sigma) / (4 * math.pi * rho * r * r)
            integrand = np.where(g > 0, g * np.log(np.where(g > 0, g, 1)), 0) - g + 1
            reference.append(-2 * math.pi * rho * np.sum(integrand * r * r) * dr)
        ours = [result["entropy_kB"][i] for i in (0, 5, 17)]
        self.assertTrue(np.allclose(ours, reference, rtol=2e-3, atol=2e-3), (ours, reference))

    def test_spatial_correlation_matches_pair_sums(self):
        rng = np.random.default_rng(8)
        cell = np.eye(3) * 12.0
        pos = rng.uniform(0, 12, (150, 3))
        a, b = rng.normal(size=150), rng.normal(size=150)
        b = 0.6 * a + 0.8 * b
        result = run("spatial-correlation", positions=pos.tolist(), cell=cell.tolist(), values_a=a.tolist(), values_b=b.tolist(),
                     grid_spacing_A=0.25, bins=6, r_max_A=6.0)
        # Brute force with the grid's own nearest-point positions.
        n = 48
        grid = np.round(pos / 12.0 * n) % n * (12.0 / n)
        d = grid[:, None, :] - grid[None, :, :]
        d -= np.round(d / 12.0) * 12.0
        dist = np.linalg.norm(d, axis=2)
        edges = np.linspace(0, 6.0, 7)
        for k in range(6):
            mask = (dist >= edges[k]) & (dist < edges[k + 1])
            if mask.sum() == 0:
                continue
            expected = (a[:, None] * b[None, :])[mask].mean()
            self.assertAlmostEqual(result["correlation"][k], expected, places=9)

    def test_spatial_binning_matches_numpy(self):
        rng = np.random.default_rng(9)
        cell = np.array([[10.0, 0, 0], [3.0, 9.0, 0], [0, 0, 8.0]])
        frac = rng.uniform(-0.5, 1.5, (500, 3))
        values = rng.normal(size=500)
        result = run("spatial-binning", positions=(frac @ cell).tolist(), cell=cell.tolist(), values=values.tolist(),
                     bins=[4, 3, 2], operation="sum")
        wrapped = frac - np.floor(frac)
        reference, _ = np.histogramdd(wrapped, bins=(4, 3, 2), range=((0, 1), (0, 1), (0, 1)), weights=values)
        self.assertTrue(np.allclose(np.array(result["values"]).reshape(4, 3, 2), reference))
        counts, _ = np.histogramdd(wrapped, bins=(4, 3, 2), range=((0, 1), (0, 1), (0, 1)))
        mean = run("spatial-binning", positions=(frac @ cell).tolist(), cell=cell.tolist(), values=values.tolist(), bins=[4, 3, 2], operation="mean")
        expected = np.where(counts > 0, reference / np.maximum(counts, 1), np.nan)
        self.assertTrue(np.allclose(np.array(mean["values"], dtype=float).reshape(4, 3, 2), expected, equal_nan=True))

    def test_displacements_in_a_triclinic_cell(self):
        rng = np.random.default_rng(10)
        cell = np.array([[8.0, 0, 0], [4.0, 7.0, 0], [2.0, 2.0, 9.0]])
        ref = rng.uniform(0, 1, (30, 3)) @ cell
        step = rng.normal(0, 0.8, (30, 3))
        cur = ref + step + rng.integers(-1, 2, (30, 3)) @ cell  # jumps across the boundary
        result = run("displacements", reference=ref.tolist(), current=cur.tolist(), reference_cell=cell.tolist(), pbc=[True] * 3)
        self.assertTrue(np.allclose(result["displacement_A"], step, atol=1e-9))

    def test_elastic_strain_recovers_a_random_deformation(self):
        rng = np.random.default_rng(11)
        pos, cell = lattice(np.eye(3) * 3.6, FCC, (3, 3, 3))
        sym = rng.normal(0, 0.01, (3, 3))
        f = np.eye(3) + (sym + sym.T) / 2
        r = rotation(rng.normal(size=3), 0.7)
        deformation = r @ f
        result = run("elastic-strain", **periodic(pos @ deformation.T, cell @ deformation.T), lattice_constant_A=3.6)
        expected = 0.5 * (deformation @ deformation.T - np.eye(3))
        voigt = [expected[0, 0], expected[1, 1], expected[2, 2], expected[1, 2], expected[0, 2], expected[0, 1]]
        for row in result["strain_voigt"]:
            self.assertTrue(np.allclose(row, voigt, atol=1e-10))

    def test_regions_against_numpy(self):
        rng = np.random.default_rng(12)
        pts = rng.uniform(-2, 12, (2000, 3))
        cell = np.eye(3) * 10
        commands = ("region a block 1 4 INF 6 -INF 5 units box; region b sphere 5 5 5 2.5 side out units box; "
                    "region c cylinder y 4 6 3 2 EDGE units box; region d cone x 5 5 1 3 2 8 units box; "
                    "region e plane 5 5 5 1 1 0 units box; region f prism 1 6 2 7 3 8 2 1 0.5 units box")
        result = run("lammps-regions", positions=pts.tolist(), cell=cell.tolist(), regions=commands, export_mesh=False)
        x, y, z = pts.T
        expected = {
            "a": (x >= 1) & (x <= 4) & (y <= 6) & (z <= 5),
            "b": np.linalg.norm(pts - 5, axis=1) > 2.5 - 1e-12,
            "c": (np.hypot(x - 4, z - 6) <= 3) & (y >= 2) & (y <= 10),
            "d": (x >= 2) & (x <= 8) & (np.hypot(y - 5, z - 5) <= 1 + (x - 2) / 6 * 2),
            "e": (x - 5) + (y - 5) >= 0,
        }
        prism = np.array([[5.0, 0, 0], [2.0, 5.0, 0], [1.0, 0.5, 5.0]])
        f = (pts - [1, 2, 3]) @ np.linalg.inv(prism)
        expected["f"] = ((f >= 0) & (f <= 1)).all(1)
        for name, mask in expected.items():
            ours = np.array(result["inside"][name])
            self.assertLessEqual(int((ours != mask).sum()), 0, name)

    def test_trajectory_unwrap_and_average_against_numpy(self):
        import tempfile
        rng = np.random.default_rng(13)
        cell = np.array([[6.0, 0, 0], [2.0, 6.0, 0], [0, 0, 6.0]])
        steps = rng.normal(0, 0.3, (40, 5, 3))
        path = np.cumsum(steps, axis=0) + rng.uniform(0, 6, (1, 5, 3))
        frac = path @ np.linalg.inv(cell)
        wrapped = (frac - np.floor(frac)) @ cell
        directory = Path(tempfile.mkdtemp())
        lines = []
        for frame in wrapped:
            lines.append(f"{len(frame)}\nLattice=\"{' '.join(map(str, cell.flatten()))}\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"")
            lines += [f"Cu {x:.12f} {y:.12f} {z:.12f}" for x, y, z in frame]
        trajectory = directory / "walk.extxyz"
        trajectory.write_text("\n".join(lines) + "\n")
        out = native.run("unwrap-trajectory", trajectory_file=str(trajectory))
        unwrapped = np.array([[[a.x, a.y, a.z] for a in f.atoms] for f in out.frames])
        self.assertTrue(np.allclose(unwrapped - unwrapped[0], path - path[0], atol=1e-8))
        average = native.run("trajectory-average", trajectory_file=str(trajectory))
        mean = np.array([[a.x, a.y, a.z] for a in average.frames[0].atoms])
        self.assertTrue(np.allclose(mean - unwrapped[0], (path - path[0]).mean(0), atol=1e-8))


def cubic_misorientation(r1, r2):
    """Smallest rotation angle (degrees) between two cubic orientations."""
    delta = r1.T @ r2
    best = 0.0
    import itertools
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            s = np.zeros((3, 3))
            for row, (col, sign) in enumerate(zip(perm, signs)):
                s[row, col] = sign
            if np.linalg.det(s) < 0:
                continue
            best = max(best, (np.trace(delta @ s) - 1) / 2)
    return math.degrees(math.acos(min(1.0, best)))


@unittest.skipIf(np is None, "numpy is required")
class GrainsFaultsAndWater(unittest.TestCase):
    def test_tricrystal_misorientations(self):
        rng = np.random.default_rng(20)
        block, _ = lattice(np.eye(3) * 3.6, FCC, (14, 14, 14))
        block -= block.mean(0)
        rotations = [rotation(rng.normal(size=3), rng.uniform(0.3, 1.0)) for _ in range(3)]
        # Three sectors around the z axis, each a differently rotated crystal.
        atoms, owner = [], []
        for g, r in enumerate(rotations):
            pts = block @ r.T
            angle = np.arctan2(pts[:, 1], pts[:, 0]) % (2 * math.pi)
            keep = (angle >= g * 2 * math.pi / 3) & (angle < (g + 1) * 2 * math.pi / 3) & (np.linalg.norm(pts[:, :2], axis=1) < 20) & (np.abs(pts[:, 2]) < 8)
            atoms.append(pts[keep])
            owner += [g] * int(keep.sum())
        atoms = np.vstack(atoms) + rng.normal(0, 0.03, (sum(len(a) for a in atoms), 3))
        result = run("grain-segmentation", positions=atoms.tolist(), min_grain_size=100)
        self.assertEqual(result["grain_count"], 3)
        expected = sorted(cubic_misorientation(rotations[i], rotations[j]) for i, j in ((0, 1), (0, 2), (1, 2)))
        found = sorted(b["misorientation_deg"] for b in result["boundaries"])
        self.assertEqual(len(found), 3)
        self.assertTrue(np.allclose(found, expected, atol=1.5), (found, expected))
        # Grain ids follow the sectors (away from boundaries and the outer surface).
        ids = np.array(result["grain_id"])
        owner = np.array(owner)
        for g in range(3):
            interior = ids[(owner == g)]
            values, counts = np.unique(interior, return_counts=True)
            self.assertGreater(counts.max() / counts.sum(), 0.9)

    def test_planar_faults_with_noise_and_extrinsic_fault(self):
        def stack(sequence, a=2.55, n=6):
            s3, d = math.sqrt(3), a * math.sqrt(2 / 3)
            shift = {"A": 0, "B": 1, "C": 2}
            pts = [(i * a + j * a / 2 + shift[c] * a / 2, j * a * s3 / 2 + shift[c] * a / (2 * s3), l * d)
                   for l, c in enumerate(sequence) for i in range(n) for j in range(n)]
            cell = [[n * a, 0, 0], [n * a / 2, n * a * s3 / 2, 0], [0, 0, len(sequence) * d]]
            return np.array(pts), np.array(cell)
        rng = np.random.default_rng(21)
        pos, cell = stack("ABCABCABABCABCABC")  # intrinsic fault
        result = run("fcc-planar-faults", **periodic(pos + rng.normal(0, 0.05, pos.shape), cell))
        self.assertEqual(result["counts"]["intrinsic stacking fault"], 72)
        # Extrinsic fault ...ABC|B|ABC...: two single hcp planes one layer apart, each twin-like.
        pos, cell = stack("ABCABCBABCABCABC")
        result = run("fcc-planar-faults", **periodic(pos, cell))
        self.assertEqual(result["counts"]["coherent twin boundary"] + result["counts"]["other hcp"], 72)

    def test_chill_plus_on_water_with_hydrogens_through_the_pipeline(self):
        from atomforge.pipeline import Pipeline
        import atomforge as af
        # Ice Ic oxygens with two hydrogens each (positions need not be exact for Chill+).
        oo = 2.76
        a = oo * 4 / math.sqrt(3)
        oxygens, cell = lattice(np.eye(3) * a, DIAMOND, (2, 2, 2))
        s = af.Structure()
        s.cell = cell.tolist()
        for o in oxygens:
            s.add_atom("O", *o)
            s.add_atom("H", *(o + [0.6, 0.6, 0.4]))
            s.add_atom("H", *(o + [-0.6, 0.6, -0.4]))
        result = Pipeline.parse("compute-property chill-plus").run(s)
        prop = np.array(result.property, dtype=float)
        symbols = np.array([at.symbol for at in result.structure.atoms])
        self.assertTrue(np.all(prop[symbols == "O"] == 2))  # cubic ice
        self.assertTrue(np.all(np.isnan(prop[symbols == "H"])))


@unittest.skipIf(np is None, "numpy is required")
class SurfacesBondsRegionsAndTrajectories(unittest.TestCase):
    def test_surface_area_of_a_dense_sphere_converges(self):
        # Atoms on a fine spherical shell (radius 10) around a filled interior:
        # the probe surface is that sphere, up to the grid and the gaps between atoms.
        n = 6000
        k = np.arange(n) + 0.5
        polar, azimuth = np.arccos(1 - 2 * k / n), math.pi * (1 + 5 ** 0.5) * k
        shell = 10 * np.column_stack([np.sin(polar) * np.cos(azimuth), np.sin(polar) * np.sin(azimuth), np.cos(polar)])
        grid = np.mgrid[-9:9.01:1.5, -9:9.01:1.5, -9:9.01:1.5].reshape(3, -1).T
        ball = np.vstack([shell, grid[np.linalg.norm(grid, axis=1) < 9]])
        result = run("surface-mesh", positions=ball.tolist(), probe_radius_A=3.0, grid_spacing_A=0.25, export_mesh=False)
        self.assertAlmostEqual(result["surface_area_A2"] / (4 * math.pi * 100), 1.0, delta=0.03)
        self.assertAlmostEqual(result["solid_volume_A3"] / (4 / 3 * math.pi * 1000), 1.0, delta=0.03)
        self.assertEqual(result["solid_regions"], 1)
        with self.assertRaises(ValueError):
            run("surface-mesh", positions=ball[:10].tolist(), cell=(np.eye(3) * 30).tolist(), pbc=[True, True, False])

    def test_bond_orders_of_molecules(self):
        from ase.build import molecule
        for name, expected in (("C2H6", 1.0), ("C2H4", 2.0), ("C2H2", 3.0), ("C6H6", 1.5)):
            mol = molecule(name)
            result = run("bond-orders", positions=mol.positions.tolist(), types=mol.numbers.tolist())
            cc = [b["rounded_order"] for b in result["bonds"] if all(mol.numbers[k] == 6 for k in b["atoms"])]
            self.assertTrue(cc and all(o == expected for o in cc), (name, cc))
            ch = [b["rounded_order"] for b in result["bonds"] if sorted(mol.numbers[k] for k in b["atoms"]) == [1, 6]]
            self.assertTrue(all(o == 1 for o in ch), (name, ch))

    def test_region_units_and_errors(self):
        pts = np.array([[1.0, 1.0, 1.0], [5.0, 5.0, 5.0], [9.0, 9.0, 9.0]])
        result = run("lammps-regions", positions=pts.tolist(), regions="region r block 0 1 0 1 0 1", lattice_scale=4.0, export_mesh=False)
        self.assertEqual(result["inside"]["r"], [True, False, False])  # LAMMPS default units lattice
        for bad in ("region x union 2 a b", "region x block 0 1", "fix 1 all nve", "region x block 0 1 0 1 0 1 move v_x NULL NULL"):
            with self.assertRaises(ValueError):
                run("lammps-regions", positions=pts.tolist(), regions=bad)

    def test_time_series_and_frame_difference(self):
        import tempfile
        rng = np.random.default_rng(23)
        pos, cell = lattice(np.eye(3) * 3.6, FCC, (3, 3, 3))
        frames = [pos + rng.normal(0, 0.02 * k, pos.shape) for k in range(6)]
        directory = Path(tempfile.mkdtemp())
        lines = []
        for frame in frames:
            lines.append(f"{len(frame)}\nLattice=\"{' '.join(map(str, cell.flatten()))}\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"")
            lines += [f"Cu {x:.12f} {y:.12f} {z:.12f}" for x, y, z in frame]
        path = directory / "heat.extxyz"
        path.write_text("\n".join(lines) + "\n")
        series = run("time-series", trajectory_file=str(path), analysis="centrosymmetry", cutoff_A=3.0)
        csp = series["mean_property"]
        self.assertAlmostEqual(csp[0], 0.0, places=12)
        self.assertTrue(all(b > a for a, b in zip(csp, csp[1:])), csp)  # more disorder, larger CSP
        diff = run("frame-difference", trajectory_file=str(path), reference_frame=0, frame=-1)
        step = frames[-1] - frames[0]
        self.assertTrue(np.allclose(diff["difference"], np.linalg.norm(step, axis=1), atol=1e-9))


@unittest.skipIf(np is None, "numpy is required")
class PipelineEdgeCases(unittest.TestCase):
    def setUp(self):
        import atomforge as af
        self.af = af

    def structure(self, positions, cell=None, symbol="Cu"):
        s = self.af.Structure()
        for p in positions:
            s.add_atom(symbol, *p)
        if cell is not None:
            s.cell = np.asarray(cell).tolist()
        return s

    def test_overlaps_in_a_triclinic_periodic_cell(self):
        from atomforge.pipeline import Pipeline
        rng = np.random.default_rng(24)
        cell = np.array([[8.0, 0, 0], [5.0, 7.0, 0], [2.0, 3.0, 9.0]])
        pts = rng.uniform(0, 1, (80, 3)) @ cell
        result = Pipeline.parse("select-overlapping 1.6 atoms=both").run(self.structure(pts, cell))
        expected = set()
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                if np.linalg.norm(minimum_image(pts[j] - pts[i], cell)) < 1.6:
                    expected |= {i, j}
        self.assertEqual(result.selected_count, len(expected))

    def test_formula_frozen_and_statistics(self):
        from atomforge.pipeline import Pipeline
        pos, cell = lattice(np.eye(3) * 3.6, FCC, (2, 2, 2))
        result = Pipeline.parse('compute-expression "z / 0" name=bad | freeze-property | compute-expression "x" | '
                                'select-expression "frozen != frozen" | reduce-property').run(self.structure(pos, cell))
        # 0/0 is NaN (atoms at z = 0), the rest infinite; NaN != NaN selects exactly the z = 0 atoms.
        self.assertEqual(result.selected_count, int((pos[:, 2] == 0).sum()))
        notes = result.stages[-1]["notes"][0]
        self.assertIn("n=32", notes)

    def test_python_script_failures_are_reported(self):
        import sys as system
        import tempfile
        from atomforge.pipeline import Pipeline, PipelineError
        script = Path(tempfile.mkdtemp()) / "broken.py"
        script.write_text("raise RuntimeError('deliberate failure')\n", encoding="utf-8")
        pos, cell = lattice(np.eye(3) * 3.6, FCC, (1, 1, 1))
        steps = Pipeline()
        steps.add("python-script", script=str(script), python=system.executable)
        with self.assertRaises(PipelineError) as raised:
            steps.run(self.structure(pos, cell))
        self.assertIn("deliberate failure", str(raised.exception))

    def test_set_pbc_controls_wrapping_and_neighbours(self):
        from atomforge.pipeline import Pipeline
        pos = np.array([[0.2, 0.2, 0.2], [0.2, 0.2, 9.9], [0.2, 0.2, 10.5]])
        s = self.structure(pos, np.eye(3) * 10)
        wrapped = Pipeline.parse("set-pbc c=false | wrap").run(s)
        self.assertAlmostEqual(wrapped.structure.atoms[2].z, 10.5)
        near = Pipeline.parse("select-sphere 0.2 0.2 0.2 0.5").run(s)
        self.assertEqual(near.selected_count, 3)  # 9.9 and 10.5 are 0.3 away across the periodic c boundary
        open_c = Pipeline.parse("set-pbc c=false | select-sphere 0.2 0.2 0.2 0.5").run(s)
        self.assertEqual(open_c.selected_count, 1)

    def test_every_step_runs_on_an_empty_structure_or_reports_why(self):
        from atomforge.pipeline import Pipeline, PipelineError, modifiers
        empty = self.af.Structure()
        empty.cell = (np.eye(3) * 5).tolist()
        skip = {"merge", "displacements", "python-script", "insert-interstitials", "relax", "nvt-dynamics", "npt-dynamics", "standardize-cell"}
        for m in modifiers():
            if m["id"] in skip or m["id"].startswith("build-"):
                continue
            try:
                Pipeline.parse(m["id"]).run(empty)
            except PipelineError as error:
                self.assertNotIn("std::", str(error), m["id"])  # a readable message, not a crash


if __name__ == "__main__":
    unittest.main(verbosity=2)
