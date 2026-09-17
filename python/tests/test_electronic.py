"""Native/Python integration and independent volumetric format fixtures."""

import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from atomforge.electronic import Grid, load_volume, miller_section, model_density, scattering_factor


class ElectronicTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def grid(self, n=4, periodic=True, value=2):
        return Grid((n, n, n), ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
                    [value] * n**3, periodic=periodic, unit="e/A^3").with_sites([(1, 0, 0, 0)])

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_charge_transfer_and_boolean_masks(self):
        total, fragment = self.grid(value=3), self.grid(value=1)
        difference = total.density_difference(fragment, fragment, weights=[1, 2])
        self.assertEqual(difference.charge_summary(), {"accumulation": 0, "depletion": 0, "net": 0})
        signed = Grid((2, 2, 2), ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                      [-2, 2] * 4, unit="e/A^3")
        positive, negative = signed.threshold_mask(0, 2), signed.threshold_mask(-2, 0)
        self.assertEqual(positive.invert_mask().values, negative.values)
        self.assertEqual(positive.boolean(negative, "union").values, [1] * 8)
        self.assertEqual(positive.boolean(negative, "intersection").values, [0] * 8)
        self.assertEqual(positive.boolean(negative, "xor").values, [1] * 8)
        self.assertAlmostEqual(signed.apply_mask(positive).integrate(), 1)
        accumulation, depletion = signed.split_density()
        self.assertEqual((accumulation - depletion).values, signed.values)
        self.assertEqual(signed.charge_summary(), {"accumulation": 1, "depletion": 1, "net": 0})
        self.assertAlmostEqual(signed.cumulative_charge(0)[-1][1], signed.integrate())
        self.assertAlmostEqual(total.cumulative_charge(2)[-1][1], total.integrate())
        with self.assertRaises(ValueError):
            signed.apply_mask(total)
        with self.assertRaises(ValueError):
            signed.boolean(negative)
        with self.assertRaises(ValueError):
            total.density_difference(fragment, weights=[])
        with self.assertRaises(ValueError):
            total.threshold_mask(2, 1)

    def test_vasp_normalization_order_spin_augmentation(self):
        text = "test\n-8\n1 0 0\n0 1 0\n0 0 1\nH\n1\nSelective dynamics\nDirect\n.5 .25 0 T T T\n\n2 2 2\n8 16 24 32 40 48 56 64\naugmentation occupancies 1 3\n1 2 3\n\n2 2 2\n0 8 0 8 0 8 0 8\n"
        fields = load_volume(self.write("CHGCAR", text)).fields
        self.assertEqual(len(fields), 2)
        self.assertEqual(fields[0].shape, (2, 2, 2))
        self.assertEqual(fields[0].values, list(range(1, 9)))
        self.assertAlmostEqual(fields[0].integrate(), 36)
        self.assertEqual(fields[0].sites, [(1, 1.0, .5, 0.0)])
        self.assertEqual(fields[1].name, "magnetization")
        self.assertEqual(fields[1].values, [0, 1] * 4)
        with self.assertRaises(ValueError):
            load_volume(self.write("CHGCAR-bad", text.replace("1 2 3\n", "1 2\n")))
        potential = load_volume(self.write("LOCPOT", text.split("augmentation")[0])).fields[0]
        self.assertEqual(potential.unit, "eV")
        self.assertEqual(potential.values[0], 8)

    def test_vasp_component_scaling_and_four_channels(self):
        header = "spin\n2 3 4\n1 0 0\n0 1 0\n0 0 1\nH\n1\nCartesian\n.5 .5 .5\n"
        text = header + "".join("\n2 2 2\n" + " ".join([str(c)] * 8) + "\n" for c in (24, 48, 72, 96))
        fields = load_volume(self.write("CHGCAR", text)).fields
        self.assertEqual([g.name for g in fields], ["total", "mx", "my", "mz"])
        self.assertEqual(fields[0].sites, [(1, 1, 1.5, 2)])
        self.assertEqual([g.values[0] for g in fields], [1, 2, 3, 4])

    def test_cube_order_origin_multiple_orbitals_and_units(self):
        # z varies fastest; alternating orbital channels, Fortran exponents.
        text = "cube\nfixture\n-1 1 2 3\n2 .5 0 0\n2 0 .5 0\n2 0 0 .5\n1 1 1 2 3\n2 5\n7\n" + " ".join("{}D0 {}D0".format(i, 10+i) for i in range(8))
        fields = load_volume(self.write("orbitals.cube", text)).fields
        self.assertEqual([g.name for g in fields], ["orbital 5", "orbital 7"])
        self.assertEqual(fields[0].values, [0, 4, 2, 6, 1, 5, 3, 7])
        self.assertEqual(fields[1].values, [10, 14, 12, 16, 11, 15, 13, 17])
        self.assertAlmostEqual(fields[0].origin[0], .529177210903)
        density = load_volume(self.write("density.cube", text), quantity="density").fields[0]
        self.assertAlmostEqual(density.values[1], 4/.529177210903**3)
        with self.assertRaises(ValueError):
            load_volume(self.write("extra.cube", text + " 99"))

    def test_cube_nval_and_explicit_angstrom(self):
        text = "cube\nfixture\n1 1 2 3 2\n-2 1 0 0\n-2 0 1 0\n-2 0 0 1\n1 1 1 2 3\n" + "1 2 " * 8
        fields = load_volume(self.write("two.cube", text), cube_coordinates="angstrom").fields
        self.assertEqual(fields[0].origin, (1, 2, 3))
        self.assertEqual(fields[1].values, [2] * 8)

    def test_xsf_order_endpoints_multiple_fields(self):
        text = "CRYSTAL\nPRIMVEC\n2 0 0\n0 2 0\n0 0 2\nPRIMCOORD\n1 1\nH 0 0 0\nBEGIN_BLOCK_DATAGRID_3D\ndata\n"
        for label in ("charge", "spin"):
            text += "BEGIN_DATAGRID_3D_{}\n2 2 2\n1 2 3\n2 0 0\n.5 2 0\n0 0 2\n1 2 3 4 5 6 7 8\nEND_DATAGRID_3D\n".format(label)
        text += "END_BLOCK_DATAGRID_3D\n"
        fields = load_volume(self.write("two.xsf", text), quantity="density").fields
        self.assertEqual(len(fields), 2)
        self.assertEqual(fields[0].values, list(range(1, 9)))
        self.assertFalse(fields[0].periodic)
        self.assertAlmostEqual(fields[0].integrate(), 36)
        self.assertEqual(fields[0].origin, (1, 2, 3))

    def test_xsf_molecular_atoms_and_unnamed_grid(self):
        text = "ATOMS\nH 1 2 3\n8 4 5 6\nBEGIN_BLOCK_DATAGRID_3D\ndata\nBEGIN_DATAGRID_3D\n2 2 2\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n1 1 1 1 1 1 1 1\nEND_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D\n"
        g = load_volume(self.write("molecule.xsf", text)).fields[0]
        self.assertEqual(g.name, "field")
        self.assertEqual(g.sites, [(1, 1, 2, 3), (8, 4, 5, 6)])

    def test_roundtrip_all_formats_periodic_boundary(self):
        original = self.grid()
        for format, name in (("cube", "out.cube"), ("xsf", "out.xsf"), ("vasp", "CHGCAR")):
            path = self.root / name
            original.save(path, format)
            loaded = load_volume(path, quantity="density").fields[0]
            self.assertAlmostEqual(loaded.integrate(), 16)
            if format != "vasp":
                self.assertEqual(loaded.shape, (5, 5, 5))
                loaded = loaded.as_periodic()
            self.assertEqual(loaded.shape, original.shape)
            for a, b in zip(loaded.values, original.values):
                self.assertAlmostEqual(a, b)

    def test_roundtrip_finite_nonconstant(self):
        grid = Grid((3, 4, 5), ((2, 0, 0), (.5, 3, 0), (0, .2, 4)), range(60), origin=(1, 2, 3))
        for format in ("cube", "xsf"):
            path = self.root / ("finite." + format)
            grid.save(path)
            restored = load_volume(path).fields[0]
            self.assertEqual(restored.values, grid.values)
            self.assertAlmostEqual(restored.integrate(), grid.integrate())
        with self.assertRaises(ValueError):
            grid.save(self.root / "invalid.vasp")
        with self.assertRaises(ValueError):
            grid.as_periodic()

    def test_arithmetic_derivatives_energy_and_resampling(self):
        a, b = self.grid(), self.grid(value=3)
        self.assertAlmostEqual((a + b).integrate(), 40)
        self.assertAlmostEqual((a - b).integrate(), -8)
        self.assertAlmostEqual((a * b).integrate(), 48)
        self.assertAlmostEqual((a / b).integrate(), 16/3)
        self.assertEqual((a / b).unit, "1")
        self.assertAlmostEqual((3 * a / 2).integrate(), 24)
        self.assertEqual(a.resample(self.grid(n=8)).shape, (8, 8, 8))
        self.assertAlmostEqual(a.smooth(.5, 1).integrate(), 16)
        for g in a.gradient():
            self.assertEqual(g.values, [0] * 64)
        self.assertEqual(a.laplacian().values, [0] * 64)
        kinetic, potential, total = a.energy_density()
        self.assertAlmostEqual(potential.values[0], -2 * kinetic.values[0])
        self.assertAlmostEqual(total.values[0], -kinetic.values[0])
        self.assertEqual(total.unit, "eV/A^3")
        self.assertEqual(a.values, [2] * 64)
        with self.assertRaises(ValueError):
            a / self.grid(value=0)
        with self.assertRaises(ValueError):
            a + self.grid(n=5)

    def test_bonding_descriptors_quadratic_density(self):
        # rho is an exact quadratic form, so the native engine's finite-
        # difference Hessian (gradient() differentiated a second time) is
        # exact at any interior point; the deep numerical correctness of the
        # C++ implementation is covered by tests/electronic_regressions.cpp,
        # so this test focuses on the Python<->native binding: operation
        # names, parameter marshaling and unit propagation.
        import math
        n, span = 6, 4.0
        h = span / (n - 1)  # finite grids include both endpoints
        A, B, C, D, rho0 = 0.01, -0.02, 0.03, 0.005, 5.0

        def analytic(x, y, z):
            return rho0 + A*x**2 + B*y**2 + C*z**2 + D*x*y

        values = [analytic(x*h, y*h, z*h) for z in range(n) for y in range(n) for x in range(n)]
        density = Grid((n, n, n), ((span, 0, 0), (0, span, 0), (0, 0, span)), values, unit="e/A^3")

        i0 = 3
        index = (i0*n + i0)*n + i0
        X0 = Y0 = Z0 = i0*h
        rho = analytic(X0, Y0, Z0)
        gx, gy, gz = 2*A*X0 + D*Y0, 2*B*Y0 + D*X0, 2*C*Z0
        # Hessian [[2A,D,0],[D,2B,0],[0,0,2C]]; middle eigenvalue of the 2x2
        # (x,y) block combined with the decoupled 2C z-eigenvalue.
        half_trace = A + B
        half_gap = math.sqrt((A - B)**2 + (D/2)**2)
        block_eigs = sorted([half_trace - half_gap, half_trace + half_gap])
        lambda2 = sorted([block_eigs[0], block_eigs[1], 2*C])[1]

        s_expected = math.sqrt(gx**2 + gy**2 + gz**2) / (2 * (3*math.pi**2)**(1/3) * rho**(4/3))
        rdg = density.reduced_density_gradient()
        self.assertAlmostEqual(rdg.values[index], s_expected, places=6)
        self.assertEqual(rdg.unit, "dimensionless")

        signed_expected = 0.0 if lambda2 == 0 else math.copysign(rho, lambda2)
        signed = density.signed_density()
        self.assertAlmostEqual(signed.values[index], signed_expected, places=6)
        self.assertEqual(signed.unit, "e/A^3")

        g, H = (gx, gy, gz), ((2*A, D, 0), (D, 2*B, 0), (0, 0, 2*C))
        k = (gx**2 + gy**2 + gz**2) / rho**2
        grad_gamma_sq = sum(((H[a][b]*rho - g[a]*g[b]) / rho**2)**2 for a in range(3) for b in range(3))
        theta = grad_gamma_sq / k**3
        dori = density.dori()
        self.assertAlmostEqual(dori.values[index], theta/(1+theta), places=5)
        self.assertEqual(dori.unit, "dimensionless")

        with self.assertRaises(ValueError):
            density.reduced_density_gradient(floor=0)
        wrong_unit = Grid(density.shape, density.cell, density.values, unit="eV")
        with self.assertRaises(ValueError):
            wrong_unit.reduced_density_gradient()
        negative = Grid((2, 2, 2), ((1, 0, 0), (0, 1, 0), (0, 0, 1)), [-1]*8, unit="e/A^3")
        with self.assertRaises(ValueError):
            negative.dori()
        # A perfectly constant density has zero gradient AND zero Hessian, so
        # DORI's theta = 0/0 is a genuine indeterminate limit, not a bug.
        with self.assertRaises(ValueError):
            self.grid().dori()

    def test_wannier_bindings_qwz_model(self):
        # Qi-Wu-Zhang two-band Chern-insulator model, H(k) = sin(2*pi*kx)*sx +
        # sin(2*pi*ky)*sy + (m+cos(2*pi*kx)+cos(2*pi*ky))*sz, written out as a
        # nearest-neighbour Wannier90 seedname_hr.dat fixture. The deep
        # numerical correctness (Kubo curvature vs. an independent plaquette
        # check, and vs. the gauge-invariant Fukui-Hatsugai-Suzuki Chern
        # number) is covered by tests/wannier_regressions.cpp; this test
        # focuses on the Python<->native binding: parsing, shapes and errors.
        from atomforge.electronic import read_hr, WannierHamiltonian

        def hr_text(m):
            blocks = {
                (0, 0, 0): [[m, 0], [0, -m]],
                (1, 0, 0): [[0.5, -0.5j], [-0.5j, -0.5]],
                (-1, 0, 0): [[0.5, 0.5j], [0.5j, -0.5]],
                (0, 1, 0): [[0.5, -0.5], [0.5, -0.5]],
                (0, -1, 0): [[0.5, 0.5], [-0.5, -0.5]],
            }
            lines = ["fixture", "2", str(len(blocks)), " ".join(["1"] * len(blocks))]
            for R, matrix in blocks.items():
                for row in range(2):
                    for col in range(2):
                        value = complex(matrix[row][col])
                        lines.append("{} {} {} {} {} {} {}".format(
                            R[0], R[1], R[2], row + 1, col + 1, value.real, value.imag))
            return "\n".join(lines) + "\n"

        path = self.root / "qwz_hr.dat"
        m = -1.0
        path.write_text(hr_text(m))
        hr = read_hr(path)
        self.assertIsInstance(hr, WannierHamiltonian)
        self.assertEqual(hr.num_wann, 2)

        kx, ky = 0.13, -0.27
        sx, sy, sz = ((0, 1), (1, 0)), ((0, -1j), (1j, 0)), ((1, 0), (0, -1))
        d = (math.sin(2 * math.pi * kx), math.sin(2 * math.pi * ky),
             m + math.cos(2 * math.pi * kx) + math.cos(2 * math.pi * ky))
        expected = [d[0] * sx[r][c] + d[1] * sy[r][c] + d[2] * sz[r][c] for r in range(2) for c in range(2)]
        h = hr.hamiltonian((kx, ky, 0))
        for actual, want in zip(h, expected):
            self.assertAlmostEqual(actual, want, places=10)

        norm = math.sqrt(sum(v * v for v in d))
        bands = hr.bands([(kx, ky, 0), (0, 0, 0)])
        self.assertEqual(len(bands), 2)
        self.assertEqual(len(bands[0]), 2)
        self.assertAlmostEqual(sorted(bands[0])[0], -norm, places=10)
        self.assertAlmostEqual(sorted(bands[0])[1], norm, places=10)

        curvature = hr.berry_curvature((kx, ky, 0))
        self.assertEqual(len(curvature), 2)
        self.assertAlmostEqual(curvature[0] + curvature[1], 0, places=8)  # 2-band identity

        chern = hr.chern_number(band=0, grid=24)
        self.assertAlmostEqual(abs(chern), 1, places=6)  # m=-1 is in the nontrivial (-2,0) regime
        path.write_text(hr_text(3.0))  # |m|>2 is the topologically trivial regime
        self.assertAlmostEqual(read_hr(path).chern_number(band=0, grid=24), 0, places=6)

        with self.assertRaises(ValueError):
            read_hr(self.root / "does_not_exist_hr.dat")
        with self.assertRaises(ValueError):
            hr.berry_curvature((kx, ky, 0), plane=(0, 0))
        with self.assertRaises(ValueError):
            hr.chern_number(band=5)
        # At m=-2 the gap closes exactly at k=0: Kubo curvature must diverge.
        path.write_text(hr_text(-2.0))
        with self.assertRaises(ValueError):
            read_hr(path).berry_curvature((0, 0, 0))

    def test_lobster_bindings_cohpcar_and_icohplist(self):
        # Deep numerical correctness (column layout, spin blocks, orbital-
        # resolved rejection) is covered by tests/lobster_regressions.cpp;
        # this test focuses on the Python<->native binding: dict shapes and
        # value marshaling for a small non-spin-polarized fixture.
        from atomforge.electronic import read_cohpcar, read_icohplist

        cohpcar_path = self.root / "COHPCAR.lobster"
        cohpcar_path.write_text(
            "created by AtomForge test fixture\n"
            "3 1 -10.0 10.0 -2.5\n"
            "No.1  COHP  ICOHP\n"
            "No.1:Fe1->Fe9(2.45)\n"
            "No.2:Fe1->O12(1.9)\n"
            "-1.0 0.1 0.01 0.2 0.02 0.3 0.03\n"
            "0.0 0.4 0.04 0.5 0.05 0.6 0.06\n"
            "1.0 0.7 0.07 0.8 0.08 0.9 0.09\n")
        data = read_cohpcar(cohpcar_path)
        self.assertFalse(data["spin_polarized"])
        self.assertAlmostEqual(data["fermi_energy_eV"], -2.5)
        self.assertEqual(data["energies_eV"], [-1.0, 0.0, 1.0])
        self.assertEqual(data["bonds"], [(0, 8, 2.45), (0, 11, 1.9)])
        self.assertEqual(data["average"][0]["cohp"], [0.1, 0.4, 0.7])
        self.assertEqual(data["bond_curves"][0][0]["cohp"], [0.2, 0.5, 0.8])
        self.assertEqual(data["bond_curves"][0][1]["icohp"], [0.03, 0.06, 0.09])
        with self.assertRaises(ValueError):
            read_cohpcar(self.root / "does_not_exist_COHPCAR.lobster")

        icohplist_path = self.root / "ICOHPLIST.lobster"
        icohplist_path.write_text(
            "COHP#  atom1  atom2  distance  translation  ICOHP(eV)\n"
            "1 Fe1 Fe9 2.45000 0 0 0 -1.23456\n"
            "2 Fe1 O12 1.90000 0 0 0 -0.54321\n\n")
        listing = read_icohplist(icohplist_path)
        self.assertFalse(listing["spin_polarized"])
        self.assertEqual(len(listing["entries"]), 2)
        first = listing["entries"][0]
        self.assertEqual((first["atom1"], first["atom2"]), (0, 8))
        self.assertAlmostEqual(first["length_A"], 2.45)
        self.assertAlmostEqual(first["icohp_up"], -1.23456)
        self.assertAlmostEqual(first["icohp_down"], 0.0)
        with self.assertRaises(ValueError):
            read_icohplist(self.root / "does_not_exist_ICOHPLIST.lobster")

    def test_betti_curve_bindings_hollow_shell(self):
        # Deep correctness (Euler characteristic via the doubled-coordinate
        # cubical complex, connected-component/cavity detection, several
        # known test shapes including a solid ring with betti1=1) is covered
        # by tests/topology_regressions.cpp; this test focuses on the
        # Python<->native binding: a hollow 3x3x3 shell (one enclosed
        # cavity) in a 5x5x5 finite grid, homotopy equivalent to a sphere.
        n = 5
        values = []
        for z in range(n):
            for y in range(n):
                for x in range(n):
                    inside = 1 <= x <= 3 and 1 <= y <= 3 and 1 <= z <= 3
                    center = (x, y, z) == (2, 2, 2)
                    values.append(1.0 if (inside and not center) else 0.0)
        grid = Grid((n, n, n), ((1, 0, 0), (0, 1, 0), (0, 0, 1)), values, periodic=False, unit="raw")

        rows = grid.betti_curve([0.5])
        self.assertEqual(len(rows), 1)
        threshold, betti0, betti1, betti2 = rows[0]
        self.assertAlmostEqual(threshold, 0.5)
        self.assertAlmostEqual(betti0, 1)
        self.assertAlmostEqual(betti1, 0)
        self.assertAlmostEqual(betti2, 1)  # the sealed 1-voxel cavity

        rows = grid.betti_curve([0.5, 2.0])
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(rows[1][2], 0)  # nothing is solid above the shell's density

        with self.assertRaises(ValueError):
            grid.betti_curve([])
        periodic = Grid((n, n, n), ((1, 0, 0), (0, 1, 0), (0, 0, 1)), values, periodic=True, unit="raw")
        with self.assertRaises(ValueError):
            periodic.betti_curve([0.5])

    def test_bader_partition_bindings_two_gaussians(self):
        # Deep numerical correctness (charge conservation, per-voxel
        # nearest-maximum assignment, periodic wraparound) is covered by
        # tests/bader_regressions.cpp; this test focuses on the
        # Python<->native binding: shapes, dict fields and error paths for
        # the same two-Gaussian fixture.
        n = 8
        c1, c2, sigma = (2.0, 3.3, 3.4), (5.0, 3.3, 3.4), 0.7
        values = []
        for z in range(n):
            for y in range(n):
                for x in range(n):
                    d1 = math.sqrt(sum((a - b) ** 2 for a, b in zip((x, y, z), c1)))
                    d2 = math.sqrt(sum((a - b) ** 2 for a, b in zip((x, y, z), c2)))
                    values.append(math.exp(-d1 * d1 / (2 * sigma ** 2)) + math.exp(-d2 * d2 / (2 * sigma ** 2)))
        grid = Grid((n, n, n), ((8, 0, 0), (0, 8, 0), (0, 0, 8)), values, periodic=False, unit="e/A^3")

        partition = grid.bader_partition()
        self.assertEqual(partition.num_basins, 2)
        first, second = partition.basin(0), partition.basin(1)
        self.assertAlmostEqual(first["charge"], second["charge"], places=6)
        self.assertAlmostEqual(first["charge"] + second["charge"], grid.integrate(), places=6)
        self.assertEqual(len(partition.basin_ids), n ** 3)
        self.assertEqual(set(partition.basin_ids), {0, 1})

        populations = partition.populations([first["maximum"], second["maximum"]])
        self.assertAlmostEqual(populations[0], first["charge"], places=10)
        self.assertAlmostEqual(populations[1], second["charge"], places=10)

        with self.assertRaises(ValueError):
            partition.basin(5)
        with self.assertRaises(ValueError):
            partition.populations([])
        with self.assertRaises(ValueError):
            self.grid(value=-1, periodic=False).bader_partition()
        with self.assertRaises(ValueError):
            Grid(grid.shape, grid.cell, grid.values, unit="eV").bader_partition()

    def test_profiles_sections_contours_and_surface_export(self):
        grid = Grid((4, 4, 4), ((3, 0, 0), (0, 3, 0), (0, 0, 3)), [x for z in range(4) for y in range(4) for x in range(4)])
        rows = grid.line_profile((0, 0, 0), (3, 0, 0), 7)
        for distance, value in rows:
            self.assertAlmostEqual(distance, value)
        self.assertEqual(grid.planar_average(0), [(0, 0), (1, 1), (2, 2), (3, 3)])
        self.assertAlmostEqual(grid.macroscopic_average(0, 3)[1][1], 1)
        section = grid.section((0, 0, 0), (3, 0, 0), (0, 3, 0), (4, 4))
        self.assertEqual(section.section_values[0], [0, 1, 2, 3])
        segments = section.contours(1.5)
        self.assertTrue(segments)
        for segment in segments:
            for x, _, _ in segment:
                self.assertAlmostEqual(x, 1.5)
        surface = grid.isosurface(1.5, color=grid)
        self.assertTrue(surface.vertices)
        for x, _, _, color in surface.vertices:
            self.assertAlmostEqual(x, 1.5)
            self.assertAlmostEqual(color, 1.5)
        path = self.root / "surface.obj"
        surface.save(path)
        self.assertIn("\nf ", path.read_text())
        surface.save(self.root / "surface.ply")
        self.assertIn("property double scalar", (self.root / "surface.ply").read_text())
        self.assertEqual(miller_section(self.grid(), (0, 0, 1), shape=(5, 5)).shape, (5, 5, 2))

    def test_integration_peaks_voronoi_spheres(self):
        g = self.grid()
        self.assertEqual(g.peaks(), [])
        basins = g.voronoi_integrate([(0, 0, 0), (1, 0, 0)])
        self.assertAlmostEqual(sum(v[0] for v in basins), 16)
        self.assertAlmostEqual(sum(v[1] for v in basins), 8)
        self.assertAlmostEqual(g.integrate_sphere((0, 0, 0), 10), 16)
        values = [0] * 64
        values[0] = 10
        peak = Grid(g.shape, g.cell, values, periodic=True)
        self.assertEqual(peak.peaks(), [(0, 0, 0, 10)])
        # Three coincident sites share each sample equally.
        for integral, volume in g.voronoi_integrate([(0, 0, 0)] * 3):
            self.assertAlmostEqual(integral, 16/3)
            self.assertAlmostEqual(volume, 8/3)

    def test_fourier_and_models(self):
        g = self.grid()
        factors = g.structure_factors([(0, 0, 0), (1, 0, 0), (-1, 0, 0)])
        self.assertAlmostEqual(factors[0][3], 16)
        self.assertAlmostEqual(factors[1][3], 0)
        self.assertAlmostEqual(g.fourier_synthesis(factors).integrate(), 16)
        self.assertAlmostEqual(g.patterson().values[0], 32)
        atomic = g.atomic_structure_factors([(0, 0, 0), (1, 0, 0)], [[1], [2]], occupancies=[.5], b_factors=[0])
        self.assertEqual(atomic, [(0, 0, 0, .5, 0), (1, 0, 0, 1, 0)])
        self.assertAlmostEqual(scattering_factor(1, [2], [1], 3), 2/math.e + 3)
        density = model_density(g, coefficients={1: ([1], [1], 0)})
        self.assertAlmostEqual(density.integrate(), 1)
        nuclear = model_density(g, scattering_lengths={1: -3})
        self.assertAlmostEqual(nuclear.integrate(), -3)
        self.assertEqual(nuclear.unit, "scattering-length/A^3")
        with self.assertRaises(ValueError):
            g.fourier_synthesis([(1, 0, 0, 1, 0)])
        with self.assertRaises(ValueError):
            g.fourier_synthesis([(0, 0, 0, 1, 0)] * 2)

    def test_ewald_python_parity(self):
        sites = [(1, 0, 0, 0), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1), (1, 1, 1, 1)]
        g = self.grid().with_sites(sites)
        result = g.ewald([1]*4 + [-1]*4, 1.2, 8, 16)
        expected = -1.747564594633182 * 27.211386245988 * .529177210903
        self.assertAlmostEqual(result["energy"], 4 * expected, places=8)
        self.assertAlmostEqual(result["potentials"][0], expected, places=8)
        with self.assertRaises(ValueError):
            self.grid().ewald([1], 1, 8, 12)

    def test_invalid_inputs_and_native_error_recovery(self):
        g = self.grid()
        operations = [lambda: g.smooth(-1), lambda: g.planar_average(4), lambda: g.line_profile((0, 0, 0), (1, 0, 0), 1),
                      lambda: g.macroscopic_average(0, 2), lambda: g.voronoi_integrate([]), lambda: g * float("nan"),
                      lambda: g.integrate_sphere((0, 0, 0), -1), lambda: g.isosurface(float("inf")),
                      lambda: Grid((2, 2, 2), ((1, 0, 0), (2, 0, 0), (0, 0, 1)), [1]*8)]
        for operation in operations:
            with self.assertRaises(ValueError):
                operation()
            self.assertEqual(g.integrate(), 16)

    def test_cli(self):
        path = self.root / "CHGCAR"
        self.grid().save(path)
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
        result = subprocess.run([sys.executable, "-m", "atomforge.electronic", str(path), "--operation", "integrate"], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(float(result.stdout), 16)
        invalid = subprocess.run([sys.executable, "-m", "atomforge.electronic", str(path), "--channel", "-1"], env=env, capture_output=True, text=True)
        self.assertNotEqual(invalid.returncode, 0)


if __name__ == "__main__":
    unittest.main()
