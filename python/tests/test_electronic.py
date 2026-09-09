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
