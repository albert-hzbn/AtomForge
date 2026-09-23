"""Native-backed regression tests for the public builder interface."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import atomforge as af


@unittest.skipUnless(os.environ.get('ATOMFORGE_PATH'), 'requires the native executable')
class BuilderTests(unittest.TestCase):
    def test_all_builders_and_in_memory_source(self):
        with tempfile.TemporaryDirectory(prefix='atomforge paths ') as folder:
            root = Path(folder)
            cu = af.build('bulk', ['--a', '3.61', '--atom', 'Cu 0 0 0']).structure
            self.assertEqual(len(cu), 4)
            host = cu.repeat(5, 5, 5)
            alloy = af.build('sss', ['--frac', 'Cu=0.7,Ni=0.3', '--seed', '42'],
                             source=host, output=root / 'alloy sample.cif')
            self.assertEqual(len(alloy.structure), 500)
            self.assertEqual(sum(a.symbol == 'Ni' for a in alloy.structure.atoms), 150)
            self.assertTrue(alloy.output.is_file())
            self.assertIn('Saved', alloy.stdout)
            sphere = af.build('nano', ['--shape', 'sphere', '--radius', '10'], source=cu)
            self.assertEqual(len(sphere.structure), 336)
            gb = af.build('gb', ['--sigma', '5', '--uca', '3', '--ucb', '3', '--overlap', '1.5'], source=cu)
            self.assertEqual(len(gb.structure), 60)
            poly = af.build('poly', ['--sizex', '25', '--sizey', '25', '--sizez', '25',
                                     '--grains', '4', '--seed', '7'], source=cu)
            self.assertEqual(len(poly.structure), 1329)
            glass = af.build('amorphous', ['--element', 'Si 40', '--element', 'O 80',
                                           '--density', '1.5', '--seed', '7'])
            self.assertEqual(len(glass.structure), 120)
            mesh = root / 'mesh with spaces.obj'
            mesh.write_text('v 1 0 0\nv -1 0 0\nv 0 1 0\nv 0 -1 0\nv 0 0 1\nv 0 0 -1\n'
                            'f 1 3 5\nf 3 2 5\nf 2 4 5\nf 4 1 5\nf 3 1 6\nf 2 3 6\nf 4 2 6\nf 1 4 6\n')
            filled = af.build('custom', ['--mesh', str(mesh), '--scale', '12'], source=cu)
            self.assertEqual(len(filled.structure), 224)
            defect = af.build('dislocation', ['--character', 'edge', '--shape', 'cylinder',
                                              '--cyl-radius', '5', '--core', '1.2', '--cutoff', '8'], source=host)
            self.assertEqual(len(defect.structure), 500)
            self.assertEqual(len(host), 500)
            self.assertEqual(len(cu), 4)
            self.assertIsNone(sphere.output)

    def test_anisotropic_dislocation(self):
        # Real Cu single-crystal elastic constants (GPa); A = 2*C44/(C11-C12)
        # = 3.2, strongly anisotropic, so this also exercises the Stroh
        # sextic solver's non-degenerate path (no noise workaround needed).
        cu = af.build('bulk', ['--a', '3.61', '--atom', 'Cu 0 0 0']).structure
        host = cu.repeat(10, 10, 10)
        aniso = af.build('dislocation', [
            '--character', 'edge', '--shape', 'cylinder', '--cyl-radius', '15',
            '--anisotropic', '--elastic-c11', '168.4', '--elastic-c12', '121.4', '--elastic-c44', '75.4',
        ], source=host)
        self.assertEqual(len(aniso.structure), 4000)
        self.assertIn('Inserted edge dislocation', aniso.stdout)
        self.assertNotIn('lattice family changed', aniso.stdout)
        for atom in aniso.structure.atoms:
            self.assertTrue(all(map(lambda v: v == v, (atom.x, atom.y, atom.z))))  # no NaN

        # Missing elastic constants must be rejected with a clear error, not
        # silently ignored or crash.
        with self.assertRaises(subprocess.CalledProcessError) as ctx:
            af.build('dislocation', ['--character', 'edge', '--anisotropic'], source=host)
        self.assertIn('elastic-c11', ctx.exception.stderr)

    def test_dislocation_dipole(self):
        # A dipole's net Burgers vector is zero, so unlike a single
        # dislocation it is periodicity-compatible; this exercises the
        # dipole superposition path (a second, opposite-sign dislocation
        # offset in the slip plane) independent of the anisotropic solver.
        cu = af.build('bulk', ['--a', '3.61', '--atom', 'Cu 0 0 0']).structure
        host = cu.repeat(10, 10, 10)
        dipole = af.build('dislocation', [
            '--character', 'edge', '--shape', 'halfplane', '--dipole', '--dipole-offset', '18 0',
        ], source=host)
        self.assertEqual(len(dipole.structure), 4000)
        self.assertIn('dislocation dipole', dipole.stdout)
        self.assertIn('Validation passed', dipole.stdout)

    def test_prepare_drag(self):
        # Uses a non-periodic (molecule-style) pair of configurations so the
        # comparison isn't entangled with AtomForge's native loader always
        # wrapping periodic structures back into their primary [0, cell)
        # cell on load (an app-wide, deliberate behavior used by every
        # native CLI operation, exercised precisely and directly -- with no
        # file round-trip -- by the native drag_prep_regressions test).
        initial = af.Structure()
        initial.add_atom('Cu', 0.0, 0.0, 0.0)
        initial.add_atom('Cu', 5.0, 0.0, 0.0)
        final = af.Structure()
        final.add_atom('Cu', 1.0, 2.0, -0.5)
        final.add_atom('Cu', 5.3, 0.0, 0.0)

        midpoint, directions = af.prepare_drag(initial, final, zeta=0.5)
        self.assertEqual(len(midpoint), 2)
        self.assertEqual(len(directions), 2)
        self.assertAlmostEqual(midpoint.atoms[0].x, 0.5, places=3)
        self.assertAlmostEqual(midpoint.atoms[0].y, 1.0, places=3)
        self.assertAlmostEqual(midpoint.atoms[0].z, -0.25, places=3)
        self.assertAlmostEqual(midpoint.atoms[1].x, 5.15, places=3)
        self.assertAlmostEqual(directions[0][0], 1.0, places=3)
        self.assertAlmostEqual(directions[0][1], 2.0, places=3)
        self.assertAlmostEqual(directions[0][2], -0.5, places=3)

        start, _ = af.prepare_drag(initial, final, zeta=0.0)
        self.assertAlmostEqual(start.atoms[0].x, 0.0, places=4)
        end, _ = af.prepare_drag(initial, final, zeta=1.0)
        self.assertAlmostEqual(end.atoms[0].x, 1.0, places=3)
        self.assertAlmostEqual(end.atoms[0].y, 2.0, places=3)

        with self.assertRaises(ValueError):
            af.prepare_drag(initial, final, zeta=1.5)

    def test_fit_dislocation(self):
        # Native, exact-line-direction validation (line = [0,0,1] for a
        # manually-specified screw dislocation) lives in
        # tests/dislocation_fit_regressions.cpp; this just exercises the
        # CLI/Python plumbing end to end.
        cu = af.build('bulk', ['--a', '3.61', '--atom', 'Cu 0 0 0']).structure
        host = cu.repeat(13, 13, 13)
        dislo = af.build('dislocation', ['--character', 'screw', '--shape', 'cylinder', '--cyl-radius', '15',
                                         '--manual-vectors', '--line', '0 0 1', '--burgers', '0 0 1'],
                          source=host).structure

        fit = af.fit_dislocation(dislo, host, (0, 0, 1), cutoff=3.0, no_pbc=True)
        self.assertIn('line_position', fit)
        self.assertIn('burgers_vector', fit)
        # This weighted-centroid/integral estimator is validated precisely
        # (position and Burgers direction) for an edge dislocation in the
        # native dislocation_fit_regressions.cpp test; here (a screw
        # configuration) just confirm the CLI/Python plumbing returns a
        # finite, non-trivial result end to end.
        burgers_mag = sum(c * c for c in fit['burgers_vector']) ** 0.5
        self.assertGreater(burgers_mag, 1e-6)
        self.assertTrue(all(abs(c) < 1e6 for c in fit['line_position']))

    def test_rotated_cell_keeps_coordinate_frame(self):
        source = af.Structure()
        source.cell = [[0, 4, 0], [-4, 0, 0], [0, 0, 4]]
        source.add_atom('Cu', -1, 2, 1)
        result = af.build('sss', ['--frac', 'Cu=1'], source=source).structure
        self.assertEqual(len(result), 1)
        for got, expected in zip(result.cell, source.cell):
            for value, target in zip(got, expected):
                self.assertAlmostEqual(value, target, places=5)
        self.assertAlmostEqual(result.atoms[0].x, -1, places=5)
        self.assertAlmostEqual(result.atoms[0].y, 2, places=5)

    def test_wulff_and_invalid_facets(self):
        cu = af.build('bulk', ['--a', '3.61', '--atom', 'Cu 0 0 0']).structure
        result = af.wulff(cu, [(1, 0, 0, 1), (1, 1, 1, 1)], radius=20)
        self.assertEqual(len(result.structure), 3552)
        self.assertIn('Wulff', result.stdout)
        for facet in [(0, 0, 0, 1), (1, 0, 0, -1), (1, 0, 0, float('nan'))]:
            with self.subTest(facet=facet), self.assertRaises(subprocess.CalledProcessError) as caught:
                af.wulff(cu, [facet])
            self.assertIn('facet', caught.exception.stderr)
        self.assertIn('--facet', af.builder_help('nano'))

    def test_input_validation_and_native_diagnostics(self):
        with self.assertRaises(TypeError):
            af.build('bulk', '--a 3.61')
        with self.assertRaises(ValueError):
            af.build('bulk', ['--output', 'unexpected.xyz'])
        with self.assertRaises(ValueError):
            af.build('invalid')
        with self.assertRaises(ValueError):
            af.build('bulk', timeout=0)
        with self.assertRaises(subprocess.CalledProcessError) as caught:
            af.build('bulk', ['--a', 'nan'])
        self.assertIn('invalid numeric', caught.exception.stderr)
        with self.assertRaises(FileNotFoundError):
            af.build('nano', source='missing-source-for-test.cif')


if __name__ == '__main__':
    unittest.main()
