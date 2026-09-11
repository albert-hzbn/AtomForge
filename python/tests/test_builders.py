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
