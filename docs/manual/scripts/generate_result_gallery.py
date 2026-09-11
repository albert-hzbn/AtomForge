"""Plot the published native builder examples in 3D and as XY projections."""
import os
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[3]
MANUAL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
if os.environ.get('MANUAL_PYTHON_DEPS'):
    sys.path.insert(0, os.environ['MANUAL_PYTHON_DEPS'])
import atomforge as af
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

EXAMPLES = [
    ('bulk', 'cu_fcc.cif', 'FCC Cu: space group 225, a = 3.61 \u00c5'),
    ('alloy', 'cu_ni_alloy.cif', 'Cu 70% / Ni 30%, seed 42'),
    ('csl', 'cu_sigma5.cif', 'Sigma 5 [001], 3 + 3 repeats, overlap 1.5 \u00c5'),
    ('polycrystal', 'cu_polycrystal.cif', '25 \u00c5 cube, four grains, seed 7'),
    ('amorphous', 'sio2_pack.xyz', '40 Si + 80 O, density 1.5 g/cm3, seed 7'),
    ('custom', 'cu_mesh.xyz', 'Octahedron mesh, scale 12 \u00c5/model unit'),
]

def main():
    palette = {'Cu': '#ad652d', 'Ni': '#3b8996', 'Si': '#597aab', 'O': '#c94646'}
    for name, filename, title in EXAMPLES:
        structure = af.load(str(MANUAL / 'examples' / filename))
        xyz = np.array([(a.x, a.y, a.z) for a in structure.atoms])
        if len(xyz) == 0 or not np.isfinite(xyz).all():
            raise ValueError(f'Invalid example {filename}')
        colors = [palette.get(a.symbol, '#727889') for a in structure.atoms]
        fig = plt.figure(figsize=(9.2, 3.8), facecolor='white')
        ax = fig.add_subplot(121, projection='3d')
        ax.scatter(*xyz.T, c=colors, s=12 if len(xyz)>50 else 150, edgecolors='#423b35', linewidths=.15)
        ax.set_box_aspect(np.maximum(np.ptp(xyz, axis=0), 1))
        ax.view_init(elev=24, azim=35)
        ax.set_axis_off()
        ax.set_title('3D atomic geometry', fontsize=11)
        ax = fig.add_subplot(122)
        ax.scatter(xyz[:,0], xyz[:,1], c=colors, s=10 if len(xyz)>50 else 120, edgecolors='#423b35', linewidths=.15)
        ax.set_aspect('equal', adjustable='box')
        span = np.ptp(xyz[:, :2], axis=0)
        center = (xyz[:, :2].max(axis=0) + xyz[:, :2].min(axis=0)) / 2
        half_span = np.maximum(span, max(span.max() * .2, 1)) * .58
        ax.set_xlim(center[0] - half_span[0], center[0] + half_span[0])
        ax.set_ylim(center[1] - half_span[1], center[1] + half_span[1])
        ax.set_xlabel('x (\u00c5)'); ax.set_ylabel('y (\u00c5)')
        ax.set_title('2D projection along Z (all atoms)', fontsize=11)
        ax.spines[['top','right']].set_visible(False)
        fig.suptitle(f'{title} | {len(xyz):,} atoms', fontsize=12)
        fig.tight_layout()
        fig.savefig(MANUAL / 'figures/tutorials' / f'{name}-output.png', dpi=180)
        plt.close(fig)
        print(f'{name}: {len(xyz)} atoms')

if __name__ == '__main__':
    main()
