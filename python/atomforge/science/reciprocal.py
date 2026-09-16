"""BXSF reciprocal-cell band energies and Fermi-surface meshes.

BANDGRID uses z-fastest (C) order, unlike the real-space XSF DATAGRID.
Coordinates retain the reciprocal-vector units supplied by the producer.
Surfaces span the reciprocal parallelepiped; they are not folded into the
Wigner-Seitz first Brillouin zone.
"""

from dataclasses import dataclass
from pathlib import Path
import math
import re


@dataclass
class BandGrid:
    shape: tuple
    vectors: tuple
    origin: tuple
    bands: dict
    fermi_energy: float
    energy_unit: str
    reciprocal_unit: str

    def surface(self, band, energy=None):
        """Triangulate E(k)=energy; the default is the file's Fermi energy."""
        from ..electronic import Grid
        level = self.fermi_energy if energy is None else float(energy)
        values = self.bands[band]
        if not math.isfinite(level) or not min(values) <= level <= max(values):
            raise ValueError("Isosurface energy is outside the band's sampled range")
        nx, ny, nz = self.shape
        # Native geometric triangulation is coordinate-unit independent.
        ordered = [values[(x*ny+y)*nz+z]
                   for z in range(nz) for y in range(ny) for x in range(nx)]
        grid = Grid(self.shape, self.vectors, ordered, origin=self.origin,
                    periodic=False, unit=self.energy_unit)
        return grid.isosurface(level)


def read_bandgrid(path, *, energy_unit="eV", reciprocal_unit="1/A"):
    """Read one three-dimensional BXSF block, preserving band labels.

    BXSF does not declare universal physical units. Supply the producing
    program's units explicitly when they differ from eV and inverse Angstrom.
    Endpoint samples are retained, as required for general BXSF grids.
    """
    text = Path(path).read_text(encoding="utf-8")
    text = re.sub(r"#[^\n]*", "", text)
    fermi = re.search(r"Fermi\s+Energy\s*:\s*([^\s]+)", text, re.I)
    if fermi is None:
        raise ValueError("BXSF is missing Fermi Energy")
    starts = list(re.finditer(r"(?im)^\s*BEGIN_BANDGRID_3D[^\n]*\n", text))
    if len(starts)!=1:
        raise ValueError("Expected exactly one three-dimensional BANDGRID")
    tokens = iter(text[starts[0].end():].split())
    def real(token):
        value=float(token.replace("D","E").replace("d","e"))
        if not math.isfinite(value):
            raise ValueError("Non-finite BXSF value")
        return value
    try:
        count=int(next(tokens)); shape=tuple(int(next(tokens)) for _ in range(3))
        samples=math.prod(shape)
        if count<1 or any(n<2 for n in shape) or count*samples>100000000:
            raise ValueError("Unsupported BANDGRID size")
        origin=tuple(real(next(tokens)) for _ in range(3))
        vectors=tuple(tuple(real(next(tokens)) for _ in range(3)) for _ in range(3))
        bands={}
        for _ in range(count):
            if next(tokens).upper()!="BAND:":
                raise ValueError("Expected BAND: label")
            label=int(next(tokens))
            if label in bands:
                raise ValueError("Repeated band label")
            bands[label]=[real(next(tokens)) for _ in range(samples)]
        if next(tokens).upper()!="END_BANDGRID_3D":
            raise ValueError("Unexpected extra data or missing BANDGRID terminator")
        if next(tokens).upper()!="END_BLOCK_BANDGRID_3D":
            raise ValueError("Missing BANDGRID block terminator")
    except StopIteration as error:
        raise ValueError("Truncated BANDGRID") from error
    return BandGrid(shape,vectors,origin,bands,real(fermi.group(1)),energy_unit,reciprocal_unit)


def plot_surface(grid, band, output, *, energy=None):
    """Save a labeled reciprocal-space surface plot as PNG, SVG or PDF."""
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    mesh=grid.surface(band,energy)
    triangles=np.asarray(mesh.vertices)[:,:3].reshape(-1,3,3)
    if not len(triangles):
        raise ValueError("This level does not produce a surface")
    figure=Figure(); FigureCanvasAgg(figure); axes=figure.add_subplot(projection="3d")
    axes.add_collection3d(Poly3DCollection(triangles,alpha=.8,facecolor="steelblue"))
    points=triangles.reshape(-1,3)
    axes.auto_scale_xyz(*points.T)
    for setter,label in ((axes.set_xlabel,"kx"),(axes.set_ylabel,"ky"),(axes.set_zlabel,"kz")):
        setter(label+" ("+grid.reciprocal_unit+")")
    axes.set_title("Band {} — reciprocal-cell surface".format(band))
    figure.savefig(output,bbox_inches="tight")
