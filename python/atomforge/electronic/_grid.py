"""Typed, immutable grid operations backed by the same engine as the desktop UI."""

import ctypes as ct
import math
from pathlib import Path

from ._native import Handle, doubles, library


def _vector(value, size=3):
    value = tuple(value)
    if len(value) != size or not all(math.isfinite(x) for x in value):
        raise ValueError("Expected {} finite components".format(size))
    return value


def _flatten(rows, columns):
    return [v for row in rows for v in _vector(row, columns)]


class Grid:
    """Scalar grid in Angstrom; values use x-fastest order (z, y, x loops).

    ``cell`` contains three lattice/span vectors, as in ``Structure.cell``.
    Periodic grids omit endpoints; finite grids include them. Operations return
    new grids. ``unit='e/A^3'`` enables density-derived energy calculations.
    """

    def __init__(self, shape, cell, values, origin=(0, 0, 0), periodic=False, unit="raw"):
        shape = tuple(shape)
        if len(shape) != 3 or any(not isinstance(n, int) or n < 2 for n in shape):
            raise ValueError("shape must contain three integers >= 2")
        if math.prod(shape) > 100000000:
            raise ValueError("Grid exceeds 100 million samples")
        flat_cell = _flatten(cell, 3)
        if len(flat_cell) != 9:
            raise ValueError("cell must have three vectors")
        data = doubles(values)
        if len(data) != math.prod(shape):
            raise ValueError("Value count differs from shape")
        self._handle = Handle(library().af_electronic_create(
            (ct.c_int * 3)(*shape), doubles(flat_cell), doubles(_vector(origin)),
            int(bool(periodic)), data, len(data), str(unit).encode("utf-8"),
        ))

    @classmethod
    def _from_handle(cls, handle):
        grid = cls.__new__(cls)
        grid._handle = handle
        return grid

    def _geometry(self):
        shape, cell, origin = (ct.c_int * 3)(), (ct.c_double * 9)(), (ct.c_double * 3)()
        periodic = self._handle.lib.af_electronic_geometry(self._handle.pointer, shape, cell, origin)
        return tuple(shape), tuple(tuple(cell[i:i + 3]) for i in (0, 3, 6)), tuple(origin), bool(periodic)

    @property
    def shape(self):
        return self._geometry()[0]

    @property
    def cell(self):
        return self._geometry()[1]

    @property
    def origin(self):
        return self._geometry()[2]

    @property
    def periodic(self):
        return self._geometry()[3]

    @property
    def values(self):
        """Return an independent list of scalar samples."""
        return self._handle.data()

    @property
    def unit(self):
        return self._handle.lib.af_electronic_label(self._handle.pointer, 0, 1).decode("utf-8")

    @property
    def name(self):
        return self._handle.lib.af_electronic_label(self._handle.pointer, 0, 0).decode("utf-8")

    @property
    def sites(self):
        """Return (atomic_number, x, y, z) rows; unknown species have number 0."""
        count = self._handle.lib.af_electronic_sites(self._handle.pointer, None)
        data = (ct.c_double * (count * 4))()
        self._handle.lib.af_electronic_sites(self._handle.pointer, data)
        return [(int(data[i]), *data[i + 1:i + 4]) for i in range(0, len(data), 4)]

    def with_sites(self, sites):
        data = doubles(_flatten(sites, 4))
        return Grid._from_handle(Handle(self._handle.lib.af_electronic_with_sites(self._handle.pointer, data, len(data) // 4)))

    def _calculate(self, operation, parameters=(), other=None):
        return self._handle.calculate(operation, parameters, other._handle if other is not None else None)

    def _field(self, operation, parameters=(), other=None):
        return Grid._from_handle(self._calculate(operation, parameters, other))

    def _fields(self, operation, parameters=()):
        return Volume(self._calculate(operation, parameters)).fields

    def __add__(self, other):
        return self._field("add", other=other)

    def __sub__(self, other):
        return self._field("subtract", other=other)

    def __mul__(self, other):
        return self._field("multiply", other=other) if isinstance(other, Grid) else self._field("scale", [other])

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Grid):
            return self._field("divide", other=other)
        if other == 0:
            raise ValueError("Division by zero")
        return self * (1 / other)

    def resample(self, target):
        return self._field("resample", other=target)

    def density_difference(self, *references, weights=None):
        """Return self - sum(weight * reference), in e/A^3.

        Use fragment densities calculated in the same cell and geometry.
        Alignment is required; resample explicitly when appropriate.
        """
        if not references:
            raise ValueError("Provide at least one reference density")
        weights = [1] * len(references) if weights is None else list(weights)
        if len(weights) != len(references):
            raise ValueError("Provide one weight per reference")
        result = self
        for reference, weight in zip(references, weights):
            result = result._field("density_difference", [weight], reference)
        return result

    def threshold_mask(self, low, high):
        """Binary mask for the inclusive interval [low, high]."""
        return self._field("threshold_mask", [low, high])

    def boolean(self, other, operation="intersection"):
        """Combine binary masks: union, intersection, difference (self-other), xor."""
        if operation not in ("union", "intersection", "difference", "xor"):
            raise ValueError("Unknown Boolean mask operation")
        return self._field(operation, other=other)

    def invert_mask(self):
        """Complement a binary mask within its cell."""
        return self.threshold_mask(0, 1).boolean(self, "difference")

    def apply_mask(self, mask):
        """Keep density inside a binary mask, preserving the density unit."""
        return self._field("apply_mask", other=mask)

    def split_density(self):
        """Return accumulation and positive depletion-magnitude fields."""
        return self._fields("split_density")

    def charge_summary(self):
        """Electron accumulation, depletion magnitude and net change.

        Apply to a difference density. These integrals describe redistribution;
        they do not assign charge transfer to individual atoms or fragments.
        """
        values = self._calculate("charge_summary").rows()[0]
        return dict(zip(("accumulation", "depletion", "net"), values))

    def cumulative_charge(self, axis=2):
        """(Normal distance in A, cumulative electrons) from the lower cell face.

        The final row includes the full cell. Periodic profiles depend on the
        chosen cell origin; finite grids use trapezoidal quadrature.
        """
        return self._calculate("cumulative_charge", [axis]).rows()

    def as_periodic(self, tolerance=1e-8):
        """Drop matching endpoint planes after checking all three boundaries.

        Only use when the finite grid spans a complete periodic cell. This
        conversion checks values, but cannot establish physical periodicity.
        """
        return self._field("periodic", [tolerance])

    def gradient(self):
        """Return Cartesian x/y/z derivatives, including triclinic metric."""
        return self._fields("gradient")

    def laplacian(self):
        return self._field("laplacian")

    def energy_density(self, floor=1e-12):
        """Return approximate kinetic, potential and total density in eV/A^3.

        Nonnegative electron density is required; samples <= floor (e/A^3)
        are masked to zero. These are gradient-expansion/local-virial fields,
        not the kinetic energy evaluated from Kohn-Sham orbitals.
        """
        return self._fields("energy", [floor])

    def smooth(self, sigma, radius=2):
        return self._field("smooth", [sigma, radius])

    def line_profile(self, start, end, count=200):
        return self._calculate("line", [*_vector(start), *_vector(end), count]).rows()

    def planar_average(self, axis=2):
        """Return normal distance in Angstrom and area-averaged field."""
        return self._calculate("planar", [axis]).rows()

    def macroscopic_average(self, axis=2, window=3):
        """Centered box average over an odd number of planar samples."""
        return self._calculate("macro", [axis, window]).rows()

    def section(self, origin, u, v, shape=(100, 100)):
        """Sample origin + s*u + t*v for 0<=s,t<=1 (Cartesian Angstrom).

        Returns a two-layer grid with identical layers; the first layer is
        the section. Use section_values/contours for 2D data, not its volume.
        """
        return self._field("section", [*_vector(origin), *_vector(u), *_vector(v), *_vector(shape, 2)])

    @property
    def section_values(self):
        nx, ny, _ = self.shape
        data = self.values
        return [data[y * nx:(y + 1) * nx] for y in range(ny)]

    def contours(self, level):
        """Return endpoint pairs for contour segments on grid z=0."""
        rows = self._calculate("contours", [level]).rows()
        return list(zip(rows[::2], rows[1::2]))

    def integrate(self):
        return self._calculate("integrate").rows()[0][0]

    def integrate_sphere(self, center, radius):
        return self._calculate("sphere", [*_vector(center), radius]).rows()[0][0]

    def voronoi_integrate(self, sites=None):
        """Return (integral, volume) per site; equidistant samples share weight."""
        if sites is None:
            sites = [row[1:] for row in self.sites]
        return self._calculate("voronoi", _flatten(sites, 3)).rows()

    def peaks(self, limit=0):
        """Grid-local maxima as (x,y,z,value); deterministic adjacent tie break."""
        return self._calculate("peaks", [limit]).rows()

    def structure_factors(self, indices):
        """Return (h,k,l,real,imag); phase convention is exp(+2*pi*i*h.r)."""
        return self._calculate("factors", _flatten(indices, 3)).rows()

    def fourier_synthesis(self, reflections):
        """Supply complete (h,k,l,real,imag) rows, including conjugate pairs."""
        return self._field("synthesis", _flatten(reflections, 5))

    def patterson(self):
        return self._field("patterson")

    def atomic_structure_factors(self, indices, scattering, occupancies=None, b_factors=None):
        """Explicit complex scattering[reflection][site], occupancy and B(A^2)."""
        count = len(self.sites)
        occupancies = [1] * count if occupancies is None else list(occupancies)
        b_factors = [0] * count if b_factors is None else list(b_factors)
        if len(occupancies) != count or len(b_factors) != count:
            raise ValueError("Expected occupancy and B for each site")
        indices, scattering = list(indices), list(scattering)
        if len(indices) != len(scattering):
            raise ValueError("Expected scattering row for each reflection")
        params = [v for pair in zip(occupancies, b_factors) for v in pair]
        for h, factors in zip(indices, scattering):
            params.extend(_vector(h))
            factors = list(factors)
            if len(factors) != count:
                raise ValueError("Expected scattering factor for each site")
            for f in factors:
                f = complex(f)
                params.extend((f.real, f.imag))
        return self._calculate("atomic_factors", params).rows()

    def ewald(self, charges, alpha, real_cutoff, reciprocal_cutoff):
        """Neutral periodic point charges; return energy(eV), potentials(V).

        alpha is in 1/A, real_cutoff in A and reciprocal_cutoff in 1/A
        (wavevector including 2*pi). Check convergence by increasing cutoffs.
        """
        if not self.periodic:
            raise ValueError("Ewald requires a periodic cell")
        rows = self._calculate("ewald", [alpha, real_cutoff, reciprocal_cutoff, *charges]).rows()
        return {"energy": rows[0][0], "potentials": [r[0] for r in rows[1:]]}

    def isosurface(self, level, color=None):
        return Surface(self._calculate("isosurface", [level], color))

    def save(self, path, format=None):
        format = format or ({".cube": "cube", ".cub": "cube", ".xsf": "xsf"}.get(Path(path).suffix.lower(), "vasp"))
        self._handle.save(path, format)


class Volume:
    """An imported collection of scalar channels and their atomic sites."""

    def __init__(self, handle):
        self._handle = handle

    @property
    def fields(self):
        lib, pointer = self._handle.lib, self._handle.pointer
        return tuple(Grid._from_handle(Handle(lib.af_electronic_select(pointer, i))) for i in range(lib.af_electronic_fields(pointer)))

    def save(self, path, format):
        self._handle.save(path, format)


class Surface:
    """Triangles with per-vertex scalar coloring; vertices are not welded."""

    def __init__(self, handle):
        self._handle = handle

    @property
    def vertices(self):
        return self._handle.rows()

    def save(self, path):
        """Write OBJ (scalar comments) or .ply (a scalar vertex property)."""
        self._handle.save(path, "obj")


def load_volume(path, quantity="auto", cube_coordinates="bohr"):
    """Read VASP, Cube or XSF. Cube/XSF 'auto' keeps scalar units raw.

    quantity: density, potential, elf, raw, auto. Cube density/potential are
    converted from atomic units. XSF scalars are assumed e/A^3 or eV only
    when explicitly requested; producer-specific units must be converted.
    """
    return Volume(Handle(library().af_electronic_load(
        str(path).encode("utf-8"), quantity.encode("ascii"), cube_coordinates.encode("ascii"),
    )))
