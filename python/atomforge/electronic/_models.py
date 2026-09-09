"""Explicit scattering models; no inferred oxidation states or element tables."""

import math


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def reciprocal_vectors(cell):
    """Reciprocal vectors without 2*pi, in inverse Angstrom."""
    a, b, c = cell
    volume = _dot(a, _cross(b, c))
    if not math.isfinite(volume) or abs(volume) < 1e-12:
        raise ValueError("Singular cell")
    return tuple(tuple(v / volume for v in vector) for vector in (_cross(b, c), _cross(c, a), _cross(a, b)))


def scattering_factor(s, a, b, c=0.0):
    """Gaussian X-ray form factor sum(a_i*exp(-b_i*s^2))+c.

    Supply published coefficients and respect their stated valid s=sin(theta)/
    wavelength range. No anomalous-dispersion correction is added implicitly.
    """
    a, b = list(a), list(b)
    if len(a) != len(b) or not a:
        raise ValueError("Expected equally sized, nonempty coefficient arrays")
    if not all(math.isfinite(v) for v in [s, c, *a, *b]) or s < 0:
        raise ValueError("Scattering parameters must be finite and s nonnegative")
    return sum(ai * math.exp(-bi * s * s) for ai, bi in zip(a, b)) + c


def model_density(geometry, coefficients=None, scattering_lengths=None, occupancies=None, b_factors=None):
    """Fourier model electron or nuclear scattering-length density.

    Exactly one mapping keyed by atomic number is required: coefficients
    ``{Z: (a, b, c)}`` for X-rays or scattering_lengths ``{Z: length}`` for
    neutrons. Length units are caller-defined. Even-grid Nyquist planes are
    omitted to retain a real, conjugate-symmetric finite Fourier expansion.
    These are free-atom scattering models, not self-consistent densities.
    """
    from ._grid import Grid

    if (coefficients is None) == (scattering_lengths is None):
        raise ValueError("Supply either X-ray coefficients or neutron scattering lengths")
    if not geometry.periodic:
        raise ValueError("Model density requires periodic geometry")
    ranges = [range(-((n - 1) // 2), (n - 1) // 2 + 1) for n in geometry.shape]
    indices = [(h, k, l) for l in ranges[2] for k in ranges[1] for h in ranges[0]]
    reciprocal = reciprocal_vectors(geometry.cell)
    scattering = []
    for h in indices:
        q = tuple(sum(h[i] * reciprocal[i][j] for i in range(3)) for j in range(3))
        s = math.sqrt(_dot(q, q)) / 2
        if coefficients is not None:
            factors = [scattering_factor(s, *coefficients[z]) for z, *_ in geometry.sites]
        else:
            factors = [scattering_lengths[z] for z, *_ in geometry.sites]
        scattering.append(factors)
    reflections = geometry.atomic_structure_factors(indices, scattering, occupancies, b_factors)
    result = geometry.fourier_synthesis(reflections)
    unit = "e/A^3" if coefficients is not None else "scattering-length/A^3"
    return Grid(result.shape, result.cell, result.values, result.origin, True, unit).with_sites(geometry.sites)


def miller_section(grid, hkl, offset=0.0, width=5.0, height=5.0, shape=(100, 100)):
    """Section centered on h*x+k*y+l*z=offset (fractional cell coordinates).

    Width/height are Angstrom. Finite inputs must contain the whole section.
    The orthonormal in-plane orientation is chosen deterministically.
    """
    hkl = tuple(hkl)
    if len(hkl) != 3 or any(int(v) != v for v in hkl) or not any(hkl):
        raise ValueError("Miller indices must be a nonzero integer triple")
    if not all(math.isfinite(v) for v in (offset, width, height)) or min(width, height) <= 0:
        raise ValueError("Positive finite section dimensions are required")
    reciprocal = reciprocal_vectors(grid.cell)
    normal = tuple(sum(hkl[i] * reciprocal[i][j] for i in range(3)) for j in range(3))
    norm2 = _dot(normal, normal)
    center = tuple(grid.origin[i] + normal[i] * offset / norm2 for i in range(3))
    normal = tuple(v / math.sqrt(norm2) for v in normal)
    seed = (1, 0, 0) if abs(normal[0]) < 0.8 else (0, 1, 0)
    u = _cross(normal, seed)
    u = tuple(v / math.sqrt(_dot(u, u)) for v in u)
    v = _cross(normal, u)
    u, v = tuple(x * width for x in u), tuple(x * height for x in v)
    origin = tuple(center[i] - (u[i] + v[i]) / 2 for i in range(3))
    return grid.section(origin, u, v, shape)
