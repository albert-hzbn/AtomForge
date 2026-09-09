"""Dependency-free unit-cell and coordinate transformations."""
from __future__ import annotations

import math
from typing import List, Tuple


def _params_from_cell(cell: List[List[float]]) -> Tuple[float, float, float, float, float, float]:
    """Return lengths (Angstrom) and alpha/beta/gamma angles (degrees)."""
    a, b, c = (math.sqrt(sum(v*v for v in row)) for row in cell)

    def angle(first: int, second: int, denominator: float) -> float:
        cosine = (sum(x*y for x, y in zip(cell[first], cell[second])) / denominator
                  if denominator else 0.0)
        return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

    return a, b, c, angle(1, 2, b*c), angle(0, 2, a*c), angle(0, 1, a*b)


def _cell_from_params(a: float, b: float, c: float,
                      alpha: float, beta: float, gamma: float
                      ) -> List[List[float]]:
    """Build a 3×3 cell matrix from lattice parameters (lengths Å, angles °)."""
    values = (a, b, c, alpha, beta, gamma)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("cell lengths and angles must be finite")
    if min(a, b, c) <= 0:
        raise ValueError("cell lengths must be positive")
    if not all(0.0 < angle < 180.0 for angle in (alpha, beta, gamma)):
        raise ValueError("cell angles must be between 0 and 180 degrees")

    rad = math.pi / 180.0
    ca = math.cos(alpha * rad)
    cb = math.cos(beta  * rad)
    cg = math.cos(gamma * rad)
    sg = math.sin(gamma * rad)
    ax = a
    bx = b * cg
    by = b * sg
    cx = c * cb
    cy = c * (ca - cb * cg) / sg if sg > 1e-10 else 0.0
    cz_sq = c * c - cx * cx - cy * cy
    if cz_sq <= 1e-12:
        raise ValueError("cell parameters do not define a non-degenerate cell")
    cz = math.sqrt(cz_sq)
    return [[ax, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]]


def _frac_to_cart(fx: float, fy: float, fz: float,
                  cell: List[List[float]]) -> Tuple[float, float, float]:
    a, b, c = cell
    x = fx * a[0] + fy * b[0] + fz * c[0]
    y = fx * a[1] + fy * b[1] + fz * c[1]
    z = fx * a[2] + fy * b[2] + fz * c[2]
    return x, y, z


def _cart_to_frac(x: float, y: float, z: float,
                  cell: List[List[float]]) -> Tuple[float, float, float]:
    """Invert a (generally non-orthogonal) cell matrix to get fractional coords."""
    a, b, c = cell
    # Build 3×3 and invert
    m = [list(col) for col in zip(a, b, c)]
    # Determinant
    det = (m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
         - m[0][1] * (m[1][0]*m[2][2] - m[1][2]*m[2][0])
         + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]))
    if abs(det) < 1e-12:
        raise ValueError("cannot convert coordinates with a singular cell")
    inv = [
        [(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det,
         (m[0][2]*m[2][1]-m[0][1]*m[2][2])/det,
         (m[0][1]*m[1][2]-m[0][2]*m[1][1])/det],
        [(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det,
         (m[0][0]*m[2][2]-m[0][2]*m[2][0])/det,
         (m[0][2]*m[1][0]-m[0][0]*m[1][2])/det],
        [(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det,
         (m[0][1]*m[2][0]-m[0][0]*m[2][1])/det,
         (m[0][0]*m[1][1]-m[0][1]*m[1][0])/det],
    ]
    fx = inv[0][0]*x + inv[0][1]*y + inv[0][2]*z
    fy = inv[1][0]*x + inv[1][1]*y + inv[1][2]*z
    fz = inv[2][0]*x + inv[2][1]*y + inv[2][2]*z
    return fx, fy, fz
