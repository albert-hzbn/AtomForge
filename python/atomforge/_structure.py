"""Core data types: Atom and Structure."""

import copy
import math
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class Atom:
    """A single atom with Cartesian coordinates (Angstrom) and display colour."""

    symbol: str
    x: float
    y: float
    z: float
    r: float = field(default=1.0, repr=False)
    g: float = field(default=1.0, repr=False)
    b: float = field(default=1.0, repr=False)

    def __repr__(self) -> str:
        return f"Atom({self.symbol!r}, x={self.x:.4f}, y={self.y:.4f}, z={self.z:.4f})"

    def copy(self) -> "Atom":
        return copy.copy(self)


class Structure:
    """
    A collection of atoms with an optional periodic unit cell.

    Attributes
    ----------
    atoms : list[Atom]
        All atoms in Cartesian coordinates (Angstrom).
    cell : list[list[float]] | None
        3×3 lattice vectors ``[[a1,a2,a3], [b1,b2,b3], [c1,c2,c3]]``
        in Angstrom, or ``None`` for a non-periodic (cluster) structure.
    """

    def __init__(self) -> None:
        self.atoms: List[Atom] = []
        self.cell: Optional[List[List[float]]] = None

    # ── Dunder helpers ──────────────────────────────────────────────────────────

    def __len__(self) -> int:
        return len(self.atoms)

    def __repr__(self) -> str:
        periodicity = "periodic" if self.cell else "non-periodic"
        species: dict = {}
        for a in self.atoms:
            species[a.symbol] = species.get(a.symbol, 0) + 1
        comp = " ".join(f"{s}{n}" for s, n in sorted(species.items()))
        return f"Structure({len(self.atoms)} atoms, {comp}, {periodicity})"

    def copy(self) -> "Structure":
        s = Structure()
        s.atoms = [a.copy() for a in self.atoms]
        s.cell = [list(row) for row in self.cell] if self.cell else None
        return s

    # ── Atom editing ────────────────────────────────────────────────────────────

    def add_atom(self, symbol: str, x: float, y: float, z: float,
                 r: float = -1.0, g: float = -1.0, b: float = -1.0) -> Atom:
        """Append an atom and return it.  Default colour is the element's CPK colour."""
        from ._io import _default_color
        cr, cg, cb = _default_color(symbol)
        atom = Atom(symbol, float(x), float(y), float(z),
                    r=cr if r < 0 else r,
                    g=cg if g < 0 else g,
                    b=cb if b < 0 else b)
        self.atoms.append(atom)
        return atom

    def remove_atom(self, index: int) -> Atom:
        """Remove and return the atom at *index*."""
        return self.atoms.pop(index)

    # ── Bulk transforms ─────────────────────────────────────────────────────────

    def translate(self, dx: float, dy: float, dz: float) -> "Structure":
        """Shift all atoms in-place by (dx, dy, dz) Angstrom.  Returns self."""
        for a in self.atoms:
            a.x += dx
            a.y += dy
            a.z += dz
        return self

    def scale(self, factor: float) -> "Structure":
        """Scale all Cartesian coordinates and cell vectors by *factor*.  Returns self."""
        for a in self.atoms:
            a.x *= factor
            a.y *= factor
            a.z *= factor
        if self.cell:
            self.cell = [[v * factor for v in row] for row in self.cell]
        return self

    def repeat(self, nx: int, ny: int = 1, nz: int = 1) -> "Structure":
        """Return a new supercell replicated *nx × ny × nz* times.  Requires a unit cell."""
        if not self.cell:
            raise ValueError("repeat() requires a unit cell (set structure.cell first)")
        a, b, c = self.cell
        out = Structure()
        out.cell = [
            [a[i] * nx for i in range(3)],
            [b[i] * ny for i in range(3)],
            [c[i] * nz for i in range(3)],
        ]
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for atom in self.atoms:
                        new_atom = atom.copy()
                        new_atom.x += a[0]*ix + b[0]*iy + c[0]*iz
                        new_atom.y += a[1]*ix + b[1]*iy + c[1]*iz
                        new_atom.z += a[2]*ix + b[2]*iy + c[2]*iz
                        out.atoms.append(new_atom)
        return out

    def filter_species(self, *symbols: str) -> "Structure":
        """Return a new Structure containing only atoms of the given element symbols."""
        keep = set(symbols)
        out = self.copy()
        out.atoms = [a for a in out.atoms if a.symbol in keep]
        return out

    def set_cell(self, a: float, b: float, c: float,
                 alpha: float = 90.0, beta: float = 90.0,
                 gamma: float = 90.0) -> "Structure":
        """
        Set the unit cell from lattice parameters (lengths in Å, angles in degrees).
        Uses the standard convention: **a** along x, **b** in the xy-plane.
        Returns self.
        """
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
        cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
        self.cell = [[ax, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]]
        return self

    # ── Viewing / I/O ───────────────────────────────────────────────────────────

    def view(self) -> None:
        """Open this structure in the AtomForge GUI (non-blocking)."""
        from ._viewer import view
        view(self)

    def save(self, path: str) -> None:
        """Save to file; format is inferred from the file extension."""
        from ._io import save
        save(self, path)
