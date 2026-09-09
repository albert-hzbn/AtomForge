"""Check published examples and the manual's executable Python recipes."""
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
MANUAL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
os.environ.setdefault("ATOMFORGE_PATH", str(ROOT / "build" / ("AtomForge.exe" if os.name == "nt" else "AtomForge")))
import atomforge as af
from atomforge.electronic import Grid, load_volume, miller_section


def main():
    examples = MANUAL / "examples"
    scratch = MANUAL / "build" / "recipe-checks"
    scratch.mkdir(parents=True, exist_ok=True)
    records = json.loads((examples / "builder-results.json").read_text())
    for record in records:
        assert len(af.load(str(examples / record["file"]))) == record["atoms"]
    shapes = json.loads((examples / "shape-results.json").read_text())
    assert len({record["shape"] for record in shapes}) == 7
    for record in shapes:
        structure = af.load(str(examples / record["file"]))
        assert len(structure) == record["atoms"] > 0
        assert all(atom.symbol == "Cu" for atom in structure.atoms)
        assert all(math.isfinite(value) for atom in structure.atoms
                   for value in (atom.x, atom.y, atom.z))
    cu = af.load(str(examples / "cu_fcc.cif"))
    assert len(cu.repeat(5, 5, 5)) == 500
    host = af.load(str(examples / "cu_host.cif"))
    dislocation = af.load(str(examples / "cu_dislocation.cif"))
    assert len(dislocation) == len(host) == 500
    assert dislocation.cell is not None
    assert all(math.isfinite(value) for atom in dislocation.atoms for value in (atom.x, atom.y, atom.z))
    assert any(math.dist((a.x, a.y, a.z), (b.x, b.y, b.z)) > 0.01
               for a, b in zip(host.atoms, dislocation.atoms))
    copy = host.copy().translate(0.1, 0, 0).scale(1.01)
    copy.save(str(scratch / "edited.cif"))
    assert len(af.load(str(scratch / "edited.cif"))) == 500
    alloy = af.load(str(examples / "cu_ni_alloy.cif"))
    assert len(alloy.filter_species("Cu")) == 350
    assert len(alloy.filter_species("Ni")) == 150

    # Analytic periodic fixture: mean density 1 in a 4^3 A^3 cell.
    n = 8
    values = [1 + 0.2 * math.cos(2*math.pi*x/n)
              for z in range(n) for y in range(n) for x in range(n)]
    grid = Grid((n, n, n), ((4, 0, 0), (0, 4, 0), (0, 0, 4)),
                values, periodic=True, unit="e/A^3")
    grid = grid.with_sites([(1, 0, 0, 0), (1, 2, 2, 2)])
    assert math.isclose(grid.integrate(), 64, abs_tol=1e-10)
    reference = grid * 0.5
    delta = grid.density_difference(reference, reference)
    assert max(map(abs, delta.values)) == 0
    assert delta.charge_summary() == {"accumulation": 0, "depletion": 0, "net": 0}
    mask = grid.threshold_mask(1.0, 1.2)
    complement = mask.invert_mask()
    assert math.isclose(grid.apply_mask(mask).integrate() + grid.apply_mask(complement).integrate(), 64, abs_tol=1e-10)
    assert math.isclose(mask.boolean(complement, "union").integrate(), 64, abs_tol=1e-10)
    assert mask.boolean(complement, "intersection").integrate() == 0
    assert math.isclose(grid.cumulative_charge(2)[-1][1], 64, abs_tol=1e-10)
    for derived in [*grid.gradient(), grid.laplacian(), grid.smooth(.2, 2), *grid.energy_density()]:
        assert all(map(math.isfinite, derived.values))
    assert len(grid.line_profile((0, 0, 0), (1, 1, 1), 100)) == 100
    assert len(grid.planar_average(2)) == n
    assert len(grid.macroscopic_average(2, 5)) == n
    assert grid.integrate_sphere((1, 1, 1), .5) > 0
    assert math.isclose(sum(row[0] for row in grid.voronoi_integrate()), 64, abs_tol=1e-10)
    assert len(grid.peaks(limit=10)) <= 10
    section = grid.section((0, 0, 2), grid.cell[0], grid.cell[1], (20, 20))
    assert len(section.section_values) == 20
    assert section.contours(1.0)
    assert miller_section(grid, (1, 1, 1), 1.5, 2, 2, (20, 20)).section_values
    refs = grid.structure_factors([(0, 0, 0), (1, 0, 0), (-1, 0, 0)])
    rebuilt = grid.fourier_synthesis(refs)
    assert max(abs(a-b) for a, b in zip(grid.values, rebuilt.values)) < 1e-10
    assert all(map(math.isfinite, grid.patterson().values))
    electrostatic = grid.ewald([1, -1], .8, 12, 8)
    assert math.isfinite(electrostatic["energy"])
    assert len(electrostatic["potentials"]) == 2
    for extension in ("cube", "xsf"):
        path = scratch / ("density." + extension)
        grid.save(str(path))
        read = load_volume(str(path), quantity="density").fields[0].as_periodic()
        assert read.shape == grid.shape
        assert math.isclose(read.integrate(), 64, abs_tol=1e-8)
    for extension in ("obj", "ply"):
        path = scratch / ("surface." + extension)
        grid.isosurface(1.1).save(str(path))
        assert path.stat().st_size > 0
    print(f"Validated {len(records)} builder outputs, {len(shapes)} shape examples, "
          "host/alloy/dislocation checks, and Python numerical recipes.")


if __name__ == "__main__":
    main()
