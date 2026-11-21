//! Probe-sphere surface meshes of solids (the C++ science/SurfaceMesh, as
//! OVITO's Construct Surface Mesh): the region no probe sphere reaches, its
//! isosurface and area, the solid and empty regions and the surface atoms.

use af_core::format::format_g;
use af_core::glm::{determinant_d, inverse_d, mul_d};
use af_core::json::Json;
use af_electronic::{isosurface, Grid};
use glam::{DMat3, DVec3, IVec3};

use crate::analysis::R;
use crate::bonds_regions::ObjWriter;
use crate::core::*;
use crate::data::Parameters;

struct Lattice {
    n: [i32; 3],
    periodic: bool,
    origin: DVec3,
    /// Columns: grid step vectors.
    steps: DMat3,
    inverse: DMat3,
}

impl Lattice {
    fn size(&self) -> usize {
        self.n[0] as usize * self.n[1] as usize * self.n[2] as usize
    }
    fn index(&self, x: i32, y: i32, z: i32) -> usize {
        (z as usize * self.n[1] as usize + y as usize) * self.n[0] as usize + x as usize
    }
    fn position(&self, i: IVec3) -> DVec3 {
        self.origin + mul_d(&self.steps, i.as_dvec3())
    }
    fn wrap(&self, i: &mut IVec3) -> bool {
        for d in 0..3 {
            if self.periodic {
                i[d] = ((i[d] % self.n[d]) + self.n[d]) % self.n[d];
            } else if i[d] < 0 || i[d] >= self.n[d] {
                return false;
            }
        }
        true
    }
}

fn length(v: DVec3) -> f64 {
    v.dot(v).sqrt()
}

fn lround(x: f64) -> i32 {
    x.round() as i64 as i32
}

fn stencil(g: &Lattice, radius: f64) -> Vec<IVec3> {
    let mut reach = [0i32; 3];
    for d in 0..3 {
        let row = DVec3::new(g.inverse.col(0)[d], g.inverse.col(1)[d], g.inverse.col(2)[d]);
        reach[d] = (radius * length(row)).ceil() as i32 + 1;
    }
    let mut out = Vec::new();
    for x in -reach[0]..=reach[0] {
        for y in -reach[1]..=reach[1] {
            for z in -reach[2]..=reach[2] {
                if length(mul_d(&g.steps, DVec3::new(x as f64, y as f64, z as f64))) <= radius {
                    out.push(IVec3::new(x, y, z));
                }
            }
        }
    }
    out
}

fn stamp(g: &Lattice, sources: &[DVec3], radius: f64, distance: &mut [f64]) {
    let diagonal = 0.5 * (length(g.steps.x_axis) + length(g.steps.y_axis) + length(g.steps.z_axis));
    let offsets = stencil(g, radius + diagonal);
    for source in sources {
        let u = mul_d(&g.inverse, *source - g.origin);
        let centre = IVec3::new(lround(u.x), lround(u.y), lround(u.z));
        for o in &offsets {
            let unwrapped = centre + *o;
            let d = length(g.position(unwrapped) - *source);
            if d > radius {
                continue;
            }
            let mut i = unwrapped;
            if !g.wrap(&mut i) {
                continue;
            }
            let slot = &mut distance[g.index(i.x, i.y, i.z)];
            if d < *slot {
                *slot = d;
            }
        }
    }
}

const NEIGHBOURS: [IVec3; 6] = [IVec3::new(1, 0, 0), IVec3::new(-1, 0, 0), IVec3::new(0, 1, 0), IVec3::new(0, -1, 0), IVec3::new(0, 0, 1), IVec3::new(0, 0, -1)];

/// Sizes of the 6-connected components of the inside voxels, largest first.
fn components(g: &Lattice, inside: &[bool]) -> Vec<usize> {
    let mut label = vec![-1i32; g.size()];
    let mut sizes = Vec::new();
    let mut stack: Vec<IVec3> = Vec::new();
    for z in 0..g.n[2] {
        for y in 0..g.n[1] {
            for x in 0..g.n[0] {
                let start = g.index(x, y, z);
                if !inside[start] || label[start] >= 0 {
                    continue;
                }
                let id = sizes.len() as i32;
                let mut count = 0usize;
                label[start] = id;
                stack.push(IVec3::new(x, y, z));
                while let Some(v) = stack.pop() {
                    count += 1;
                    for s in &NEIGHBOURS {
                        let mut w = v + *s;
                        if !g.wrap(&mut w) {
                            continue;
                        }
                        let k = g.index(w.x, w.y, w.z);
                        if inside[k] && label[k] < 0 {
                            label[k] = id;
                            stack.push(w);
                        }
                    }
                }
                sizes.push(count);
            }
        }
    }
    sizes.sort_by(|a, b| b.cmp(a));
    sizes
}

pub fn surface_mesh(p: &Parameters) -> R {
    let all = points(p.array("positions")?, "positions")?;
    let pbc = p.pbc("pbc", [false; 3])?;
    let periodic = pbc.iter().any(|&b| b);
    if periodic && !(pbc[0] && pbc[1] && pbc[2]) {
        return Err("Use periodic boundaries along all three cell vectors (a slab: add vacuum) or none".into());
    }
    if periodic && !p.has("cell") {
        return Err("Periodic analysis requires a cell".into());
    }
    let probe = positive(p.number_or("probe_radius_A", 4.0)?, "probe_radius_A")?;
    let spacing = positive(p.number_or("grid_spacing_A", 0.5)?, "grid_spacing_A")?;
    let mut mask = vec![1.0; all.len()];
    if p.has("mask") {
        mask = finite_array(p.array("mask")?, "mask", 1)?.values;
        if mask.len() != all.len() {
            return Err("mask needs one value per atom".into());
        }
    }
    let atoms: Vec<DVec3> = all.iter().zip(&mask).filter(|(_, &m)| m != 0.0).map(|(x, _)| DVec3::new(x[0], x[1], x[2])).collect();
    if atoms.is_empty() {
        return Err("No atoms define the solid".into());
    }
    let mut g = Lattice { n: [0; 3], periodic, origin: DVec3::ZERO, steps: DMat3::IDENTITY, inverse: DMat3::IDENTITY };
    if periodic {
        let cell = cell_matrix(p.array("cell")?)?;
        let mut columns = [DVec3::ZERO; 3];
        for d in 0..3 {
            let lattice = DVec3::from(cell[d]);
            g.n[d] = 4.max((length(lattice) / spacing).ceil() as i32);
            columns[d] = lattice / g.n[d] as f64;
        }
        g.steps = DMat3::from_cols(columns[0], columns[1], columns[2]);
    } else {
        let (mut low, mut high) = (atoms[0], atoms[0]);
        for x in &atoms {
            for k in 0..3 {
                if x[k] < low[k] {
                    low[k] = x[k];
                }
                if high[k] < x[k] {
                    high[k] = x[k];
                }
            }
        }
        let pad = probe + 3.0 * spacing;
        g.origin = low - DVec3::splat(pad);
        let extent = high - low + DVec3::splat(2.0 * pad);
        let mut columns = [DVec3::ZERO; 3];
        for d in 0..3 {
            g.n[d] = (extent[d] / spacing).ceil() as i32 + 1;
            columns[d][d] = extent[d] / (g.n[d] - 1) as f64;
        }
        g.steps = DMat3::from_cols(columns[0], columns[1], columns[2]);
    }
    if g.size() > 60_000_000 {
        return Err("The grid is too fine for this structure; increase grid_spacing_A".into());
    }
    g.inverse = inverse_d(&g.steps);
    let voxel = determinant_d(&g.steps).abs();
    let mut to_atom = vec![f64::INFINITY; g.size()];
    stamp(&g, &atoms, probe, &mut to_atom);
    let allowed: Vec<bool> = to_atom.iter().map(|&d| d >= probe).collect();
    let mut boundary = Vec::new();
    for z in 0..g.n[2] {
        for y in 0..g.n[1] {
            for x in 0..g.n[0] {
                if !allowed[g.index(x, y, z)] {
                    continue;
                }
                for s in &NEIGHBOURS {
                    let mut w = IVec3::new(x, y, z) + *s;
                    if g.wrap(&mut w) && !allowed[g.index(w.x, w.y, w.z)] {
                        boundary.push(g.position(IVec3::new(x, y, z)));
                        break;
                    }
                }
            }
        }
    }
    let mut to_probe: Vec<f64> = allowed.iter().map(|&a| if a { 0.0 } else { f64::INFINITY }).collect();
    stamp(&g, &boundary, probe + 2.0 * spacing, &mut to_probe);
    let mut field = Grid { shape: g.n, periodic, origin: g.origin, ..Default::default() };
    let scale_of = |d: usize| (if periodic { g.n[d] } else { g.n[d] - 1 }) as f64;
    field.cell = DMat3::from_cols(g.steps.x_axis * scale_of(0), g.steps.y_axis * scale_of(1), g.steps.z_axis * scale_of(2));
    field.values = vec![0.0; g.size()];
    let (mut solid, mut empty) = (vec![false; g.size()], vec![false; g.size()]);
    let mut solid_count = 0usize;
    for k in 0..g.size() {
        let capped = if 2.0 * probe < to_probe[k] { 2.0 * probe } else { to_probe[k] };
        field.values[k] = capped - probe;
        solid[k] = field.values[k] > 0.0;
        empty[k] = !solid[k];
        if solid[k] {
            solid_count += 1;
        }
    }
    let mesh = isosurface(&field, 0.0, None)?;
    let mut area = 0.0;
    let mut t = 0;
    while t + 2 < mesh.vertices.len() {
        area += 0.5 * length((mesh.vertices[t + 1] - mesh.vertices[t]).cross(mesh.vertices[t + 2] - mesh.vertices[t]));
        t += 3;
    }
    let solid_regions = components(&g, &solid);
    let empty_regions = components(&g, &empty);
    let mut surface_atom = vec![false; all.len()];
    let offsets = stencil(&g, probe + 3.0 * spacing);
    for i in 0..all.len() {
        if mask[i] == 0.0 {
            continue;
        }
        let x = DVec3::new(all[i][0], all[i][1], all[i][2]);
        let u = mul_d(&g.inverse, x - g.origin);
        let centre = IVec3::new(lround(u.x), lround(u.y), lround(u.z));
        let mut nearest = f64::INFINITY;
        for o in &offsets {
            let mut w = centre + *o;
            let at = g.position(w);
            if !g.wrap(&mut w) || !allowed[g.index(w.x, w.y, w.z)] {
                continue;
            }
            let d = length(at - x);
            if d < nearest {
                nearest = d;
            }
        }
        surface_atom[i] = nearest <= probe + 1.5 * spacing;
    }
    let surface_count = surface_atom.iter().filter(|&&s| s).count();
    let total = voxel * g.size() as f64;
    let mut result = Json::object();
    result.set("surface_area_A2", area);
    result.set("solid_volume_A3", voxel * solid_count as f64);
    if periodic {
        result.set("cell_volume_A3", total);
        result.set("solid_fraction", solid_count as f64 / g.size() as f64);
        result.set("specific_surface_area_per_A", area / total);
    }
    result.set("solid_region_volumes_A3", Json::Array(solid_regions.iter().map(|&c| Json::Number(voxel * c as f64)).collect()));
    result.set("empty_region_volumes_A3", Json::Array(empty_regions.iter().map(|&c| Json::Number(voxel * c as f64)).collect()));
    result.set("solid_regions", solid_regions.len() as f64);
    result.set("empty_regions", empty_regions.len() as f64);
    result.set("surface_atom", Json::Array(surface_atom.iter().map(|&b| Json::Bool(b)).collect()));
    result.set("surface_atoms", surface_count as f64);
    result.set("triangles", (mesh.vertices.len() / 3) as f64);
    result.set("grid", Json::Array(g.n.iter().map(|&n| Json::Number(n as f64)).collect()));
    result.set("probe_radius_A", probe);
    if p.boolean("export_mesh", true)? && mesh.vertices.len() <= 3_000_000 {
        let mut obj = ObjWriter::new(&format!("AtomForge surface mesh: probe radius {} A, area {} A^2", format_g(probe, 6), format_g(area, 6)));
        obj.add("", &mesh.vertices);
        let mut files = Json::object();
        files.set("surface.obj", obj.text());
        result.set("files", files);
    }
    Ok(result)
}
