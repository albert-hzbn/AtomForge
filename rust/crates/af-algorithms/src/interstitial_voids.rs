//! Interstitial void detection (the C++ algorithms/InterstitialVoidAnalysis):
//! largest empty spheres from circumspheres of neighbouring atoms, bond
//! midpoints and clearance peaks on a grid, refined to the centres of their
//! first-neighbour shells and classified by coordination. Single-precision
//! glm arithmetic and libstdc++ ordering throughout, so the sites (and their
//! order) equal the C++ ones.

use std::collections::{BTreeSet, HashMap};

use af_core::cxx::{nth_element_by, sort_by};
use af_core::glm::{self, cross, dot, length, mul, normalize};
use af_core::Structure;
use glam::{IVec3, Mat3, Vec3};

const PLANE_TOL: f32 = 2e-3;
const DEGENERATE_TOL: f32 = 1e-8;
const EMPTY_SPHERE_TOL: f32 = 1e-4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VoidKind {
    Tetrahedral = 0,
    Octahedral,
    Irregular,
}

#[derive(Clone, Debug, PartialEq)]
pub struct VoidRegion {
    pub position: Vec3,
    pub clearance: f32,
    pub coordination: i32,
    pub volume: f32,
    pub volume_bin: i32,
    pub kind: VoidKind,
    pub neighbor_atom_indices: Vec<i32>,
    pub faces: Vec<[i32; 3]>,
    pub face_vertices: Vec<[Vec3; 3]>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VoidParams {
    pub grid_resolution: i32,
    /// 0 or negative: every unique void.
    pub max_voids: i32,
    pub min_clearance: f32,
    pub min_separation: f32,
    pub first_neighbor_shell_factor: f32,
    pub max_neighbor_count: i32,
}

impl Default for VoidParams {
    fn default() -> Self {
        Self { grid_resolution: 14, max_voids: 0, min_clearance: 0.6, min_separation: 0.5, first_neighbor_shell_factor: 1.22, max_neighbor_count: 12 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VoidResult {
    pub success: bool,
    pub message: String,
    pub regions: Vec<VoidRegion>,
}

/// std::min / std::max on floats (the first argument wins ties and NaN).
fn min_f(a: f32, b: f32) -> f32 {
    if b < a {
        b
    } else {
        a
    }
}
fn max_f(a: f32, b: f32) -> f32 {
    if a < b {
        b
    } else {
        a
    }
}

#[derive(Clone, Copy)]
struct ShellPoint {
    atom_index: i32,
    position: Vec3,
}

#[derive(Clone, Copy)]
struct LocalNeighbor {
    position: Vec3,
    d2: f32,
}

struct Geometry {
    use_pbc: bool,
    origin: Vec3,
    cell: Mat3,
    inv: Mat3,
    exact: bool,
}

impl Geometry {
    fn wrap(&self, p: Vec3) -> Vec3 {
        if !self.use_pbc {
            return p;
        }
        let mut frac = mul(&self.inv, p - self.origin);
        frac -= glm::floor(frac);
        self.origin + mul(&self.cell, frac)
    }

    /// The shortest periodic image of a separation (also trying the 26
    /// neighbouring images when `exact`).
    fn nearest_image(&self, d: Vec3) -> Vec3 {
        let base = glm::minimum_image_delta(d, self.use_pbc, &self.cell, &self.inv);
        if !self.use_pbc || !self.exact {
            return base;
        }
        let mut best = base;
        let mut best_length = dot(base, base);
        for a in -1..=1 {
            for b in -1..=1 {
                for c in -1..=1 {
                    let candidate = base + mul(&self.cell, Vec3::new(a as f32, b as f32, c as f32));
                    let l = dot(candidate, candidate);
                    if l < best_length {
                        best = candidate;
                        best_length = l;
                    }
                }
            }
        }
        best
    }

    /// The smallest plane spacing of the cell.
    fn width(&self) -> f32 {
        let mut width = f32::MAX;
        for k in 0..3 {
            let row = Vec3::new(self.inv.col(0)[k], self.inv.col(1)[k], self.inv.col(2)[k]);
            let l = length(row);
            if l > 0.0 {
                width = min_f(width, 1.0 / l);
            }
        }
        width
    }
}

fn compute_circumsphere(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Option<(Vec3, f32)> {
    let ba = b - a;
    let ca = c - a;
    let da = d - a;
    let m = Mat3::from_cols(2.0 * ba, 2.0 * ca, 2.0 * da);
    let mt = m.transpose();
    let det = glm::determinant(&mt);
    if det.abs() < 1e-9 {
        return None;
    }
    let rhs = Vec3::new(dot(b, b) - dot(a, a), dot(c, c) - dot(a, a), dot(d, d) - dot(a, a));
    let center = mul(&glm::inverse(&mt), rhs);
    let radius = length(center - a);
    radius.is_finite().then_some((center, radius))
}

fn build_first_neighbor_shell(center: Vec3, atoms: &[Vec3], g: &Geometry, shell_factor: f32, max_neighbor_count: i32) -> Vec<ShellPoint> {
    #[derive(Clone, Copy)]
    struct DistAtom {
        d2: f32,
        atom_index: i32,
        pos: Vec3,
    }
    let mut near: Vec<DistAtom> = Vec::with_capacity(atoms.len());
    let mut min_d2 = f32::MAX;
    for (ai, &p) in atoms.iter().enumerate() {
        let delta = g.nearest_image(p - center);
        let d2 = dot(delta, delta);
        min_d2 = min_f(min_d2, d2);
        near.push(DistAtom { d2, atom_index: ai as i32, pos: center + delta });
    }
    let mut shell = Vec::new();
    if near.is_empty() {
        return shell;
    }
    let factor = max_f(1.02, shell_factor);
    let limit2 = max_f(0.0, min_d2) * factor * factor;
    // In a cell narrower than the shell one atom can sit in it more than once.
    if g.use_pbc && 2.0 * limit2.sqrt() > g.width() {
        let count = near.len();
        for n in 0..count {
            for a in -1..=1 {
                for b in -1..=1 {
                    for c in -1..=1 {
                        if a == 0 && b == 0 && c == 0 {
                            continue;
                        }
                        let position = near[n].pos + mul(&g.cell, Vec3::new(a as f32, b as f32, c as f32));
                        let delta = position - center;
                        let d2 = dot(delta, delta);
                        if d2 <= limit2 {
                            near.push(DistAtom { d2, atom_index: near[n].atom_index, pos: position });
                        }
                    }
                }
            }
        }
    }
    let mut candidates: Vec<DistAtom> = near.iter().copied().filter(|da| da.d2 <= limit2).collect();
    if candidates.len() as i32 > max_neighbor_count {
        let m = max_neighbor_count as usize;
        nth_element_by(&mut candidates, m, |a, b| a.d2 < b.d2);
        candidates.truncate(m);
        sort_by(&mut candidates, |a, b| a.d2 < b.d2);
    }
    shell.extend(candidates.iter().map(|da| ShellPoint { atom_index: da.atom_index, position: da.pos }));
    // The nearest atom is always present in pathological tolerance cases.
    if shell.is_empty() {
        let mut best = 0;
        for i in 1..near.len() {
            if near[i].d2 < near[best].d2 {
                best = i;
            }
        }
        shell.push(ShellPoint { atom_index: near[best].atom_index, position: near[best].pos });
    }
    shell
}

/// Whether the origin lies in the convex hull of the unit vectors
/// (Frank-Wolfe iterations towards the hull's minimum-norm point).
fn origin_in_hull(directions: &[Vec3]) -> bool {
    if directions.len() < 2 {
        return false;
    }
    let mut x = directions[0];
    for _ in 0..500 {
        if dot(x, x) < 1e-6 {
            return true;
        }
        let mut best = 0;
        for i in 1..directions.len() {
            if dot(x, directions[i]) < dot(x, directions[best]) {
                best = i;
            }
        }
        let d = directions[best] - x;
        let gap = -dot(x, d);
        if gap < 1e-7 {
            break;
        }
        let step = min_f(1.0, gap / max_f(1e-12, dot(d, d)));
        x += step * d;
    }
    dot(x, x) < 1e-4
}

struct PlaneGroup {
    normal: Vec3,
    d: f32,
    local_vertices: Vec<usize>,
}

fn push_unique(out: &mut Vec<usize>, value: usize) {
    if !out.contains(&value) {
        out.push(value);
    }
}

/// Faces of the convex polyhedron of the first-neighbour shell and its volume.
fn polyhedral_faces(center: Vec3, shell_points: &[ShellPoint]) -> Option<(Vec<[i32; 3]>, Vec<[Vec3; 3]>, f32)> {
    let n = shell_points.len();
    if n < 4 {
        return None;
    }
    let shell: Vec<Vec3> = shell_points.iter().map(|p| p.position).collect();
    let mut planes: Vec<PlaneGroup> = Vec::with_capacity(n * 2);
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let a = shell[i];
                let mut normal = cross(shell[j] - a, shell[k] - a);
                let len2 = dot(normal, normal);
                if len2 < DEGENERATE_TOL {
                    continue;
                }
                if dot(normal, center - a) > 0.0 {
                    normal = -normal;
                }
                let supports = (0..n).all(|t| t == i || t == j || t == k || !(dot(normal, shell[t] - a) > PLANE_TOL));
                if !supports {
                    continue;
                }
                let nrm = normal / len2.sqrt();
                let d = dot(nrm, a);
                let group = planes.iter().position(|p| !((dot(p.normal, nrm) - 1.0).abs() > 1e-4) && !((p.d - d).abs() > 1e-3));
                match group {
                    None => planes.push(PlaneGroup { normal: nrm, d, local_vertices: vec![i, j, k] }),
                    Some(g) => {
                        push_unique(&mut planes[g].local_vertices, i);
                        push_unique(&mut planes[g].local_vertices, j);
                        push_unique(&mut planes[g].local_vertices, k);
                    }
                }
            }
        }
    }
    if planes.is_empty() {
        return None;
    }
    let mut faces = Vec::new();
    let mut face_vertices: Vec<[Vec3; 3]> = Vec::new();
    let mut unique: BTreeSet<[usize; 3]> = BTreeSet::new();
    for plane in &planes {
        let mut verts = plane.local_vertices.clone();
        if verts.len() < 3 {
            continue;
        }
        let mut centroid = Vec3::ZERO;
        for &li in &verts {
            centroid += shell[li];
        }
        centroid /= verts.len() as f32;
        let mut u = shell[verts[0]] - centroid;
        if dot(u, u) < DEGENERATE_TOL {
            u = cross(plane.normal, Vec3::new(1.0, 0.0, 0.0));
        }
        if dot(u, u) < DEGENERATE_TOL {
            u = cross(plane.normal, Vec3::new(0.0, 1.0, 0.0));
        }
        if dot(u, u) < DEGENERATE_TOL {
            continue;
        }
        u = normalize(u);
        let v = normalize(cross(plane.normal, u));
        sort_by(&mut verts, |&a, &b| {
            let ra = shell[a] - centroid;
            let rb = shell[b] - centroid;
            dot(ra, v).atan2(dot(ra, u)) < dot(rb, v).atan2(dot(rb, u))
        });
        let first = verts[0];
        for t in 1..verts.len().saturating_sub(1) {
            let (li0, mut li1, mut li2) = (first, verts[t], verts[t + 1]);
            if dot(cross(shell[li1] - shell[li0], shell[li2] - shell[li0]), plane.normal) < 0.0 {
                std::mem::swap(&mut li1, &mut li2);
            }
            let mut sorted_key = [li0, li1, li2];
            sorted_key.sort();
            if unique.insert(sorted_key) {
                faces.push([shell_points[li0].atom_index, shell_points[li1].atom_index, shell_points[li2].atom_index]);
                face_vertices.push([shell[li0], shell[li1], shell[li2]]);
            }
        }
    }
    if faces.is_empty() {
        return None;
    }
    let mut volume = 0.0f64;
    for tri in &face_vertices {
        let (a, b, c) = (tri[0] - center, tri[1] - center, tri[2] - center);
        volume += (dot(a, cross(b, c)) as f64).abs() / 6.0;
    }
    let volume = volume.abs() as f32;
    (volume > 1e-6).then_some((faces, face_vertices, volume))
}

fn assign_volume_bins(regions: &mut [VoidRegion]) {
    if regions.is_empty() {
        return;
    }
    let mut order: Vec<usize> = (0..regions.len()).collect();
    sort_by(&mut order, |&a, &b| regions[a].volume < regions[b].volume);
    let mut representative = -1.0f32;
    let mut current = -1;
    for idx in order {
        let v = regions[idx].volume;
        if current < 0 {
            current = 0;
            representative = v;
        } else {
            let tol = max_f(1e-4, representative * 0.05);
            if (v - representative).abs() > tol {
                current += 1;
                representative = v;
            } else {
                representative = 0.5 * (representative + v);
            }
        }
        regions[idx].volume_bin = current;
    }
}

#[derive(Clone, Copy)]
struct Candidate {
    center: Vec3,
    radius: f32,
    a: i32,
}

/// The C++ worker count: std::thread::hardware_concurrency (8 when unknown).
/// The result depends on it (the per-worker candidates are concatenated in
/// worker order), so `threads` overrides it for reproducible output.
fn worker_count(atoms: usize, threads: Option<usize>) -> usize {
    let hw = threads.unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8));
    hw.min(atoms).max(1)
}

/// Runs `work(index)` over indices interleaved across workers (worker w takes
/// w, w + W, ...) and concatenates each worker's output in worker order, as
/// the C++ per-worker vectors are.
fn interleaved<T: Send>(count: usize, workers: usize, work: impl Fn(usize, &mut Vec<T>) + Sync) -> Vec<Vec<T>> {
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|w| {
                let work = &work;
                scope.spawn(move || {
                    let mut out = Vec::new();
                    let mut i = w;
                    while i < count {
                        work(i, &mut out);
                        i += workers;
                    }
                    out
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().expect("void worker")).collect()
    })
}

pub fn detect_interstitial_voids(structure: &Structure, params: &VoidParams) -> VoidResult {
    detect_interstitial_voids_with_threads(structure, params, None)
}

/// `detect_interstitial_voids` as on a machine with `threads` hardware
/// threads (None: this machine's), for output independent of the machine.
pub fn detect_interstitial_voids_with_threads(structure: &Structure, params: &VoidParams, threads: Option<usize>) -> VoidResult {
    let mut result = VoidResult::default();
    if structure.atoms.is_empty() {
        result.message = "Source structure has no atoms.".into();
        return result;
    }
    let grid_n = params.grid_resolution.clamp(6, 48);
    let limit_voids = params.max_voids > 0;
    let max_voids = if limit_voids { params.max_voids.max(1) } else { i32::MAX };
    let min_clear = max_f(0.05, params.min_clearance);
    let min_sep = max_f(0.05, params.min_separation);
    let neighbor_shell = max_f(1.02, params.first_neighbor_shell_factor);
    let max_neighbor = params.max_neighbor_count.clamp(4, 24);

    let (use_pbc, cell, inv) = match glm::cell_matrices(structure) {
        Some((c, i)) => (true, c, i),
        None => (false, Mat3::IDENTITY, Mat3::IDENTITY),
    };
    let origin = Vec3::new(structure.cell_offset[0] as f32, structure.cell_offset[1] as f32, structure.cell_offset[2] as f32);
    let atoms: Vec<Vec3> = structure.atoms.iter().map(|a| Vec3::new(a.x as f32, a.y as f32, a.z as f32)).collect();
    let neighbor_pool = (grid_n + 10).clamp(16, 40) as usize;

    let mut g = Geometry { use_pbc, origin, cell, inv, exact: false };
    // A narrow cell: the 26 (two shells: 124) adjacent images of every atom too.
    let mut image_shifts = vec![Vec3::ZERO];
    if use_pbc {
        let width = g.width();
        let spacing = (glm::determinant(&cell).abs() / atoms.len().max(1) as f32).cbrt();
        if width < 4.4 * spacing {
            for a in -2..=2 {
                for b in -2..=2 {
                    for c in -2..=2 {
                        if a != 0 || b != 0 || c != 0 {
                            image_shifts.push(mul(&cell, Vec3::new(a as f32, b as f32, c as f32)));
                        }
                    }
                }
            }
        }
    }
    g.exact = image_shifts.len() > 1;
    let g = g;

    let atom_count = atoms.len();
    let workers = worker_count(atom_count, threads);
    let per_worker = interleaved(atom_count, workers, |ai, local: &mut Vec<Candidate>| {
        #[derive(Clone, Copy)]
        struct Near {
            first: f32,
            delta: Vec3,
        }
        let mut near: Vec<Near> = Vec::with_capacity(atom_count * image_shifts.len());
        for aj in 0..atom_count {
            let base = glm::minimum_image_delta(atoms[aj] - atoms[ai], use_pbc, &cell, &inv);
            for &shift in &image_shifts {
                let d = base + shift;
                let d2 = dot(d, d);
                if d2 < 1e-8 {
                    continue;
                }
                near.push(Near { first: d2, delta: d });
            }
        }
        if near.len() < 4 {
            return;
        }
        let pool_n = neighbor_pool.min(near.len());
        if pool_n < 4 {
            return;
        }
        if pool_n < near.len() {
            nth_element_by(&mut near, pool_n, |a, b| a.first < b.first);
            near.truncate(pool_n);
        }
        sort_by(&mut near, |a, b| a.first < b.first);
        let pool: Vec<LocalNeighbor> =
            near.iter().take(pool_n).map(|n| LocalNeighbor { position: atoms[ai] + n.delta, d2: n.first }).collect();
        let d1 = max_f(0.0, near[0].first).sqrt();
        let span2 = (d1 * 1.9) * (d1 * 1.9);
        let tetra: Vec<usize> = (0..pool.len()).filter(|&pi| pool[pi].d2 <= span2).collect();
        if tetra.len() < 3 {
            return;
        }
        for ib in 0..tetra.len() {
            for ic in ib + 1..tetra.len() {
                for id in ic + 1..tetra.len() {
                    let (b, c, d) = (tetra[ib], tetra[ic], tetra[id]);
                    let Some((cc, rr)) = compute_circumsphere(atoms[ai], pool[b].position, pool[c].position, pool[d].position) else { continue };
                    if !(rr > 0.0) || !rr.is_finite() {
                        continue;
                    }
                    let tol = rr - EMPTY_SPHERE_TOL;
                    let tol2 = tol * tol;
                    let empty = pool.iter().enumerate().all(|(pi, p)| {
                        if pi == b || pi == c || pi == d {
                            return true;
                        }
                        let dd = p.position - cc;
                        !(dot(dd, dd) < tol2)
                    });
                    if !empty {
                        continue;
                    }
                    local.push(Candidate { center: g.wrap(cc), radius: rr, a: ai as i32 });
                }
            }
        }
        // Midpoints of nearby atom pairs (edge/bridge-centred voids such as BCC octahedral sites).
        let pair_n = pool_n.min(14);
        for pb in pool.iter().take(pair_n) {
            let mid = atoms[ai] + 0.5 * (pb.position - atoms[ai]);
            let rr = max_f(0.0, pb.d2).sqrt() * 0.5;
            if rr >= min_clear {
                local.push(Candidate { center: g.wrap(mid), radius: rr, a: ai as i32 });
            }
        }
    });
    let mut centers: Vec<Candidate> = per_worker.into_iter().flatten().collect();

    // Clearance peaks on a grid recover sites the circumspheres under-sample.
    let mut bb_min = Vec3::splat(f32::MAX);
    let mut bb_max = Vec3::splat(-f32::MAX);
    for &p in &atoms {
        bb_min = glm::min(bb_min, p);
        bb_max = glm::max(bb_max, p);
    }
    let gn = grid_n as usize;
    let total = gn * gn * gn;
    let clearance_of = |p: Vec3| -> f32 {
        let mut nearest2 = f32::MAX;
        for &a in &atoms {
            let d = g.nearest_image(p - a);
            nearest2 = min_f(nearest2, dot(d, d));
        }
        max_f(0.0, nearest2).sqrt()
    };
    let sample_workers = interleaved(total, workers, |idx, out: &mut Vec<(usize, Vec3, f32)>| {
        let i = idx / (gn * gn);
        let rem = idx % (gn * gn);
        let (j, k) = (rem / gn, rem % gn);
        let frac = Vec3::new((i as f32 + 0.5) / grid_n as f32, (j as f32 + 0.5) / grid_n as f32, (k as f32 + 0.5) / grid_n as f32);
        let p = if use_pbc { origin + mul(&cell, frac) } else { bb_min + frac * (bb_max - bb_min) };
        out.push((idx, p, clearance_of(p)));
    });
    let mut samples = vec![(Vec3::ZERO, 0.0f32); total];
    for (idx, p, c) in sample_workers.into_iter().flatten() {
        samples[idx] = (p, c);
    }
    let sample_idx = |i: i32, j: i32, k: i32| (i as usize) * gn * gn + (j as usize) * gn + k as usize;
    let wrap_grid = |v: i32| {
        let r = v % grid_n;
        if r < 0 {
            r + grid_n
        } else {
            r
        }
    };
    let (begin, end) = if use_pbc { (0, grid_n) } else { (1, grid_n - 1) };
    for i in begin..end {
        for j in begin..end {
            for k in begin..end {
                let (sp, sc) = samples[sample_idx(i, j, k)];
                if sc < min_clear {
                    continue;
                }
                let mut peak = true;
                'search: for di in -1..=1 {
                    for dj in -1..=1 {
                        for dk in -1..=1 {
                            if di == 0 && dj == 0 && dk == 0 {
                                continue;
                            }
                            let (ni, nj, nk) =
                                if use_pbc { (wrap_grid(i + di), wrap_grid(j + dj), wrap_grid(k + dk)) } else { (i + di, j + dj, k + dk) };
                            if !use_pbc && (ni <= 0 || nj <= 0 || nk <= 0 || ni >= grid_n - 1 || nj >= grid_n - 1 || nk >= grid_n - 1) {
                                continue;
                            }
                            if samples[sample_idx(ni, nj, nk)].1 > sc {
                                peak = false;
                                break 'search;
                            }
                        }
                    }
                }
                if peak {
                    centers.push(Candidate { center: g.wrap(sp), radius: sc, a: -1 });
                }
            }
        }
    }

    let budget = 20000usize.max((max_voids as usize).saturating_mul(500).max(atom_count * 120));
    if centers.len() > budget {
        nth_element_by(&mut centers, budget, |a, b| a.radius > b.radius);
        centers.truncate(budget);
    }
    sort_by(&mut centers, |a, b| a.radius > b.radius);

    let shell_at = |point: Vec3| -> Vec<ShellPoint> {
        let shell = build_first_neighbor_shell(point, &atoms, &g, neighbor_shell, max_neighbor);
        if shell.len() < 4 {
            return build_first_neighbor_shell(point, &atoms, &g, max_f(1.5, neighbor_shell), 8.max(max_neighbor));
        }
        shell
    };

    // Accepted sites by bins of size minSep (fractional bins in a periodic cell).
    let mut bin_counts = [1i32; 3];
    if use_pbc {
        for (k, count) in bin_counts.iter_mut().enumerate() {
            let row = Vec3::new(inv.col(0)[k], inv.col(1)[k], inv.col(2)[k]);
            let width = if length(row) > 0.0 { 1.0 / length(row) } else { min_sep };
            *count = ((width / min_sep).floor() as i32).clamp(1, 1024);
        }
    }
    let bin_of = |p: Vec3| -> IVec3 {
        if use_pbc {
            let f = mul(&inv, p - origin);
            let mut bin = IVec3::ZERO;
            for k in 0..3 {
                let b = ((f[k] - f[k].floor()) * bin_counts[k] as f32).floor() as i32;
                bin[k] = b.max(0).min(bin_counts[k] - 1);
            }
            bin
        } else {
            let q = glm::floor(p / min_sep);
            IVec3::new(q.x as i32, q.y as i32, q.z as i32)
        }
    };
    let bin_key = |b: IVec3| -> u64 {
        (((b.x + 1048576) as u32 as u64) << 42) ^ (((b.y + 1048576) as u32 as u64) << 21) ^ ((b.z + 1048576) as u32 as u64)
    };
    let mut site_bins: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut regions: Vec<VoidRegion> = Vec::new();
    let near_site = |p: Vec3, regions: &[VoidRegion], site_bins: &HashMap<u64, Vec<usize>>| -> bool {
        let centre = bin_of(p);
        let reach: [i32; 3] = std::array::from_fn(|k| if use_pbc && bin_counts[k] < 3 { 0 } else { 1 });
        for dx in -reach[0]..=reach[0] {
            for dy in -reach[1]..=reach[1] {
                for dz in -reach[2]..=reach[2] {
                    let mut b = centre + IVec3::new(dx, dy, dz);
                    if use_pbc {
                        for k in 0..3 {
                            b[k] = ((b[k] % bin_counts[k]) + bin_counts[k]) % bin_counts[k];
                        }
                    }
                    if let Some(found) = site_bins.get(&bin_key(b)) {
                        for &index in found {
                            let d = g.nearest_image(p - regions[index].position);
                            if dot(d, d) < min_sep * min_sep {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        if use_pbc && bin_counts.iter().any(|&c| c < 3) {
            for c in regions {
                let d = g.nearest_image(p - c.position);
                if dot(d, d) < min_sep * min_sep {
                    return true;
                }
            }
        }
        false
    };

    for cand in &centers {
        let mut center = cand.center;
        if near_site(center, &regions, &site_bins) {
            continue;
        }
        let mut clearance = clearance_of(center);
        if clearance < min_clear {
            continue;
        }
        let mut shell = shell_at(center);
        if shell.len() < 4 {
            continue;
        }
        // Move to the centre of the surrounding atoms while that is roomier.
        for _ in 0..8 {
            let mut centroid = Vec3::ZERO;
            for sp in &shell {
                centroid += sp.position;
            }
            let centroid = g.wrap(centroid / shell.len() as f32);
            let refined = clearance_of(centroid);
            if refined > clearance + 1e-4 {
                let refined_shell = shell_at(centroid);
                if refined_shell.len() < 4 {
                    break;
                }
                clearance = refined;
                shell = refined_shell;
                center = centroid;
            } else {
                break;
            }
        }
        // A grid peak approximates a largest empty sphere: its four nearest atoms' circumcentre.
        if cand.a < 0 && shell.len() >= 4 {
            let mut nearest = shell.clone();
            sort_by(&mut nearest, |x, y| dot(x.position - center, x.position - center) < dot(y.position - center, y.position - center));
            if let Some((cc, rr)) = compute_circumsphere(nearest[0].position, nearest[1].position, nearest[2].position, nearest[3].position) {
                if length(cc - center) < clearance {
                    let wrapped = g.wrap(cc);
                    let refined = clearance_of(wrapped);
                    if refined > clearance + 1e-4 && (refined - rr).abs() < 1e-3 {
                        let refined_shell = shell_at(wrapped);
                        if refined_shell.len() >= 4 {
                            clearance = refined;
                            shell = refined_shell;
                            center = wrapped;
                        }
                    }
                }
            }
        }
        // A site is a critical point of the clearance: its touching atoms surround it.
        let touching: Vec<Vec3> = shell
            .iter()
            .filter_map(|sp| {
                let d = sp.position - center;
                let r = length(d);
                (r > 0.0 && r <= clearance * 1.002 + 2e-3).then(|| d / r)
            })
            .collect();
        if !origin_in_hull(&touching) {
            continue;
        }
        let Some((faces, face_vertices, volume)) = polyhedral_faces(center, &shell) else { continue };
        if volume <= 1e-6 {
            continue;
        }
        // Reject oversized false cavities.
        let mut r1 = f32::MAX;
        for sp in &shell {
            r1 = min_f(r1, length(sp.position - center));
        }
        let max_reasonable = (4.0f32 / 3.0) * 3.14159265 * max_f(0.0, r1 * 1.5).powf(3.0);
        if volume > max_reasonable {
            continue;
        }
        let coordination = shell.len() as i32;
        let kind = if coordination <= 4 {
            VoidKind::Tetrahedral
        } else if coordination <= 6 {
            VoidKind::Octahedral
        } else {
            VoidKind::Irregular
        };
        let region = VoidRegion {
            position: center,
            clearance,
            coordination,
            volume,
            volume_bin: 0,
            kind,
            neighbor_atom_indices: shell.iter().map(|sp| sp.atom_index).collect(),
            faces,
            face_vertices,
        };
        if !near_site(region.position, &regions, &site_bins) {
            site_bins.entry(bin_key(bin_of(region.position))).or_default().push(regions.len());
            regions.push(region);
        }
    }

    sort_by(&mut regions, |a, b| a.clearance > b.clearance);
    if limit_voids && regions.len() as i32 > max_voids {
        regions.truncate(max_voids as usize);
    }
    assign_volume_bins(&mut regions);
    result.success = !regions.is_empty();
    result.message = if result.success { format!("Detected {} void regions.", regions.len()) } else { "No void regions detected.".into() };
    result.regions = regions;
    result
}
