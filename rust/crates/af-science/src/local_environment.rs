//! Local-environment analyses (the C++ science/LocalEnvironment): Chill+
//! ice and hydrate types, the local entropy fingerprint, shortest-path
//! rings, displacement vectors and bond length/angle distributions.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::rc::Rc;

use af_core::json::Json;
use num_complex::Complex64;

use crate::analysis::{ji, jv, neighbors_by_center, normalized_legendre, R};
use crate::core::*;
use crate::data::Parameters;

const PI: f64 = std::f64::consts::PI;

/// Positions, cell and periodicity of a request (LocalGeometry's AtomSet).
pub(crate) struct AtomSet {
    pub positions: Vec<Vec3>,
    pub cell: Mat3,
    pub pbc: Pbc,
    #[allow(dead_code)]
    pub has_cell: bool,
}

pub(crate) fn read_atom_set(p: &Parameters, positions: &str) -> Result<AtomSet, String> {
    let points_ = points(p.array(positions)?, positions)?;
    let pbc = p.pbc("pbc", [false; 3])?;
    let has_cell = p.has("cell");
    if pbc.iter().any(|&b| b) && !has_cell {
        return Err("Periodic analysis requires a cell".into());
    }
    let cell = if has_cell { cell_matrix(p.array("cell")?)? } else { identity() };
    Ok(AtomSet { positions: points_, cell, pbc, has_cell })
}

/// The classification `key` with its names, counts and fractions.
pub(crate) fn add_classification(result: &mut Json, key: &str, types: &[i32], names: &[&str], prefix: &str) {
    let mut counts = vec![0i64; names.len()];
    for &t in types {
        counts[t as usize] += 1;
    }
    let (mut count_json, mut fractions) = (Json::object(), Json::object());
    for (t, name) in names.iter().enumerate() {
        count_json.set(name, counts[t] as f64);
        fractions.set(name, if types.is_empty() { 0.0 } else { counts[t] as f64 / types.len() as f64 });
    }
    result.set(key, ji(types));
    result.set(&if prefix.is_empty() { "type_names".to_string() } else { format!("{prefix}names") }, Json::Array(names.iter().map(|n| Json::String(n.to_string())).collect()));
    result.set(&format!("{prefix}counts"), count_json);
    result.set(&format!("{prefix}fractions"), fractions);
}

pub(crate) fn first_of_pair(n: &Neighbor) -> bool {
    n.i < n.j || (n.i == n.j && n.shift > [0, 0, 0])
}

pub fn chill_plus(p: &Parameters) -> R {
    let a = read_atom_set(p, "positions")?;
    let cutoff = positive(p.number_or("cutoff_A", 3.5)?, "cutoff_A")?;
    let neighbors = neighbor_list(&a.positions, &a.cell, a.pbc, cutoff)?;
    let groups = neighbors_by_center(&neighbors, a.positions.len());
    let n = a.positions.len();
    let zero = Complex64::new(0.0, 0.0);
    let mut q = vec![[zero; 4]; n];
    for i in 0..n {
        for &g in &groups[i] {
            let v = neighbors[g].vector;
            let r = norm(v);
            if r < 1e-10 {
                return Err("Coincident atoms make Chill+ undefined".into());
            }
            let cosine = (v[2] / r).clamp(-1.0, 1.0);
            let azimuth = v[1].atan2(v[0]);
            for m in 0..=3 {
                let l = normalized_legendre(3, m, cosine);
                let angle = m as f64 * azimuth;
                q[i][m as usize] += Complex64::new(l * angle.cos(), l * angle.sin());
            }
        }
    }
    let correlation = |i: usize, j: usize| {
        let mut dot_product = (q[i][0] * q[j][0].conj()).re;
        let (mut ni, mut nj) = (q[i][0].norm_sqr(), q[j][0].norm_sqr());
        for m in 1..=3 {
            dot_product += 2.0 * (q[i][m] * q[j][m].conj()).re;
            ni += 2.0 * q[i][m].norm_sqr();
            nj += 2.0 * q[j][m].norm_sqr();
        }
        if ni > 0.0 && nj > 0.0 {
            dot_product / (ni * nj).sqrt()
        } else {
            f64::NAN
        }
    };
    let names = ["other", "hexagonal ice", "cubic ice", "interfacial ice", "hydrate", "interfacial hydrate"];
    let (mut kind, mut staggered_count, mut eclipsed_count) = (vec![0i32; n], vec![0i32; n], vec![0i32; n]);
    for i in 0..n {
        if groups[i].len() != 4 {
            continue;
        }
        let (mut staggered, mut eclipsed) = (0, 0);
        for &g in &groups[i] {
            let c = correlation(i, neighbors[g].j);
            if c <= -0.8 {
                staggered += 1;
            } else if (-0.35..=0.25).contains(&c) {
                eclipsed += 1;
            }
        }
        staggered_count[i] = staggered;
        eclipsed_count[i] = eclipsed;
        kind[i] = if eclipsed == 4 {
            4
        } else if eclipsed == 3 {
            5
        } else if staggered == 4 {
            2
        } else if staggered == 3 && eclipsed == 1 {
            1
        } else if staggered == 3 || staggered == 2 {
            3
        } else {
            0
        };
    }
    let mut result = Json::object();
    add_classification(&mut result, "chill_type", &kind, &names, "");
    result.set("staggered_bonds", ji(&staggered_count));
    result.set("eclipsed_bonds", ji(&eclipsed_count));
    result.set("cutoff_A", cutoff);
    Ok(result)
}

pub fn local_entropy(p: &Parameters) -> R {
    let a = read_atom_set(p, "positions")?;
    let cutoff = positive(p.number_or("cutoff_A", 5.0)?, "cutoff_A")?;
    let sigma = positive(p.number_or("sigma_A", 0.2)?, "sigma_A")?;
    let local_density = p.boolean("local_density", false)?;
    let average = p.boolean("average", true)?;
    let average_cutoff = positive(p.number_or("average_cutoff_A", cutoff)?, "average_cutoff_A")?;
    let n = a.positions.len();
    if n == 0 {
        return Err("No atoms".into());
    }
    let mut global_density = f64::NAN;
    if p.has("cell") && a.pbc[0] && a.pbc[1] && a.pbc[2] {
        global_density = n as f64 / cell_volume(&a.cell);
    }
    if !local_density && !global_density.is_finite() {
        return Err("The global density needs a fully periodic cell; set local_density for other systems".into());
    }
    let reach = if cutoff + 3.0 * sigma < average_cutoff { average_cutoff } else { cutoff + 3.0 * sigma };
    let neighbors = neighbor_list(&a.positions, &a.cell, a.pbc, reach)?;
    let groups = neighbors_by_center(&neighbors, n);
    let dr = sigma / 10.0;
    let bins = (cutoff / dr).ceil() as usize;
    let radius: Vec<f64> = (0..bins).map(|b| (b as f64 + 0.5) * cutoff / bins as f64).collect();
    let width = cutoff / bins as f64;
    let gaussian_norm = 1.0 / (2.0 * PI * sigma * sigma).sqrt();
    let mut entropy = vec![0.0; n];
    let mut g = vec![0.0; bins];
    for i in 0..n {
        let mut density = global_density;
        if local_density {
            let inside = groups[i].iter().filter(|&&k| norm(neighbors[k].vector) <= cutoff).count();
            density = inside as f64 / (4.0 / 3.0 * PI * cutoff * cutoff * cutoff);
        }
        if !(density > 0.0) {
            entropy[i] = 0.0;
            continue;
        }
        g.iter_mut().for_each(|v| *v = 0.0);
        for &k in &groups[i] {
            let rij = norm(neighbors[k].vector);
            let lo = (rij - 5.0 * sigma) / width;
            let low = (if 0.0 < lo { lo } else { 0.0 }) as usize;
            let hi = (rij + 5.0 * sigma) / width;
            let high = bins.min((if 0.0 < hi { hi } else { 0.0 }) as usize + 1);
            for b in low..high {
                let d = radius[b] - rij;
                g[b] += gaussian_norm * (-d * d / (2.0 * sigma * sigma)).exp();
            }
        }
        let mut integral = 0.0;
        for b in 0..bins {
            let r = radius[b];
            let gr = g[b] / (4.0 * PI * density * r * r);
            let integrand = (if gr > 0.0 { gr * gr.ln() } else { 0.0 }) - gr + 1.0;
            integral += integrand * r * r * width;
        }
        entropy[i] = -2.0 * PI * density * integral;
    }
    let mut averaged = entropy.clone();
    if average {
        for i in 0..n {
            let (mut sum, mut count) = (entropy[i], 1.0);
            for &k in &groups[i] {
                if norm(neighbors[k].vector) <= average_cutoff {
                    sum += entropy[neighbors[k].j];
                    count += 1.0;
                }
            }
            averaged[i] = sum / count;
        }
    }
    let mut mean = 0.0;
    for s in &averaged {
        mean += s;
    }
    let mut result = Json::object();
    result.set("entropy_kB", jv(&averaged));
    if average {
        result.set("unaveraged_entropy_kB", jv(&entropy));
    }
    result.set("mean_entropy_kB", mean / n as f64);
    result.set("cutoff_A", cutoff);
    result.set("sigma_A", sigma);
    Ok(result)
}

// ---------------------------------------------------------------- rings
/// An atom and its periodic image.
type Node = [i32; 4];

#[derive(Clone, Copy)]
struct Bond {
    j: i32,
    shift: [i32; 3],
}

struct Search {
    distance: HashMap<Node, i32>,
    parents: HashMap<Node, Vec<Node>>,
    levels: Vec<Vec<Node>>,
}

struct Network {
    bonds: Vec<Vec<Bond>>,
    distances: HashMap<i32, HashMap<Node, i32>>,
}

impl Network {
    fn new(atoms: &AtomSet, cutoff: f64) -> Result<Self, String> {
        let mut bonds = vec![Vec::new(); atoms.positions.len()];
        for link in neighbor_list(&atoms.positions, &atoms.cell, atoms.pbc, cutoff)? {
            bonds[link.i].push(Bond { j: link.j as i32, shift: link.shift });
        }
        Ok(Self { bonds, distances: HashMap::new() })
    }

    fn step(from: &Node, b: &Bond) -> Node {
        [b.j, from[1] + b.shift[0], from[2] + b.shift[1], from[3] + b.shift[2]]
    }

    fn search(&self, atom: i32, depth: i32, with_parents: bool) -> Search {
        let source = [atom, 0, 0, 0];
        let mut s = Search { distance: HashMap::from([(source, 0)]), parents: HashMap::new(), levels: vec![vec![source]] };
        let mut d = 0;
        while d < depth && !s.levels.last().unwrap().is_empty() {
            let mut next = Vec::new();
            for node in s.levels.last().unwrap().clone() {
                for b in &self.bonds[node[0] as usize] {
                    let to = Self::step(&node, b);
                    let entry = *s.distance.entry(to).or_insert_with(|| {
                        next.push(to);
                        d + 1
                    });
                    if with_parents && entry == d + 1 {
                        s.parents.entry(to).or_default().push(node);
                    }
                }
            }
            s.levels.push(next);
            d += 1;
        }
        s
    }

    fn distance(&mut self, from: &Node, to: &Node, depth: i32) -> i32 {
        if !self.distances.contains_key(&from[0]) {
            let found = self.search(from[0], depth, false).distance;
            self.distances.insert(from[0], found);
        }
        let map = &self.distances[&from[0]];
        *map.get(&[to[0], to[1] - from[1], to[2] - from[2], to[3] - from[3]]).unwrap_or(&(depth + 1))
    }
}

fn canonical_ring(ring: &[Node]) -> Vec<Node> {
    let m = ring.len();
    let mut best: Vec<Node> = Vec::new();
    for start in 0..m {
        for direction in [1i32, -1] {
            let origin = ring[start];
            let sequence: Vec<Node> = (0..m)
                .map(|k| {
                    let node = ring[if direction > 0 { (start + k) % m } else { (start + m - k % m) % m }];
                    [node[0], node[1] - origin[1], node[2] - origin[2], node[3] - origin[3]]
                })
                .collect();
            if best.is_empty() || sequence < best {
                best = sequence;
            }
        }
    }
    best
}

struct PathCache<'a> {
    search: &'a Search,
    paths: HashMap<Node, Rc<Vec<Vec<Node>>>>,
    limit: usize,
    truncated: bool,
}

impl PathCache<'_> {
    /// All shortest paths from the source to a node (bounded in number).
    fn paths_to(&mut self, node: &Node) -> Rc<Vec<Vec<Node>>> {
        if let Some(found) = self.paths.get(node) {
            return found.clone();
        }
        let mut result: Vec<Vec<Node>> = Vec::new();
        if self.search.distance[node] == 0 {
            result.push(vec![*node]);
        } else {
            for parent in self.search.parents[node].clone() {
                for path in self.paths_to(&parent).iter() {
                    if result.len() >= self.limit {
                        self.truncated = true;
                        break;
                    }
                    let mut extended = path.clone();
                    extended.push(*node);
                    result.push(extended);
                }
            }
        }
        let result = Rc::new(result);
        self.paths.insert(*node, result.clone());
        result
    }
}

fn disjoint(x: &[Node], x_from: usize, x_to: usize, y: &[Node], y_from: usize, y_to: usize) -> bool {
    for i in x_from..x_to {
        for j in y_from..y_to {
            if x[i] == y[j] {
                return false;
            }
        }
    }
    true
}

pub fn find_rings(p: &Parameters) -> R {
    let a = read_atom_set(p, "positions")?;
    let cutoff = positive(p.number_or("cutoff_A", 2.0)?, "cutoff_A")?;
    let max_size = integer(p.number_or("max_size", 10.0)?, "max_size", 3)?;
    if max_size > 24 {
        return Err("max_size must be at most 24".into());
    }
    let n = a.positions.len();
    let half = (max_size / 2) as i32;
    let mut network = Network::new(&a, cutoff)?;
    let max_paths = integer(p.number_or("max_paths", 512.0)?, "max_paths", 0)?;
    let path_limit = if max_paths == 0 { usize::MAX } else { max_paths as usize };
    let mut truncated = false;
    let (mut rings, mut rejected): (BTreeSet<Vec<Node>>, BTreeSet<Vec<Node>>) = (BTreeSet::new(), BTreeSet::new());
    let shortest_path_ring = |network: &mut Network, ring: &[Node]| {
        let m = ring.len();
        for i in 0..m {
            for j in i + 2..m {
                let around = (j - i).min(m - (j - i));
                if around < 2 {
                    continue;
                }
                if network.distance(&ring[i], &ring[j], half) < around as i32 {
                    return false;
                }
            }
        }
        true
    };
    for r in 0..n {
        let source = r as i32;
        let s = network.search(source, half, true);
        let mut cache = PathCache { search: &s, paths: HashMap::new(), limit: path_limit, truncated: false };
        let mut candidates: Vec<Vec<Node>> = Vec::new();
        let mut d = 1;
        while d <= half && (d as usize) < s.levels.len() {
            for m in &s.levels[d as usize] {
                if (2 * d) as i64 <= max_size {
                    let to_m = cache.paths_to(m);
                    for i in 0..to_m.len() {
                        for j in i + 1..to_m.len() {
                            let (x, y) = (&to_m[i], &to_m[j]);
                            if !disjoint(x, 1, x.len() - 1, y, 1, y.len() - 1) {
                                continue;
                            }
                            let mut ring = x.clone();
                            for k in (1..y.len() - 1).rev() {
                                ring.push(y[k]);
                            }
                            candidates.push(ring);
                        }
                    }
                }
                if (2 * d + 1) as i64 <= max_size {
                    for b in network.bonds[m[0] as usize].clone() {
                        let other = Network::step(m, &b);
                        match s.distance.get(&other) {
                            Some(&od) if od == d && *m < other => {}
                            _ => continue,
                        }
                        let (xs, ys) = (cache.paths_to(m), cache.paths_to(&other));
                        for x in xs.iter() {
                            for y in ys.iter() {
                                if !disjoint(x, 1, x.len(), y, 1, y.len()) {
                                    continue;
                                }
                                let mut ring = x.clone();
                                for k in (1..y.len()).rev() {
                                    ring.push(y[k]);
                                }
                                candidates.push(ring);
                            }
                        }
                    }
                }
            }
            // Candidates of this level, considered in order.
            for ring in candidates.drain(..) {
                if ring.iter().any(|node| node[0] < source) {
                    continue;
                }
                let key = canonical_ring(&ring);
                if rings.contains(&key) || rejected.contains(&key) {
                    continue;
                }
                if shortest_path_ring(&mut network, &ring) {
                    rings.insert(key);
                } else {
                    rejected.insert(key);
                }
            }
            d += 1;
        }
        truncated = truncated || cache.truncated;
    }
    let mut sizes: BTreeMap<usize, i64> = BTreeMap::new();
    let mut per_atom = vec![0i32; n];
    let mut list = Vec::new();
    for ring in &rings {
        *sizes.entry(ring.len()).or_insert(0) += 1;
        for node in ring {
            per_atom[node[0] as usize] += 1;
        }
        if list.len() < 10000 {
            list.push(Json::Array(ring.iter().map(|node| Json::Number(node[0] as f64)).collect()));
        }
    }
    let mut result = Json::object();
    result.set("ring_sizes", jv(&sizes.keys().map(|&k| k as f64).collect::<Vec<_>>()));
    result.set("ring_counts", jv(&sizes.values().map(|&v| v as f64).collect::<Vec<_>>()));
    result.set("total_rings", rings.len() as f64);
    result.set("rings_per_atom", ji(&per_atom));
    let listed = list.len();
    result.set("rings", Json::Array(list));
    result.set("rings_listed", listed as f64);
    result.set("paths_truncated", truncated);
    result.set("cutoff_A", cutoff);
    result.set("max_size", max_size as f64);
    result.set("max_paths", max_paths as f64);
    Ok(result)
}

pub fn displacement_vectors(p: &Parameters) -> R {
    let reference = points(p.array("reference")?, "reference")?;
    let current = points(p.array("current")?, "current")?;
    if reference.len() != current.len() {
        return Err("reference and current need the same number of atoms".into());
    }
    let pbc = p.pbc("pbc", [false; 3])?;
    let periodic = pbc.iter().any(|&b| b);
    if periodic && !p.has("reference_cell") {
        return Err("Periodic displacements need reference_cell".into());
    }
    let reference_cell = if p.has("reference_cell") { cell_matrix(p.array("reference_cell")?)? } else { identity() };
    let current_cell = if p.has("current_cell") { cell_matrix(p.array("current_cell")?)? } else { reference_cell };
    let mapping = if p.has("affine_mapping") { p.json("affine_mapping")?.string()?.to_string() } else { "off".to_string() };
    if mapping != "off" && mapping != "to_reference" && mapping != "to_current" {
        return Err("affine_mapping must be off, to_reference or to_current".into());
    }
    if mapping != "off" && !(p.has("reference_cell") && p.has("current_cell")) {
        return Err("Affine mapping needs reference_cell and current_cell".into());
    }
    let image_cell = if mapping == "to_current" { current_cell } else { reference_cell };
    let mic = MinimumImage::new(&image_cell, pbc)?;
    let n = reference.len();
    let mut vectors = Vec::new();
    let mut magnitude = vec![0.0; n];
    let (mut sum, mut sum_squares, mut largest) = (0.0, 0.0, 0.0f64);
    for i in 0..n {
        let (mut from, mut to) = (reference[i], current[i]);
        if mapping == "to_reference" {
            to = row_times(fractional(to, &current_cell)?, &reference_cell);
        } else if mapping == "to_current" {
            from = row_times(fractional(from, &reference_cell)?, &current_cell);
        }
        let mut d = sub(to, from);
        if periodic {
            d = mic.apply(d);
        }
        vectors.push(vec3_json(d));
        magnitude[i] = norm(d);
        sum += magnitude[i];
        sum_squares += magnitude[i] * magnitude[i];
        if largest < magnitude[i] {
            largest = magnitude[i];
        }
    }
    let mut result = Json::object();
    result.set("displacement_A", Json::Array(vectors));
    result.set("magnitude_A", jv(&magnitude));
    result.set("mean_magnitude_A", if n > 0 { sum / n as f64 } else { f64::NAN });
    result.set("rms_displacement_A", if n > 0 { (sum_squares / n as f64).sqrt() } else { f64::NAN });
    result.set("max_magnitude_A", largest);
    Ok(result)
}

pub fn bond_distributions(p: &Parameters) -> R {
    let a = read_atom_set(p, "positions")?;
    let n = a.positions.len();
    let cutoff = positive(p.number("cutoff_A")?, "cutoff_A")?;
    let bins = integer(p.number_or("bins", 100.0)?, "bins", 4)?;
    let mut species = vec![0.0; n];
    if p.has("types") {
        species = finite_array(p.array("types")?, "types", 1)?.values;
        if species.len() != n {
            return Err("types needs one value per atom".into());
        }
    }
    let neighbors = neighbor_list(&a.positions, &a.cell, a.pbc, cutoff)?;
    let groups = neighbors_by_center(&neighbors, n);
    let name = |t: f64| (t.round() as i64).to_string();
    let fmin = |x: f64, y: f64| if y < x { y } else { x };
    let fmax = |x: f64, y: f64| if x < y { y } else { x };
    let length_width = cutoff / bins as f64;
    let angle_width = 180.0 / bins as f64;
    let nb = bins as usize;
    let (mut lengths, mut angles): (BTreeMap<String, Vec<f64>>, BTreeMap<String, Vec<f64>>) = (BTreeMap::new(), BTreeMap::new());
    let (mut all_lengths, mut all_angles) = (vec![0.0; nb], vec![0.0; nb]);
    let mut bond_lengths = Vec::new();
    let (mut bond_count, mut angle_count) = (0i64, 0i64);
    for link in &neighbors {
        if !first_of_pair(link) {
            continue;
        }
        let r = norm(link.vector);
        if r <= 0.0 {
            continue;
        }
        let b = (nb - 1).min((r / length_width) as usize);
        let (ti, tj) = (species[link.i], species[link.j]);
        let key = format!("{}-{}", name(fmin(ti, tj)), name(fmax(ti, tj)));
        let histogram = lengths.entry(key).or_insert_with(|| vec![0.0; nb]);
        histogram[b] += 1.0;
        all_lengths[b] += 1.0;
        bond_lengths.push(r);
        bond_count += 1;
    }
    for i in 0..n {
        let group = &groups[i];
        for x in 0..group.len() {
            for y in x + 1..group.len() {
                let (u, v) = (&neighbors[group[x]], &neighbors[group[y]]);
                let c = (dot(u.vector, v.vector) / (norm(u.vector) * norm(v.vector))).clamp(-1.0, 1.0);
                let angle = c.acos() * 180.0 / PI;
                let b = (nb - 1).min((angle / angle_width) as usize);
                let (tu, tv) = (species[u.j], species[v.j]);
                let key = format!("{}-{}-{}", name(fmin(tu, tv)), name(species[i]), name(fmax(tu, tv)));
                let histogram = angles.entry(key).or_insert_with(|| vec![0.0; nb]);
                histogram[b] += 1.0;
                all_angles[b] += 1.0;
                angle_count += 1;
            }
        }
    }
    let length_centers: Vec<f64> = (0..bins).map(|b| (b as f64 + 0.5) * length_width).collect();
    let angle_centers: Vec<f64> = (0..bins).map(|b| (b as f64 + 0.5) * angle_width).collect();
    let mut result = Json::object();
    result.set("length_A", jv(&length_centers));
    result.set("angle_deg", jv(&angle_centers));
    result.set("bond_length_counts", jv(&all_lengths));
    result.set("bond_angle_counts", jv(&all_angles));
    let (mut by_pair, mut by_triplet) = (Json::object(), Json::object());
    for (key, histogram) in &lengths {
        by_pair.set(key, jv(histogram));
    }
    for (key, histogram) in &angles {
        by_triplet.set(key, jv(histogram));
    }
    result.set("bond_length_counts_by_pair", by_pair);
    result.set("bond_angle_counts_by_triplet", by_triplet);
    result.set("bonds", bond_count as f64);
    result.set("angles", angle_count as f64);
    if !bond_lengths.is_empty() {
        let mut sum = 0.0;
        for r in &bond_lengths {
            sum += r;
        }
        result.set("mean_bond_length_A", sum / bond_lengths.len() as f64);
    }
    Ok(result)
}
