//! Cluster and void analyses (the C++ science/Clusters): bonded clusters
//! with centroids, radii of gyration and percolation; probe-accessible voids
//! on a grid with the atoms lining them.

use std::collections::{BTreeMap, HashMap, VecDeque};

use af_core::cxx::{nth_element_by, sort_by};
use af_core::json::Json;

use crate::analysis::{ji, R};
use crate::core::*;
use crate::data::Parameters;

const PI: f64 = std::f64::consts::PI;

pub(crate) fn union_root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

fn median_nearest(positions: &[Vec3], cell: &Mat3, pbc: Pbc) -> Result<f64, String> {
    let mut nearest = vec![f64::INFINITY; positions.len()];
    for n in neighbor_list(positions, cell, pbc, 6.0)? {
        let d = norm(n.vector);
        if d < nearest[n.i] {
            nearest[n.i] = d;
        }
    }
    let mut finite: Vec<f64> = nearest.into_iter().filter(|d| d.is_finite()).collect();
    if finite.is_empty() {
        return Err("No neighbours within 6 Angstrom; set atomic_radius_A".into());
    }
    let middle = finite.len() / 2;
    nth_element_by(&mut finite, middle, |a, b| a < b);
    Ok(finite[middle])
}

pub fn cluster_analysis(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    let pbc = p.pbc("pbc", [false; 3])?;
    if pbc.iter().any(|&b| b) && !p.has("cell") {
        return Err("Periodic analysis requires a cell".into());
    }
    let cell = if p.has("cell") { cell_matrix(p.array("cell")?)? } else { identity() };
    let cutoff = positive(p.number_or("cutoff_A", 3.0)?, "cutoff_A")?;
    let selected: Vec<usize> = if p.has("mask") {
        let mask = finite_array(p.array("mask")?, "mask", 1)?.values;
        if mask.len() != positions.len() {
            return Err("mask needs one value per atom".into());
        }
        (0..mask.len()).filter(|&i| mask[i] != 0.0).collect()
    } else {
        (0..positions.len()).collect()
    };
    let minimum = integer(p.number_or("min_size", 1.0)?, "min_size", 1)?;
    let chosen: Vec<Vec3> = selected.iter().map(|&i| positions[i]).collect();
    let links = neighbor_list(&chosen, &cell, pbc, cutoff)?;
    let mut parent: Vec<usize> = (0..chosen.len()).collect();
    let mut adjacency: Vec<Vec<(usize, Vec3)>> = vec![Vec::new(); chosen.len()];
    for link in &links {
        // C++17: the right-hand side of = first.
        let rj = union_root(&mut parent, link.j);
        let ri = union_root(&mut parent, link.i);
        parent[ri] = rj;
        adjacency[link.i].push((link.j, link.vector));
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for k in 0..chosen.len() {
        let root = union_root(&mut parent, k);
        groups.entry(root).or_default().push(k);
    }
    let mut clusters: Vec<Vec<usize>> = groups.into_values().filter(|m| m.len() as i64 >= minimum).collect();
    sort_by(&mut clusters, |a, b| a.len() > b.len());
    let mut cluster_id = vec![-1i32; positions.len()];
    let mut list = Vec::new();
    let mut size_counts: BTreeMap<usize, i32> = BTreeMap::new();
    for (c, members) in clusters.iter().enumerate() {
        let mut unwrapped: BTreeMap<usize, Vec3> = BTreeMap::new();
        unwrapped.insert(members[0], chosen[members[0]]);
        let mut queue = VecDeque::from([members[0]]);
        let mut percolates = false;
        while let Some(k) = queue.pop_front() {
            for &(j, vector) in &adjacency[k] {
                let reached = add(unwrapped[&k], vector);
                match unwrapped.get(&j) {
                    None => {
                        unwrapped.insert(j, reached);
                        queue.push_back(j);
                    }
                    Some(&existing) => {
                        if norm(sub(reached, existing)) > 1e-6 {
                            percolates = true;
                        }
                    }
                }
            }
        }
        let mut centroid = [0.0; 3];
        for &k in members {
            centroid = add(centroid, unwrapped[&k]);
            cluster_id[selected[k]] = c as i32;
        }
        centroid = scale(centroid, 1.0 / members.len() as f64);
        let mut gyration = 0.0;
        for &k in members {
            let d = sub(unwrapped[&k], centroid);
            gyration += dot(d, d);
        }
        gyration = (gyration / members.len() as f64).sqrt();
        let mut entry = Json::object();
        entry.set("size", members.len() as f64);
        entry.set("centroid_A", vec3_json(centroid));
        entry.set("radius_of_gyration_A", gyration);
        entry.set("percolating", percolates);
        list.push(entry);
        *size_counts.entry(members.len()).or_insert(0) += 1;
    }
    let histogram = size_counts.iter().map(|(&size, &count)| Json::Array(vec![Json::Number(size as f64), Json::Number(count as f64)])).collect();
    let mut result = Json::object();
    result.set("cluster_count", clusters.len() as f64);
    result.set("selected_atoms", selected.len() as f64);
    result.set("clusters", Json::Array(list));
    result.set("size_histogram", Json::Array(histogram));
    result.set("cluster_id", ji(&cluster_id));
    Ok(result)
}

struct Image {
    x: Vec3,
    atom: usize,
}

#[derive(Default)]
struct Void {
    points: usize,
    sum: Vec3,
    exterior: bool,
}

pub fn void_analysis(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    let pbc = p.pbc("pbc", [false; 3])?;
    if pbc.iter().any(|&b| b) && !p.has("cell") {
        return Err("Periodic analysis requires a cell".into());
    }
    let cell = if p.has("cell") { cell_matrix(p.array("cell")?)? } else { identity() };
    let mut radius = p.number_or("atomic_radius_A", 0.0)?;
    if radius <= 0.0 {
        radius = 0.5 * median_nearest(&positions, &cell, pbc)?;
    }
    let probe = p.number_or("probe_radius_A", 1.0)?;
    if !(probe >= 0.0) {
        return Err("probe_radius_A must be nonnegative".into());
    }
    let spacing = positive(p.number_or("grid_spacing_A", 0.3)?, "grid_spacing_A")?;
    let reach = radius + probe;
    let mut domain = cell;
    let mut origin = [0.0; 3];
    if !p.has("cell") {
        let (mut low, mut high) = (positions[0], positions[0]);
        for x in &positions {
            for k in 0..3 {
                if x[k] < low[k] {
                    low[k] = x[k];
                }
                if high[k] < x[k] {
                    high[k] = x[k];
                }
            }
        }
        let margin = reach + 2.0 * spacing;
        domain = [
            [high[0] - low[0] + 2.0 * margin, 0.0, 0.0],
            [0.0, high[1] - low[1] + 2.0 * margin, 0.0],
            [0.0, 0.0, high[2] - low[2] + 2.0 * margin],
        ];
        origin = [low[0] - margin, low[1] - margin, low[2] - margin];
    }
    let inverse_domain = inverse(&domain)?;
    let size: [i32; 3] = std::array::from_fn(|k| 4.max((norm(domain[k]) / spacing).ceil() as i32));
    let total = size[0] as i64 * size[1] as i64 * size[2] as i64;
    if total > 60_000_000 {
        return Err("Void grid too large; increase grid_spacing_A".into());
    }
    let mut images = Vec::new();
    let mut range = [0i32; 3];
    for k in 0..3 {
        if pbc[k] {
            range[k] = (reach * norm([inverse_domain[0][k], inverse_domain[1][k], inverse_domain[2][k]])).ceil() as i32 + 1;
        }
    }
    for (a, x) in positions.iter().enumerate() {
        let mut f = row_times(sub(*x, origin), &inverse_domain);
        for k in 0..3 {
            if pbc[k] {
                f[k] -= f[k].floor();
            }
        }
        for i in -range[0]..=range[0] {
            for j in -range[1]..=range[1] {
                for k in -range[2]..=range[2] {
                    images.push(Image { x: add(origin, row_times([f[0] + i as f64, f[1] + j as f64, f[2] + k as f64], &domain)), atom: a });
                }
            }
        }
    }
    let key = |x: Vec3| [(x[0] / reach).floor() as i64, (x[1] / reach).floor() as i64, (x[2] / reach).floor() as i64];
    let hash = |k: [i64; 3]| (k[0].wrapping_mul(73856093)) ^ (k[1].wrapping_mul(19349663)) ^ (k[2].wrapping_mul(83492791));
    let mut bins: HashMap<i64, Vec<usize>> = HashMap::new();
    for (i, image) in images.iter().enumerate() {
        bins.entry(hash(key(image.x))).or_default().push(i);
    }
    let point_at = |i: i32, j: i32, k: i32| {
        add(origin, row_times([(i as f64 + 0.5) / size[0] as f64, (j as f64 + 0.5) / size[1] as f64, (k as f64 + 0.5) / size[2] as f64], &domain))
    };
    let index = |i: i32, j: i32, k: i32| (i as usize * size[1] as usize + j as usize) * size[2] as usize + k as usize;
    let mut accessible = vec![false; total as usize];
    for i in 0..size[0] {
        for j in 0..size[1] {
            for k in 0..size[2] {
                let x = point_at(i, j, k);
                let c = key(x);
                let mut blocked = false;
                'search: for dx in -1..=1i64 {
                    for dy in -1..=1i64 {
                        for dz in -1..=1i64 {
                            let Some(found) = bins.get(&hash([c[0] + dx, c[1] + dy, c[2] + dz])) else { continue };
                            for &index2 in found {
                                let d = sub(images[index2].x, x);
                                if dot(d, d) < reach * reach {
                                    blocked = true;
                                    break 'search;
                                }
                            }
                        }
                    }
                }
                accessible[index(i, j, k)] = !blocked;
            }
        }
    }
    let mut label = vec![-1i32; total as usize];
    let voxel = cell_volume(&domain) / total as f64;
    let mut voids: Vec<Void> = Vec::new();
    let steps = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
    for i in 0..size[0] {
        for j in 0..size[1] {
            for k in 0..size[2] {
                if !accessible[index(i, j, k)] || label[index(i, j, k)] >= 0 {
                    continue;
                }
                let mut v = Void::default();
                let id = voids.len() as i32;
                let mut queue = VecDeque::from([[i, j, k, 0, 0, 0]]);
                label[index(i, j, k)] = id;
                while let Some(q) = queue.pop_front() {
                    v.points += 1;
                    v.sum = add(
                        v.sum,
                        [
                            (q[0] as f64 + 0.5) / size[0] as f64 + q[3] as f64,
                            (q[1] as f64 + 0.5) / size[1] as f64 + q[4] as f64,
                            (q[2] as f64 + 0.5) / size[2] as f64 + q[5] as f64,
                        ],
                    );
                    for s in &steps {
                        let mut g = [q[0] + s[0], q[1] + s[1], q[2] + s[2]];
                        let mut offset = [q[3], q[4], q[5]];
                        let mut outside = false;
                        for axis in 0..3 {
                            if g[axis] >= 0 && g[axis] < size[axis] {
                                continue;
                            }
                            if !pbc[axis] {
                                outside = true;
                                continue;
                            }
                            let wrap = if g[axis] < 0 { -1 } else { 1 };
                            g[axis] -= wrap * size[axis];
                            offset[axis] += wrap;
                        }
                        if outside {
                            v.exterior = true;
                            continue;
                        }
                        let n = index(g[0], g[1], g[2]);
                        if !accessible[n] || label[n] >= 0 {
                            continue;
                        }
                        label[n] = id;
                        queue.push_back([g[0], g[1], g[2], offset[0], offset[1], offset[2]]);
                    }
                }
                voids.push(v);
            }
        }
    }
    let mut order: Vec<usize> = (0..voids.len()).filter(|&v| !voids[v].exterior).collect();
    sort_by(&mut order, |&a, &b| voids[a].points > voids[b].points);
    let rank: BTreeMap<usize, i32> = order.iter().enumerate().map(|(r, &v)| (v, r as i32)).collect();
    let mut lining = vec![-1i32; positions.len()];
    let near_distance = reach + 1.5 * spacing;
    let span = (near_distance / reach).ceil() as i64;
    for i in 0..size[0] {
        for j in 0..size[1] {
            for k in 0..size[2] {
                let id = label[index(i, j, k)];
                if id < 0 {
                    continue;
                }
                let Some(&r) = rank.get(&(id as usize)) else { continue };
                let x = point_at(i, j, k);
                let c = key(x);
                for dx in -span..=span {
                    for dy in -span..=span {
                        for dz in -span..=span {
                            let Some(found) = bins.get(&hash([c[0] + dx, c[1] + dy, c[2] + dz])) else { continue };
                            for &image in found {
                                if norm(sub(images[image].x, x)) < near_distance {
                                    let owner = &mut lining[images[image].atom];
                                    if *owner < 0 {
                                        *owner = r;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    let mut list = Vec::new();
    let mut total_volume = 0.0;
    for &id in &order {
        let v = &voids[id];
        let volume = v.points as f64 * voxel;
        total_volume += volume;
        let centroid = add(origin, row_times(scale(v.sum, 1.0 / v.points as f64), &domain));
        let equivalent = (3.0 * volume / (4.0 * PI)).cbrt();
        let mut entry = Json::object();
        entry.set("accessible_volume_A3", volume);
        entry.set("centroid_A", vec3_json(centroid));
        entry.set("equivalent_radius_A", equivalent);
        entry.set("pore_radius_estimate_A", equivalent + reach);
        entry.set("grid_points", v.points as f64);
        list.push(entry);
    }
    let mut result = Json::object();
    result.set("void_count", order.len() as f64);
    result.set("voids", Json::Array(list));
    result.set("accessible_volume_A3", total_volume);
    result.set("void_fraction", if p.has("cell") { Json::Number(total_volume / cell_volume(&cell)) } else { Json::Null });
    result.set("atomic_radius_A", radius);
    result.set("probe_radius_A", probe);
    result.set("grid", Json::Array(size.iter().map(|&s| Json::Number(s as f64)).collect()));
    result.set("lining_void_id", ji(&lining));
    Ok(result)
}
