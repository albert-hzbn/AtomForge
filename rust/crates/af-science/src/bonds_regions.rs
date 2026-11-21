//! Bond orders from covalent radii and LAMMPS region membership with region
//! surfaces (the C++ science/BondsAndRegions and MeshText).

use std::collections::BTreeMap;

use af_core::format::format_g;
use af_core::json::Json;
use af_electronic::{isosurface, Grid};
use glam::{DMat3, DVec3};

use crate::analysis::{ji, jv, R};
use crate::core::*;
use crate::data::Parameters;
use crate::local_environment::first_of_pair;

const HUGE: f64 = 1e30;

/// OBJ text of named triangle sets (ostream precision 8).
pub struct ObjWriter {
    text: String,
    vertices: usize,
}

impl ObjWriter {
    pub fn new(comment: &str) -> Self {
        Self { text: format!("# {comment}\n"), vertices: 0 }
    }

    /// An object of triangles (three vertices each; shared vertices are not merged).
    pub fn add(&mut self, name: &str, triangles: &[DVec3]) {
        if !name.is_empty() {
            self.text += &format!("o {name}\n");
        }
        for v in triangles {
            self.text += &format!("v {} {} {}\n", format_g(v.x, 8), format_g(v.y, 8), format_g(v.z, 8));
        }
        let mut t = 0;
        while t + 2 < triangles.len() {
            self.text += &format!("f {} {} {}\n", self.vertices + t + 1, self.vertices + t + 2, self.vertices + t + 3);
            t += 3;
        }
        self.vertices += triangles.len();
    }

    pub fn text(&self) -> &str {
        &self.text
    }
}

fn min3(a: f64, b: f64, c: f64) -> f64 {
    let mut m = a;
    if b < m {
        m = b;
    }
    if c < m {
        m = c;
    }
    m
}

struct Region {
    id: String,
    style: String,
    args: Vec<f64>,
    axis: usize,
    outside: bool,
    origin: Vec3,
    a: Vec3,
    b: Vec3,
    c: Vec3,
}

impl Region {
    /// Signed distance, positive inside (approximate for cones and prisms).
    fn distance(&self, x: Vec3) -> Result<f64, String> {
        let args = &self.args;
        let mut d = 0.0;
        match self.style.as_str() {
            "block" => {
                d = HUGE;
                for k in 0..3 {
                    d = min3(d, x[k] - args[2 * k], args[2 * k + 1] - x[k]);
                }
            }
            "sphere" => d = args[3] - norm(sub(x, [args[0], args[1], args[2]])),
            "cylinder" | "cone" => {
                let ax = self.axis;
                let (u, v) = ((ax + 1) % 3, (ax + 2) % 3);
                let (first, second) = (u.min(v), u.max(v));
                let radial = (x[first] - args[0]).hypot(x[second] - args[1]);
                let (mut radius, mut low, mut high) = (args[2], args[3], args[4]);
                if self.style == "cone" {
                    low = args[4];
                    high = args[5];
                    let t = if high > low { ((x[ax] - low) / (high - low)).clamp(0.0, 1.0) } else { 0.0 };
                    radius = args[2] + t * (args[3] - args[2]);
                }
                d = min3(radius - radial, x[ax] - low, high - x[ax]);
            }
            "plane" => {
                let normal = [args[3], args[4], args[5]];
                d = dot(sub(x, [args[0], args[1], args[2]]), scale(normal, 1.0 / norm(normal)));
            }
            "prism" => {
                let cell = [self.a, self.b, self.c];
                let f = fractional(sub(x, self.origin), &cell)?;
                d = HUGE;
                for k in 0..3 {
                    let height = cell_volume(&cell) / norm(cross(cell[(k + 1) % 3], cell[(k + 2) % 3]));
                    d = min3(d, f[k] * height, (1.0 - f[k]) * height);
                }
            }
            _ => {}
        }
        Ok(if self.outside { -d } else { d })
    }
}

fn parse_bound(word: &str, lower: bool, edge: f64, fixed: &mut bool) -> Result<f64, String> {
    *fixed = true;
    if word == "INF" || word == "-INF" {
        return Ok(if lower { -HUGE } else { HUGE });
    }
    if word == "EDGE" {
        return Ok(edge);
    }
    *fixed = false;
    let (value, used) = af_core::format::stod(word).map_err(|_| "stod".to_string())?;
    if used != word.len() {
        return Err(format!("Not a number in a region command: {word}"));
    }
    Ok(value)
}

fn parse_regions(text: &str, lattice_scale: f64, low: Vec3, high: Vec3) -> Result<Vec<Region>, String> {
    let counts: BTreeMap<&str, usize> = [("block", 6), ("sphere", 4), ("cylinder", 6), ("cone", 7), ("plane", 6), ("prism", 9)].into_iter().collect();
    let mut regions = Vec::new();
    for raw in text.split('\n') {
        let line = match raw.find('#') {
            Some(hash) => &raw[..hash],
            None => raw,
        };
        let words: Vec<&str> = line.split_ascii_whitespace().collect();
        if words.is_empty() {
            continue;
        }
        if words[0] != "region" {
            return Err(format!("Each line must be a LAMMPS region command: {line}"));
        }
        if words.len() < 4 {
            return Err(format!("Incomplete region command: {line}"));
        }
        let mut r = Region { id: words[1].into(), style: words[2].into(), args: Vec::new(), axis: 0, outside: false, origin: [0.0; 3], a: [0.0; 3], b: [0.0; 3], c: [0.0; 3] };
        let Some(&count) = counts.get(r.style.as_str()) else {
            return Err(format!("Unsupported region style {} (block, sphere, cylinder, cone, plane, prism)", r.style));
        };
        let axial = r.style == "cylinder" || r.style == "cone";
        let mut k = 3;
        if axial {
            if words[3] != "x" && words[3] != "y" && words[3] != "z" {
                return Err(format!("The axis of a {} is x, y or z", r.style));
            }
            r.axis = (words[3].as_bytes()[0] - b'x') as usize;
            k += 1;
        }
        let needed = count - if axial { 1 } else { 0 };
        if words.len() < k + needed {
            return Err(format!("Too few values for region {}", r.id));
        }
        let mut fixed = vec![false; needed];
        for m in 0..needed {
            let (mut coordinate, mut lower) = (0usize, true);
            if (r.style == "block" || r.style == "prism") && m < 6 {
                coordinate = m / 2;
                lower = m % 2 == 0;
            } else if r.style == "cylinder" && m >= 3 {
                coordinate = r.axis;
                lower = m == 3;
            } else if r.style == "cone" && m >= 4 {
                coordinate = r.axis;
                lower = m == 4;
            }
            let mut is_fixed = false;
            r.args.push(parse_bound(words[k + m], lower, if lower { low[coordinate] } else { high[coordinate] }, &mut is_fixed)?);
            fixed[m] = is_fixed;
        }
        let mut lattice = true;
        let mut m = k + needed;
        while m < words.len() {
            if words[m] == "side" && m + 1 < words.len() {
                m += 1;
                r.outside = words[m] == "out";
            } else if words[m] == "units" && m + 1 < words.len() {
                m += 1;
                lattice = words[m] == "lattice";
            } else if words[m] == "open" && m + 1 < words.len() {
                m += 1;
            } else {
                return Err(format!("Unsupported region keyword {} (side, units and open are understood)", words[m]));
            }
            m += 1;
        }
        if lattice {
            for m in 0..r.args.len() {
                if !fixed[m] && !(r.style == "plane" && m >= 3) {
                    r.args[m] *= lattice_scale;
                }
            }
        }
        if r.style == "prism" {
            let a = &r.args;
            r.origin = [a[0], a[2], a[4]];
            r.a = [a[1] - a[0], 0.0, 0.0];
            r.b = [a[6], a[3] - a[2], 0.0];
            r.c = [a[7], a[8], a[5] - a[4]];
        }
        regions.push(r);
    }
    if regions.is_empty() {
        return Err("No region commands".into());
    }
    Ok(regions)
}

pub fn bond_orders(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    let n = positions.len();
    let species = finite_array(p.array("types")?, "types", 1)?.values;
    if species.len() != n {
        return Err("types needs one atomic number per atom".into());
    }
    let pbc = p.pbc("pbc", [false; 3])?;
    if pbc.iter().any(|&b| b) && !p.has("cell") {
        return Err("Periodic analysis requires a cell".into());
    }
    let cell = if p.has("cell") { cell_matrix(p.array("cell")?)? } else { identity() };
    let tolerance = positive(p.number_or("tolerance", 1.25)?, "tolerance")?;
    let covalent = af_core::elements::literature_covalent_radii();
    let radius = |z: f64| {
        let k = z.round() as i64 as usize;
        if k < covalent.len() && covalent[k] > 0.0 {
            covalent[k] as f64
        } else {
            1.5
        }
    };
    let mut largest = 0.0;
    for &z in &species {
        let r = radius(z);
        if largest < r {
            largest = r;
        }
    }
    let mut valence = vec![0.0; n];
    let mut coordination = vec![0i32; n];
    let mut bonds = Vec::new();
    let mut rounded: Vec<(f64, i64)> = Vec::new();
    for link in neighbor_list(&positions, &cell, pbc, 2.0 * largest * tolerance)? {
        if !first_of_pair(&link) {
            continue;
        }
        let r = norm(link.vector);
        let single = radius(species[link.i]) + radius(species[link.j]);
        if r > tolerance * single || r <= 0.0 {
            continue;
        }
        let order = 10.0f64.powf((single - r) / 0.6);
        let mut nearest = 1.0;
        for candidate in [1.0, 1.5, 2.0, 3.0] {
            if (order - candidate).abs() < (order - nearest).abs() {
                nearest = candidate;
            }
        }
        match rounded.iter_mut().find(|(k, _)| *k == nearest) {
            Some(entry) => entry.1 += 1,
            None => rounded.push((nearest, 1)),
        }
        for atom in [link.i, link.j] {
            valence[atom] += order;
            coordination[atom] += 1;
        }
        let mut bond = Json::object();
        bond.set("atoms", jv(&[link.i as f64, link.j as f64]));
        bond.set("length_A", r);
        bond.set("order", order);
        bond.set("rounded_order", nearest);
        bonds.push(bond);
    }
    rounded.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut counts = Json::object();
    for (order, count) in &rounded {
        counts.set(&format_g(*order, 6), *count as f64);
    }
    let mut result = Json::object();
    let bond_count = bonds.len();
    result.set("bonds", Json::Array(bonds));
    result.set("bond_count", bond_count as f64);
    result.set("rounded_order_counts", counts);
    result.set("valence_sum", jv(&valence));
    result.set("coordination", ji(&coordination));
    Ok(result)
}

pub fn lammps_regions(p: &Parameters) -> R {
    let positions = points(p.array("positions")?, "positions")?;
    if !p.has("regions") {
        return Err("Give the region commands in regions".into());
    }
    let text = match p.json("regions")? {
        Json::String(s) => s.replace(';', "\n"),
        _ => String::new(),
    };
    let (mut low, mut high) = ([HUGE; 3], [-HUGE; 3]);
    if p.has("cell") {
        let cell = cell_matrix(p.array("cell")?)?;
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let corner = row_times([i as f64, j as f64, k as f64], &cell);
                    for d in 0..3 {
                        if corner[d] < low[d] {
                            low[d] = corner[d];
                        }
                        if high[d] < corner[d] {
                            high[d] = corner[d];
                        }
                    }
                }
            }
        }
    } else {
        for x in &positions {
            for d in 0..3 {
                if x[d] - 1.0 < low[d] {
                    low[d] = x[d] - 1.0;
                }
                if high[d] < x[d] + 1.0 {
                    high[d] = x[d] + 1.0;
                }
            }
        }
        if positions.is_empty() {
            return Err("No atoms and no cell to define the domain".into());
        }
    }
    let regions = parse_regions(&text, positive(p.number_or("lattice_scale", 1.0)?, "lattice_scale")?, low, high)?;
    let (mut per_region, mut counts) = (Json::object(), Json::object());
    let mut any = vec![false; positions.len()];
    for r in &regions {
        let mut inside = vec![false; positions.len()];
        let mut count = 0i64;
        for (i, x) in positions.iter().enumerate() {
            inside[i] = r.distance(*x)? >= 0.0;
            if inside[i] {
                count += 1;
                any[i] = true;
            }
        }
        per_region.set(&r.id, Json::Array(inside.iter().map(|&b| Json::Bool(b)).collect()));
        counts.set(&r.id, count as f64);
    }
    let mut result = Json::object();
    result.set("inside", per_region);
    result.set("atoms_inside", counts);
    result.set("inside_any", Json::Array(any.iter().map(|&b| Json::Bool(b)).collect()));
    if p.boolean("export_mesh", true)? {
        let resolution = integer(p.number_or("resolution", 64.0)?, "resolution", 8)?;
        let mut longest = 0.0;
        for d in 0..3 {
            if longest < high[d] - low[d] {
                longest = high[d] - low[d];
            }
        }
        let step = longest / resolution as f64;
        let mut obj = ObjWriter::new("AtomForge LAMMPS regions");
        for r in &regions {
            let mut field = Grid { periodic: false, origin: DVec3::new(low[0] - step, low[1] - step, low[2] - step), ..Default::default() };
            let mut columns = [DVec3::ZERO; 3];
            for d in 0..3 {
                field.shape[d] = ((high[d] - low[d]) / step).ceil() as i32 + 3;
                columns[d][d] = step * (field.shape[d] - 1) as f64;
            }
            field.cell = DMat3::from_cols(columns[0], columns[1], columns[2]);
            field.values = vec![0.0; field.shape[0] as usize * field.shape[1] as usize * field.shape[2] as usize];
            for z in 0..field.shape[2] {
                for y in 0..field.shape[1] {
                    for x in 0..field.shape[0] {
                        let q = field.position(x, y, z);
                        let mut d = r.distance([q.x, q.y, q.z])?;
                        for k in 0..3 {
                            d = min3(d, q[k] - low[k], high[k] - q[k]);
                        }
                        let index = (z as usize * field.shape[1] as usize + y as usize) * field.shape[0] as usize + x as usize;
                        field.values[index] = d;
                    }
                }
            }
            obj.add(&r.id, &isosurface(&field, 0.0, None)?.vertices);
        }
        let mut files = Json::object();
        files.set("regions.obj", obj.text());
        result.set("files", files);
    }
    Ok(result)
}
