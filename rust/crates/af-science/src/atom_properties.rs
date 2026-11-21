//! Per-atom values of tool results, for colouring and selections (the C++
//! science/AtomProperties), with the viridis colour map.

use std::collections::BTreeMap;

use af_core::json::Json;

#[derive(Clone, Debug, PartialEq)]
pub struct AtomProperty {
    pub name: String,
    pub values: Vec<f64>,
}

pub mod property_name {
    pub const COORDINATION: &str = "Coordination";
    pub const VORONOI_COORDINATION: &str = "Voronoi coordination";
    pub const ORIENTATION_ANGLE: &str = "Orientation angle (degrees)";
    pub const PTM_ORDERING: &str = "PTM ordering (0 none, 1 pure, 2 L1_0, 3 L1_2 minority, 4 L1_2 majority, 5 B2, 6 zincblende/wurtzite, 7 other)";
}

fn numbers(value: Option<&Json>) -> Vec<f64> {
    match value {
        Some(Json::Array(items)) => items
            .iter()
            .map(|item| match item {
                Json::Number(n) => *n,
                Json::Bool(b) => {
                    if *b {
                        1.0
                    } else {
                        0.0
                    }
                }
                _ => f64::NAN,
            })
            .collect(),
        _ => Vec::new(),
    }
}

fn element(tensor: &Json, i: usize, j: usize) -> f64 {
    let Json::Array(rows) = tensor else { return f64::NAN };
    if rows.len() != 3 {
        return f64::NAN;
    }
    let Json::Array(row) = &rows[i] else { return f64::NAN };
    if row.len() != 3 {
        return f64::NAN;
    }
    match row[j] {
        Json::Number(n) => n,
        _ => f64::NAN,
    }
}

fn add(properties: &mut Vec<AtomProperty>, name: &str, values: Vec<f64>) {
    if !values.is_empty() {
        properties.push(AtomProperty { name: name.to_string(), values });
    }
}

type Builder = fn(&Json, &mut Vec<AtomProperty>);

const PTM_STRUCTURE: &str = "PTM structure (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico, 5 sc, 6 cubic diamond, 7 hex diamond)";
const STRUCTURE_TYPE: &str = "Structure type (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)";

fn local_strain(r: &Json, p: &mut Vec<AtomProperty>) {
    add(p, "D2min (A^2)", numbers(r.find("d2min_A2")));
    let (mut volumetric, mut shear) = (Vec::new(), Vec::new());
    if let Some(Json::Array(strains)) = r.find("green_lagrange_strain") {
        for e in strains {
            let (xx, yy, zz) = (element(e, 0, 0), element(e, 1, 1), element(e, 2, 2));
            let (xy, xz, yz) = (element(e, 0, 1), element(e, 0, 2), element(e, 1, 2));
            volumetric.push(xx + yy + zz);
            shear.push((xy * xy + xz * xz + yz * yz + ((xx - yy) * (xx - yy) + (yy - zz) * (yy - zz) + (xx - zz) * (xx - zz)) / 6.0).sqrt());
        }
    }
    add(p, "Von Mises shear strain", shear);
    add(p, "Volumetric strain trace(E)", volumetric);
    add(p, property_name::COORDINATION, numbers(r.find("coordination")));
}

fn bond_order(r: &Json, p: &mut Vec<AtomProperty>) {
    if let Some(Json::Object(members)) = r.find("order") {
        for (name, values) in members {
            add(p, name, numbers(Some(values)));
        }
    }
    add(p, property_name::COORDINATION, numbers(r.find("coordination")));
}

fn lammps_regions(r: &Json, p: &mut Vec<AtomProperty>) {
    add(p, "Inside a region (1 yes)", numbers(r.find("inside_any")));
    if let Some(Json::Object(members)) = r.find("inside") {
        for (name, values) in members {
            add(p, &format!("Inside region {name} (1 yes)"), numbers(Some(values)));
        }
    }
}

fn elastic_strain(r: &Json, p: &mut Vec<AtomProperty>) {
    add(p, "Volumetric elastic strain", numbers(r.find("volumetric_strain")));
    add(p, "Von Mises elastic shear strain", numbers(r.find("shear_strain")));
    let names = ["Elastic strain xx", "Elastic strain yy", "Elastic strain zz", "Elastic strain yz", "Elastic strain xz", "Elastic strain xy"];
    if let Some(Json::Array(rows)) = r.find("strain_voigt") {
        for (c, name) in names.iter().enumerate() {
            let values = rows
                .iter()
                .map(|row| match row {
                    Json::Array(items) if items.len() == 6 => match items[c] {
                        Json::Number(n) => n,
                        _ => f64::NAN,
                    },
                    _ => f64::NAN,
                })
                .collect();
            add(p, name, values);
        }
    }
}

fn builders() -> BTreeMap<&'static str, Builder> {
    BTreeMap::from([
        ("symmetry", (|r: &Json, p: &mut Vec<AtomProperty>| {
            add(p, "Symmetry orbit (equivalent atoms)", numbers(r.find("orbit_id")));
            add(p, "Wyckoff position (0 = a, 1 = b, ...)", numbers(r.find("wyckoff_index")));
        }) as Builder),
        ("cluster-analysis", |r, p| add(p, "Cluster id (-1 unselected)", numbers(r.find("cluster_id")))),
        ("void-analysis", |r, p| add(p, "Lining void id (-1 none)", numbers(r.find("lining_void_id")))),
        ("dislocation-lines", |r, p| {
            add(p, "Dislocation line id (-1 none)", numbers(r.find("line_id")));
            add(p, STRUCTURE_TYPE, numbers(r.find("structure_type")));
        }),
        ("structure-type", |r, p| add(p, STRUCTURE_TYPE, numbers(r.find("structure_type")))),
        ("diamond-structure", |r, p| {
            add(
                p,
                "Diamond type (0 other, 1 cubic, 2 cubic 1st nb, 3 cubic 2nd nb, 4 hexagonal, 5 hex 1st nb, 6 hex 2nd nb)",
                numbers(r.find("diamond_type")),
            )
        }),
        ("centrosymmetry", |r, p| add(p, "Centrosymmetry (A^2)", numbers(r.find("centrosymmetry_A2")))),
        ("local-strain", local_strain),
        ("bond-order", bond_order),
        ("wigner-seitz", |r, p| {
            add(p, "Distance to assigned site (A)", numbers(r.find("distance_A")));
            add(p, "Assigned site index", numbers(r.find("site_index")));
        }),
        ("ptm", |r, p| {
            add(p, PTM_STRUCTURE, numbers(r.find("structure_type")));
            add(p, "PTM RMSD", numbers(r.find("rmsd")));
            add(p, property_name::ORIENTATION_ANGLE, numbers(r.find("orientation_angle_deg")));
            add(p, "Local lattice constant (A)", numbers(r.find("lattice_constant_A")));
            add(p, property_name::PTM_ORDERING, numbers(r.find("ordering")));
        }),
        ("voronoi", |r, p| {
            add(p, "Atomic volume (A^3)", numbers(r.find("atomic_volume_A3")));
            add(p, property_name::VORONOI_COORDINATION, numbers(r.find("coordination")));
            add(p, "Max face order", numbers(r.find("max_face_order")));
        }),
        ("bond-orders", |r, p| {
            add(p, "Bond-valence sum", numbers(r.find("valence_sum")));
            add(p, "Bonded neighbours", numbers(r.find("coordination")));
        }),
        ("lammps-regions", lammps_regions),
        ("ambient-occlusion", |r, p| add(p, "Ambient brightness", numbers(r.find("brightness")))),
        ("surface-mesh", |r, p| add(p, "Surface atom (1 yes)", numbers(r.find("surface_atom")))),
        ("cna", |r, p| add(p, "CNA structure (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)", numbers(r.find("structure_type")))),
        ("vorotop", |r, p| {
            add(p, "VoroTop structure (0 other, 1 fcc, 2 bcc, 3 hcp with the built-in filter)", numbers(r.find("structure_type")));
            add(p, "Voronoi topology id", numbers(r.find("topology_id")));
        }),
        ("grain-segmentation", |r, p| {
            add(p, "Grain id (0 none)", numbers(r.find("grain_id")));
            add(p, PTM_STRUCTURE, numbers(r.find("structure_type")));
        }),
        ("elastic-strain", elastic_strain),
        ("fcc-planar-faults", |r, p| {
            add(p, "Planar fault (0 non-hcp, 1 stacking fault, 2 twin, 3 multilayer, 4 other hcp)", numbers(r.find("fault_type")))
        }),
        ("chill-plus", |r, p| {
            add(
                p,
                "Chill+ type (0 other, 1 hex ice, 2 cubic ice, 3 interfacial ice, 4 hydrate, 5 interfacial hydrate)",
                numbers(r.find("chill_type")),
            )
        }),
        ("local-entropy", |r, p| add(p, "Local entropy (kB)", numbers(r.find("entropy_kB")))),
        ("find-rings", |r, p| add(p, "Rings through atom", numbers(r.find("rings_per_atom")))),
        ("displacements", |r, p| add(p, "Displacement magnitude (A)", numbers(r.find("magnitude_A")))),
    ])
}

/// The per-atom properties of a tool's result.
pub fn per_atom_properties(tool: &str, result: &Json) -> Vec<AtomProperty> {
    let mut properties = Vec::new();
    if !matches!(result, Json::Object(_)) {
        return properties;
    }
    if let Some(builder) = builders().get(tool) {
        builder(result, &mut properties);
    }
    properties
}

/// The tools with per-atom properties (sorted).
pub fn tools_with_atom_properties() -> Vec<String> {
    builders().keys().map(|k| k.to_string()).collect()
}

fn clamp(x: f64, low: f64, high: f64) -> f64 {
    if x < low {
        low
    } else if high < x {
        high
    } else {
        x
    }
}

/// The viridis colour map (polynomial fit).
pub fn viridis(t: f64) -> [f32; 3] {
    let t = clamp(if t.is_finite() { t } else { 0.0 }, 0.0, 1.0);
    let r = 0.2777273272234177
        + t * (0.1050930431085774 + t * (-0.3308618287255563 + t * (-4.634230498983486 + t * (6.228269936347081 + t * (4.776384997670288 + t * -5.435455855934631)))));
    let g = 0.005407344544966578
        + t * (1.404613529898575 + t * (0.214847559468213 + t * (-5.799100973351585 + t * (14.17993336680509 + t * (-13.74514537774601 + t * 4.645852612178535)))));
    let b = 0.3340998053353061
        + t * (1.384590162594685 + t * (0.09509516302823659 + t * (-19.33244095627987 + t * (56.69055260068105 + t * (-65.35303263337234 + t * 26.3124352495832)))));
    [clamp(r, 0.0, 1.0) as f32, clamp(g, 0.0, 1.0) as f32, clamp(b, 0.0, 1.0) as f32]
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PropertyRange {
    pub low: f64,
    pub high: f64,
}

impl Default for PropertyRange {
    fn default() -> Self {
        Self { low: 0.0, high: 1.0 }
    }
}

/// The range of the finite values ([0, 1] when there are none).
pub fn finite_range(values: &[f64]) -> PropertyRange {
    let mut range = PropertyRange { low: f64::INFINITY, high: f64::NEG_INFINITY };
    for &v in values {
        if v.is_finite() {
            if v < range.low {
                range.low = v;
            }
            if range.high < v {
                range.high = v;
            }
        }
    }
    if !(range.low <= range.high) {
        return PropertyRange::default();
    }
    range
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PropertyDisplay {
    pub auto_range: bool,
    pub range: PropertyRange,
    pub hide_outside: bool,
    pub hide_invalid: bool,
}

impl Default for PropertyDisplay {
    fn default() -> Self {
        Self { auto_range: true, range: PropertyRange::default(), hide_outside: false, hide_invalid: false }
    }
}

/// Colours (viridis over the range) and visibility of atoms by a property.
pub fn colour_by_property(values: &[f64], display: &PropertyDisplay) -> (Vec<[f32; 3]>, Vec<bool>) {
    let range = if display.auto_range { finite_range(values) } else { display.range };
    let span = range.high - range.low;
    let mut colours = vec![[0.55f32; 3]; values.len()];
    let mut visible = vec![true; values.len()];
    for (i, &v) in values.iter().enumerate() {
        if !v.is_finite() {
            visible[i] = !display.hide_invalid;
            continue;
        }
        colours[i] = viridis(if span > 0.0 { (v - range.low) / span } else { 0.5 });
        if display.hide_outside && (v < range.low || v > range.high) {
            visible[i] = false;
        }
    }
    (colours, visible)
}
