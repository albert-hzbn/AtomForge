//! Analysis > Electronic Structure > Electronic Post-Processing (the C++
//! ui/ElectronicPostProcessing with ui/ElectronicSliceViewport and the
//! graphics/ElectronicViewport previews): volumetric data from VASP, Cube and
//! XSF files, 38 calculations in three tool groups (charge transfer, general
//! analysis, bonding and topology), results as new fields or tables, an
//! isosurface or volume preview in 3D, a 2D section with contours, exports
//! (grids, tables, meshes, PNG and SVG), workspaces and calculation undo.
//!
//! The 3D previews are rendered on the CPU here (the C++ uses OpenGL): the
//! volume by the same ray march as its fragment shader, the isosurface as
//! depth-sorted, Gouraud-lit triangles with the same lighting.

use std::collections::{BTreeMap, VecDeque};
use std::sync::mpsc::{channel, Receiver};

use af_core::format::format_g;
use af_core::glm::{inverse_d, mul_d};
use af_electronic::bader::bader_on_grid;
use af_electronic::charge::{apply_mask, boolean_mask, charge_summary, cumulative_charge, density_difference, split_density, threshold_mask};
use af_electronic::density::{integrate, integrate_sphere, peaks, voronoi_integrate};
use af_electronic::field::{contours, dori, energy_density, gradient, laplacian, line_profile, macroscopic_average, planar_average, reduced_density_gradient, section, signed_density, smooth};
use af_electronic::fourier::{ewald, fourier_synthesis, patterson, structure_factors, Complex, Reflection};
use af_electronic::project::{load_project, save_project, Workspace, FORMAT};
use af_electronic::topology::betti_curve;
use af_electronic::volume::{arithmetic, periodic_endpoints, resample, scale, Site};
use af_electronic::{isosurface, load_volume, save_mesh, save_volume, Grid, Mesh, Volume};
use eframe::egui;
use glam::{DVec2, DVec3, IVec3, Mat3, Vec3};

use crate::analysis_panels::{field, file_name, line_plot, section as heading};
use crate::electronic_slice::{axis_plane, clip_segment_to_cell, clip_to_cell, display_range, palette_colour, sample_slice_plane, slice_plane, SlicePlane};
use crate::file_browser::{FileBrowser, Purpose};
use crate::plot_widget::unit_display;
use crate::structure_dialog::{Input, Outcome, StructureDialog};

const GENERAL_TOOLS: [&str; 25] = [
    "Total integral",
    "Add reference",
    "Subtract reference (difference density)",
    "Multiply reference",
    "Divide by reference",
    "Scale",
    "Gaussian smoothing",
    "Cartesian gradient",
    "Laplacian",
    "Energy-density conversion",
    "Line profile",
    "Planar average",
    "Macroscopic average",
    "Plane section",
    "Contour segments (z=0)",
    "Peak search",
    "Voronoi site integration",
    "Sphere integration",
    "Structure factors",
    "Fourier synthesis",
    "Patterson density",
    "Ewald site potentials",
    "Isosurface",
    "Resample onto reference",
    "Verify periodic endpoint planes",
];
const CHARGE_TOOLS: [&str; 8] = ["Density difference", "Threshold mask", "Boolean masks", "Apply mask", "Accumulation / depletion", "Charge summary", "Cumulative profile", "Invert mask"];
const BONDING_TOOLS: [&str; 5] = ["Reduced density gradient", "Signed density (NCI)", "DORI", "Betti curve", "Bader partition"];
const QUANTITIES: [&str; 5] = ["Auto (Cube/XSF: raw)", "Density", "Potential", "ELF", "Raw"];
const QUANTITY_NAMES: [&str; 5] = ["auto", "density", "potential", "elf", "raw"];
const EXPORT_FORMATS: [&str; 8] = ["Grid: XSF", "Grid: Cube", "Grid: VASP", "Table: CSV", "Surface: OBJ", "Surface: PLY", "View: 3D PNG", "View: 2D SVG"];
const EXPORT_EXTENSIONS: [&str; 8] = [".xsf", ".cube", ".vasp", ".csv", ".obj", ".ply", ".png", ".svg"];
const AXES: [&str; 3] = ["a*", "b*", "c*"];
const BOOLEAN_OPERATIONS: [&str; 4] = ["Union", "Intersection", "Difference (field minus operand)", "XOR"];
const BOOLEAN_NAMES: [&str; 4] = ["union", "intersection", "difference", "xor"];
const QUALITIES: [&str; 4] = ["Fast", "Balanced", "High", "Very high"];
const PALETTES: [&str; 3] = ["Spectrum", "Blue-white-red", "Sequential blue"];
const SLICE_RESOLUTIONS: [i32; 4] = [65, 129, 257, 513];
const VOLUME_RESOLUTIONS: [i32; 4] = [64, 96, 192, 256];
const RAY_STEPS: [i32; 4] = [128, 256, 512, 1024];
const BACKGROUND: [f32; 3] = [0.96, 0.97, 0.985];
const MESH_LIMIT: usize = 12_000_000;

/// A calculation's result (the C++ Output).
#[derive(Clone, Debug, Default)]
pub struct Output {
    pub volume: Volume,
    pub table: Vec<f64>,
    pub columns: i32,
    pub heading: String,
    pub mesh: Mesh,
    pub loaded: bool,
    pub reference_loaded: bool,
    pub append_mesh: bool,
    pub surface_ready: bool,
    pub colour_low: f64,
    pub colour_high: f64,
    pub colour_unit: String,
    pub source_path: String,
}

/// Everything one calculation needs (copied, so it can run in the background).
#[derive(Clone)]
struct Job {
    op: i32,
    source: Grid,
    reference: Grid,
    sites: Vec<Site>,
    count: i32,
    axis: i32,
    window: i32,
    radius: i32,
    mask_low: f64,
    mask_high: f64,
    reference_weight: f64,
    boolean_operation: usize,
    scalar: f64,
    sigma: f64,
    alpha: f64,
    real_cutoff: f64,
    reciprocal_cutoff: f64,
    start: DVec3,
    end: DVec3,
    u: DVec3,
    v: DVec3,
    reflections: String,
    charges: String,
    thresholds: String,
    bonding_floor: f64,
    coloured: bool,
    append_surface: bool,
    existing_surface: bool,
    existing_unit: String,
}

fn parse_reflections(text: &str, synthesis: bool) -> Result<(Vec<IVec3>, Vec<Reflection>), String> {
    let (mut indices, mut reflections) = (Vec::new(), Vec::new());
    for row in text.lines() {
        if row.trim_matches(|c| c == ' ' || c == '\t' || c == '\r').is_empty() {
            continue;
        }
        let mut words = row.split_whitespace();
        let mut int = || words.next().and_then(|w| w.parse::<i32>().ok());
        let h = match (int(), int(), int()) {
            (Some(a), Some(b), Some(c)) => IVec3::new(a, b, c),
            _ => return Err("Malformed reflection row".into()),
        };
        let (mut re, mut im) = (0.0, 0.0);
        if synthesis {
            match (words.next().and_then(|w| w.parse::<f64>().ok()), words.next().and_then(|w| w.parse::<f64>().ok())) {
                (Some(a), Some(b)) => (re, im) = (a, b),
                _ => return Err("Malformed reflection row".into()),
            }
        }
        if words.next().is_some() {
            return Err("Extra reflection columns".into());
        }
        indices.push(h);
        reflections.push(Reflection { hkl: h, value: Complex { re, im } });
    }
    Ok((indices, reflections))
}

fn run(job: Job) -> Result<Output, String> {
    let op = job.op;
    let source = &job.source;
    let reference = &job.reference;
    if op == 22 && !display_range(&source.values)?.contains(job.scalar) {
        return Err("3D isovalue is outside the selected field's range.".into());
    }
    let mut out = Output { volume: Volume { sites: job.sites.clone(), fields: Vec::new() }, ..Default::default() };
    let profile = |out: &mut Output, rows: Vec<DVec2>| {
        out.columns = 2;
        out.heading = "distance_A,value".into();
        for r in rows {
            out.table.extend([r.x, r.y]);
        }
    };
    match op {
        0 => {
            out.columns = 1;
            out.heading = "integral (field unit * A^3)".into();
            out.table = vec![integrate(source)?];
        }
        1..=4 => out.volume.fields.push(arithmetic(source, reference, ["add", "subtract", "multiply", "divide"][op as usize - 1])?),
        5 => out.volume.fields.push(scale(source, job.scalar)?),
        25 => out.volume.fields.push(density_difference(source, reference, job.reference_weight)?),
        26 => out.volume.fields.push(threshold_mask(source, job.mask_low, job.mask_high)?),
        27 => out.volume.fields.push(boolean_mask(source, reference, BOOLEAN_NAMES[job.boolean_operation])?),
        28 => out.volume.fields.push(apply_mask(source, reference)?),
        29 => out.volume.fields.extend(split_density(source)?),
        30 => {
            out.columns = 3;
            out.heading = "accumulation_e,depletion_e,net_e".into();
            out.table = charge_summary(source)?.to_vec();
        }
        31 => {
            profile(&mut out, cumulative_charge(source, job.axis)?);
            out.heading = "distance_A,cumulative_e".into();
        }
        32 => out.volume.fields.push(boolean_mask(&threshold_mask(source, 0.0, 1.0)?, source, "difference")?),
        6 => out.volume.fields.push(smooth(source, job.sigma, job.radius)?),
        7 => out.volume.fields.extend(gradient(source)?),
        8 => out.volume.fields.push(laplacian(source)?),
        9 => out.volume.fields.extend(energy_density(source, 1e-12)?),
        10 => profile(&mut out, line_profile(source, job.start, job.end, job.count)?),
        11 => profile(&mut out, planar_average(source, job.axis)?),
        12 => profile(&mut out, macroscopic_average(source, job.axis, job.window)?),
        13 => out.volume.fields.push(section(source, job.start, job.u, job.v, job.count, job.count)?),
        14 => {
            out.columns = 3;
            out.heading = "x_A,y_A,z_A (successive pairs are segment endpoints)".into();
            for p in contours(source, job.scalar)? {
                out.table.extend([p.x, p.y, p.z]);
            }
        }
        15 => {
            out.columns = 4;
            out.heading = "x_A,y_A,z_A,value".into();
            for p in peaks(source, 0)? {
                out.table.extend([p.position.x, p.position.y, p.position.z, p.value]);
            }
        }
        16 => {
            let positions: Vec<DVec3> = job.sites.iter().map(|s| s.position).collect();
            out.columns = 2;
            out.heading = "integral,volume_A^3 (atom order)".into();
            for b in voronoi_integrate(source, &positions)? {
                out.table.extend([b.integral, b.volume]);
            }
        }
        17 => {
            out.columns = 1;
            out.heading = "sphere_integral".into();
            out.table = vec![integrate_sphere(source, job.start, job.scalar)?];
        }
        18 | 19 => {
            let (indices, reflections) = parse_reflections(&job.reflections, op == 19)?;
            if op == 19 {
                out.volume.fields.push(fourier_synthesis(source, &reflections)?);
            } else {
                out.columns = 5;
                out.heading = "h,k,l,real,imag".into();
                for r in structure_factors(source, &indices)? {
                    out.table.extend([r.hkl.x as f64, r.hkl.y as f64, r.hkl.z as f64, r.value.re, r.value.im]);
                }
            }
        }
        20 => out.volume.fields.push(patterson(source)?),
        21 => {
            if !source.periodic {
                return Err("Ewald requires periodic geometry".into());
            }
            let mut charges = Vec::new();
            for w in job.charges.split_whitespace() {
                charges.push(w.parse::<f64>().map_err(|_| "Invalid charge list".to_string())?);
            }
            let positions: Vec<DVec3> = job.sites.iter().map(|s| s.position).collect();
            let e = ewald(&source.cell, &positions, &charges, job.alpha, job.real_cutoff, job.reciprocal_cutoff)?;
            out.columns = 1;
            out.heading = "cell_energy_eV_then_site_potentials_V".into();
            out.table = vec![e.energy];
            out.table.extend(e.potentials);
        }
        22 => {
            let colour_grid = if job.coloured { reference } else { source };
            if job.append_surface && job.existing_surface && job.existing_unit != colour_grid.unit {
                return Err("Appended surfaces must use the same coloring units.".into());
            }
            out.mesh = isosurface(source, job.scalar, if job.coloured { Some(reference) } else { None })?;
            out.append_mesh = job.append_surface;
            out.surface_ready = true;
            out.colour_low = colour_grid.values.iter().copied().fold(f64::INFINITY, f64::min);
            out.colour_high = colour_grid.values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            out.colour_unit = colour_grid.unit.clone();
        }
        23 => out.volume.fields.push(resample(source, reference)?),
        24 => out.volume.fields.push(periodic_endpoints(source, 1e-8)?),
        33 => out.volume.fields.push(reduced_density_gradient(source, job.bonding_floor)?),
        34 => out.volume.fields.push(signed_density(source, job.bonding_floor)?),
        35 => out.volume.fields.push(dori(source, job.bonding_floor)?),
        36 => {
            let thresholds: Vec<f64> = job.thresholds.split_whitespace().map_while(|w| w.parse().ok()).collect();
            if thresholds.is_empty() {
                return Err("Supply at least one threshold".into());
            }
            out.columns = 4;
            out.heading = "threshold,betti0,betti1,betti2".into();
            for r in betti_curve(source, &thresholds)? {
                out.table.extend([r.x, r.y, r.z, r.w]);
            }
        }
        37 => {
            let partition = bader_on_grid(source)?;
            let mut basins = source.clone();
            basins.name = "Bader basin index".into();
            basins.unit = "basin index".into();
            for (v, &b) in basins.values.iter_mut().zip(&partition.basin) {
                *v = b as f64;
            }
            basins.validate()?;
            out.volume.fields.push(basins);
            out.columns = 6;
            out.heading = "basin,charge_e,volume_A3,max_x_A,max_y_A,max_z_A".into();
            for (i, b) in partition.basins.iter().enumerate() {
                out.table.extend([i as f64, b.charge, b.volume, b.maximum.x, b.maximum.y, b.maximum.z]);
            }
        }
        _ => return Err(format!("Unknown operation {op}")),
    }
    Ok(out)
}

/// The 2D section view's own settings (the C++ ElectronicSliceViewport).
struct SliceView {
    axis: i32,
    position: f32,
    point: [f32; 3],
    normal: [f32; 3],
    level: f32,
    zoom: f32,
    angle: f32,
    pan: [f32; 2],
    dirty: bool,
    show_contour: bool,
    estimate_level: bool,
    resolution: i32,
    slice: Option<Grid>,
    visible: Vec<f64>,
    contours: Vec<DVec3>,
}

impl Default for SliceView {
    fn default() -> Self {
        Self {
            axis: 2,
            position: 0.5,
            point: [0.5; 3],
            normal: [1.0; 3],
            level: 0.1,
            zoom: 1.0,
            angle: 0.0,
            pan: [0.0; 2],
            dirty: true,
            show_contour: false,
            estimate_level: true,
            resolution: 129,
            slice: None,
            visible: Vec::new(),
            contours: Vec::new(),
        }
    }
}

impl SliceView {
    fn invalidate(&mut self) {
        self.dirty = true;
        self.estimate_level = true;
    }

    fn reset(&mut self) {
        self.zoom = 1.0;
        self.angle = 0.0;
        self.pan = [0.0; 2];
    }

    fn set_resolution(&mut self, r: i32) {
        if self.resolution != r {
            self.resolution = r;
            self.dirty = true;
        }
    }

    fn plane(&self, grid: &Grid) -> Result<SlicePlane, String> {
        if self.axis == 3 {
            let p = self.point.map(|v| v as f64);
            let n = self.normal.map(|v| v as f64);
            return slice_plane(grid, DVec3::from_array(p), DVec3::from_array(n));
        }
        Ok(axis_plane(grid, self.axis as usize, self.position as f64))
    }

    fn settings(&self) -> BTreeMap<String, f64> {
        let v = [
            ("axis", self.axis as f64),
            ("position", self.position as f64),
            ("px", self.point[0] as f64),
            ("py", self.point[1] as f64),
            ("pz", self.point[2] as f64),
            ("nx", self.normal[0] as f64),
            ("ny", self.normal[1] as f64),
            ("nz", self.normal[2] as f64),
            ("level", self.level as f64),
            ("zoom", self.zoom as f64),
            ("angle", self.angle as f64),
            ("panx", self.pan[0] as f64),
            ("pany", self.pan[1] as f64),
            ("contour", self.show_contour as i32 as f64),
        ];
        v.into_iter().map(|(k, x)| (k.to_string(), x)).collect()
    }

    fn restore(&mut self, values: &BTreeMap<String, f64>) {
        let read = |name: &str, fallback: f64| values.get(name).copied().unwrap_or(fallback);
        self.axis = read("axis", 2.0).clamp(0.0, 3.0) as i32;
        self.position = read("position", 0.5).clamp(0.0, 1.0) as f32;
        self.point = [read("px", 0.5) as f32, read("py", 0.5) as f32, read("pz", 0.5) as f32];
        self.normal = [read("nx", 1.0) as f32, read("ny", 1.0) as f32, read("nz", 1.0) as f32];
        self.level = read("level", 0.1) as f32;
        self.zoom = read("zoom", 1.0).clamp(0.1, 100.0) as f32;
        self.angle = read("angle", 0.0) as f32;
        self.pan = [read("panx", 0.0) as f32, read("pany", 0.0) as f32];
        self.show_contour = read("contour", 0.0) != 0.0;
        self.dirty = true;
        self.estimate_level = false;
    }

    /// Resamples the section when its settings changed (bounded preview resolution).
    fn update(&mut self, grid: &Grid) -> Result<SlicePlane, String> {
        let geometry = self.plane(grid)?;
        if geometry.boundary.len() < 3 {
            return Ok(geometry);
        }
        if self.dirty || self.slice.is_none() {
            let a = if self.axis == 0 { 1 } else { 0 };
            let b = if self.axis == 2 { 1 } else { 2 };
            let resolution = self.resolution.clamp(33, 513);
            let extra = if grid.periodic { 1 } else { 0 };
            let nu = if self.axis == 3 { resolution } else { (grid.shape[a] + extra).clamp(2, resolution) };
            let nv = if self.axis == 3 { resolution } else { (grid.shape[b] + extra).clamp(2, resolution) };
            let slice = sample_slice_plane(grid, &geometry, nu, nv)?;
            self.visible.clear();
            let inverse = inverse_d(&grid.cell);
            for y in 0..nv {
                for x in 0..nu {
                    let f = mul_d(&inverse, slice.position(x, y, 0) - grid.origin);
                    if self.axis != 3 || (f.cmpge(DVec3::splat(-1e-9)).all() && f.cmple(DVec3::splat(1.0 + 1e-9)).all()) {
                        self.visible.push(slice.values[slice.index(x, y, 0)]);
                    }
                }
            }
            for &p in &geometry.boundary {
                self.visible.push(grid.sample(p)?);
            }
            self.slice = Some(slice);
        }
        if self.estimate_level {
            self.level = display_range(&self.visible)?.suggested as f32;
            self.estimate_level = false;
        }
        Ok(geometry)
    }

    fn recontour(&mut self) -> Result<(), String> {
        if let Some(s) = &self.slice {
            self.contours = contours(s, self.level as f64)?;
        }
        self.dirty = false;
        Ok(())
    }

    /// The section's screen mapping: an orthonormal basis of the plane, rotated, zoomed and panned.
    fn screen(&self, slice: &Grid, size: egui::Vec2) -> impl Fn(DVec3) -> egui::Vec2 {
        let ex = slice.cell.x_axis.normalize();
        let ey = (slice.cell.y_axis - ex * ex.dot(slice.cell.y_axis)).normalize();
        let center = slice.origin + (slice.cell.x_axis + slice.cell.y_axis) * 0.5;
        let extent = (slice.cell.x_axis + slice.cell.y_axis).length().max((slice.cell.x_axis - slice.cell.y_axis).length());
        let scale = 0.85 * size.x.min(size.y) * self.zoom / extent as f32;
        let angle = self.angle * std::f32::consts::PI / 180.0;
        let (cs, sn) = (angle.cos(), angle.sin());
        let pan = self.pan;
        move |point: DVec3| {
            let p = point - center;
            let (x, y) = (p.dot(ex) as f32, p.dot(ey) as f32);
            egui::vec2(size.x * (0.5 + pan[0]) + scale * (cs * x - sn * y), size.y * (0.5 + pan[1]) - scale * (sn * x + cs * y))
        }
    }

    /// The coloured cells (clipped to the cell for arbitrary planes): (polygon, value).
    fn cells(&self, grid: &Grid) -> Vec<(Vec<DVec3>, f64)> {
        let Some(s) = &self.slice else { return Vec::new() };
        let mut out = Vec::new();
        for y in 0..s.shape[1] - 1 {
            for x in 0..s.shape[0] - 1 {
                let value = (s.values[s.index(x, y, 0)] + s.values[s.index(x + 1, y, 0)] + s.values[s.index(x, y + 1, 0)] + s.values[s.index(x + 1, y + 1, 0)]) * 0.25;
                let mut polygon = vec![s.position(x, y, 0), s.position(x, y + 1, 0), s.position(x + 1, y + 1, 0), s.position(x + 1, y, 0)];
                if self.axis == 3 {
                    polygon = clip_to_cell(grid, polygon);
                }
                if polygon.len() >= 3 {
                    out.push((polygon, value));
                }
            }
        }
        out
    }

    fn segments(&self, grid: &Grid) -> Vec<(DVec3, DVec3)> {
        let mut out = Vec::new();
        for pair in self.contours.chunks_exact(2) {
            if self.axis == 3 {
                if let Some(s) = clip_segment_to_cell(grid, pair[0], pair[1]) {
                    out.push(s);
                }
            } else {
                out.push((pair[0], pair[1]));
            }
        }
        out
    }

    /// The section as SVG (1200 x 1200).
    fn svg(&self, grid: &Grid, low: f32, high: f32, palette: i32) -> Result<String, String> {
        let geometry = self.plane(grid)?;
        if geometry.boundary.len() < 3 {
            return Err("The section does not intersect the cell".into());
        }
        let mut view = SliceView { slice: Some(sample_slice_plane(grid, &geometry, self.resolution, self.resolution)?), ..self.clone_settings() };
        let slice = view.slice.clone().unwrap();
        let screen = self.screen(&slice, egui::vec2(1200.0, 1200.0));
        let num = |v: f32| format_g(v as f64, 9);
        let mut out = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1200\" height=\"1200\" viewBox=\"0 0 1200 1200\">\n<title>AtomForge density section</title><defs><clipPath id=\"frame\"><rect width=\"1200\" height=\"1200\"/></clipPath></defs><rect width=\"1200\" height=\"1200\" fill=\"#f5f7fb\"/><g clip-path=\"url(#frame)\">\n");
        for (polygon, value) in view.cells(grid) {
            let t = if high > low { ((value - low as f64) / (high as f64 - low as f64)) as f32 } else { 0.5 };
            let c = palette_colour(t, palette);
            out += &format!("<polygon fill=\"rgb({},{},{})\" points=\"", (c.x * 255.0) as i32, (c.y * 255.0) as i32, (c.z * 255.0) as i32);
            for p in polygon {
                let q = screen(p);
                out += &format!("{},{} ", num(q.x), num(q.y));
            }
            out += "\"/>\n";
        }
        if self.show_contour {
            view.contours = contours(&slice, self.level as f64)?;
            for (a, b) in view.segments(grid) {
                let (a, b) = (screen(a), screen(b));
                out += &format!("<path stroke=\"#1e2837\" fill=\"none\" d=\"M{},{} L{},{}\"/>\n", num(a.x), num(a.y), num(b.x), num(b.y));
            }
        }
        out += "</g></svg>\n";
        Ok(out)
    }

    fn clone_settings(&self) -> SliceView {
        let mut v = SliceView::default();
        v.restore(&self.settings());
        v.resolution = self.resolution;
        v
    }
}

/// The 3D preview (the C++ ElectronicViewport, rendered on the CPU).
#[derive(Default)]
struct View3D {
    quality: i32,
    /// The field resampled onto at most `VOLUME_RESOLUTIONS[quality]` points per axis.
    volume: Vec<f32>,
    volume_shape: [usize; 3],
    inverse_cell: Mat3,
    volume_center: DVec3,
    volume_radius: f64,
    /// Mesh vertices (normalised positions), normals and scalars; three per triangle.
    positions: Vec<Vec3>,
    normals: Vec<Vec3>,
    scalars: Vec<f32>,
    mesh_center: DVec3,
    mesh_radius: f64,
    /// Back-to-front triangle order for (yaw, pitch).
    order: Vec<u32>,
    sorted_for: Option<(f32, f32)>,
    texture: Option<(egui::TextureHandle, String)>,
}

#[derive(Clone, Copy)]
struct Camera {
    yaw: f32,
    pitch: f32,
    zoom: f32,
    pan: [f32; 2],
}

impl Camera {
    fn rotation(&self) -> Mat3 {
        Mat3::from_rotation_x(self.pitch) * Mat3::from_rotation_y(self.yaw)
    }
    fn extent(&self) -> f32 {
        1.25 / self.zoom.max(0.1)
    }
}

fn ramp(t: f32, palette: i32) -> Vec3 {
    palette_colour(t, palette)
}

impl View3D {
    fn set_volume(&mut self, grid: &Grid) -> Result<(), String> {
        grid.validate()?;
        let resolution = VOLUME_RESOLUTIONS[self.quality.clamp(0, 3) as usize];
        let extra = if grid.periodic { 1 } else { 0 };
        let shape = [0, 1, 2].map(|i| (grid.shape[i] + extra).clamp(2, resolution) as usize);
        let mut values = vec![0.0f32; shape[0] * shape[1] * shape[2]];
        for z in 0..shape[2] {
            for y in 0..shape[1] {
                for x in 0..shape[0] {
                    let f = DVec3::new(x as f64 / (shape[0] - 1) as f64, y as f64 / (shape[1] - 1) as f64, z as f64 / (shape[2] - 1) as f64);
                    values[(z * shape[1] + y) * shape[0] + x] = grid.sample(grid.origin + mul_d(&grid.cell, f))? as f32;
                }
            }
        }
        let mut radius = 0.0f64;
        for z in [-1.0, 1.0] {
            for y in [-1.0, 1.0] {
                for x in [-1.0, 1.0] {
                    radius = radius.max(mul_d(&grid.cell, DVec3::new(x, y, z) * 0.5).length());
                }
            }
        }
        let inverse = inverse_d(&(grid.cell / radius));
        self.inverse_cell = Mat3::from_cols(inverse.x_axis.as_vec3(), inverse.y_axis.as_vec3(), inverse.z_axis.as_vec3());
        self.volume_center = grid.origin + mul_d(&grid.cell, DVec3::splat(0.5));
        self.volume_radius = radius;
        self.volume = values;
        self.volume_shape = shape;
        self.texture = None;
        Ok(())
    }

    fn set_mesh(&mut self, mesh: &Mesh) {
        self.positions.clear();
        self.normals.clear();
        self.scalars.clear();
        self.sorted_for = None;
        self.texture = None;
        if mesh.vertices.is_empty() {
            return;
        }
        let (mut low, mut high) = (mesh.vertices[0], mesh.vertices[0]);
        for &p in &mesh.vertices {
            low = low.min(p);
            high = high.max(p);
        }
        let center = (low + high) * 0.5;
        let radius = 1e-12f64.max((high - low).length() * 0.5);
        self.mesh_center = center;
        self.mesh_radius = radius;
        self.positions = mesh.vertices.iter().map(|&p| ((p - center) / radius).as_vec3()).collect();
        self.scalars = mesh.colors.iter().map(|&c| c as f32).collect();
        let key = |p: Vec3| ((p.x as f64 * 1e6).round() as i64, (p.y as f64 * 1e6).round() as i64, (p.z as f64 * 1e6).round() as i64);
        let mut normals: std::collections::HashMap<(i64, i64, i64), Vec3> = std::collections::HashMap::new();
        for t in self.positions.chunks_exact(3) {
            let n = (t[1] - t[0]).cross(t[2] - t[0]);
            for &p in t {
                *normals.entry(key(p)).or_insert(Vec3::ZERO) += n;
            }
        }
        self.normals = self
            .positions
            .iter()
            .map(|&p| {
                let n = normals[&key(p)];
                if n.length() > 1e-12 {
                    n.normalize()
                } else {
                    Vec3::Z
                }
            })
            .collect();
    }

    fn sort(&mut self, camera: Camera) {
        if self.sorted_for == Some((camera.yaw, camera.pitch)) {
            return;
        }
        let r = camera.rotation();
        let n = self.positions.len() / 3;
        let depths: Vec<f32> = (0..n).map(|i| (r * ((self.positions[i * 3] + self.positions[i * 3 + 1] + self.positions[i * 3 + 2]) / 3.0)).z).collect();
        let mut order: Vec<u32> = (0..n as u32).collect();
        order.sort_by(|&a, &b| depths[a as usize].total_cmp(&depths[b as usize]));
        self.order = order;
        self.sorted_for = Some((camera.yaw, camera.pitch));
    }

    /// Normalised device coordinates of a point (ElectronicViewport::project).
    fn project(&self, point: DVec3, volume: bool, aspect: f32, camera: Camera) -> egui::Vec2 {
        let (center, radius) = if volume { (self.volume_center, self.volume_radius) } else { (self.mesh_center, self.mesh_radius) };
        let p = ((point - center) / radius).as_vec3();
        let eye = camera.rotation() * p;
        let extent = camera.extent();
        egui::vec2((eye.x + camera.pan[0]) / (extent * aspect), (eye.y + camera.pan[1]) / extent)
    }

    /// The lit, depth-sorted surface triangles in normalised device
    /// coordinates, with linear RGBA vertex colours.
    fn surface_triangles(&mut self, aspect: f32, camera: Camera, opacity: f32, low: f32, high: f32, palette: i32, specular: f32, shininess: f32) -> Vec<[(egui::Vec2, [f32; 4]); 3]> {
        self.sort(camera);
        let r = camera.rotation();
        let extent = camera.extent();
        let key_light = Vec3::new(-0.4, 0.6, 1.0).normalize();
        let fill_light = Vec3::new(0.7, 0.1, 0.6).normalize();
        let shade = |i: usize| -> (egui::Vec2, [f32; 4]) {
            let eye = r * self.positions[i] + Vec3::new(camera.pan[0], camera.pan[1], -3.0);
            let mut n = (r * self.normals[i]).normalize();
            let v = (-eye).normalize();
            if n.dot(v) < 0.0 {
                n = -n;
            }
            let t = if high > low { ((self.scalars[i] - low) / (high - low)).clamp(0.0, 1.0) } else { 0.5 };
            let diffuse = 0.27 + 0.57 * n.dot(key_light).max(0.0) + 0.16 * n.dot(fill_light).max(0.0);
            let highlight = specular * n.dot((key_light + v).normalize()).max(0.0).powf(shininess);
            let c = (ramp(t, palette) * diffuse + Vec3::splat(highlight)).clamp(Vec3::ZERO, Vec3::ONE);
            (egui::vec2(eye.x / (extent * aspect), eye.y / extent), [c.x, c.y, c.z, opacity])
        };
        self.order.iter().map(|&t| [shade(t as usize * 3), shade(t as usize * 3 + 1), shade(t as usize * 3 + 2)]).collect()
    }

    /// The ray-marched volume (the C++ volume fragment shader), RGBA rows top first.
    fn render_volume(&self, width: usize, height: usize, camera: Camera, opacity: f32, threshold: f32, low: f32, high: f32, palette: i32) -> Vec<u8> {
        let steps = RAY_STEPS[self.quality.clamp(0, 3) as usize];
        let inverse_rotation = camera.rotation().transpose();
        let (aspect, extent) = (width as f32 / height as f32, camera.extent());
        let direction = self.inverse_cell * (inverse_rotation * Vec3::new(0.0, 0.0, -1.0));
        let shape = self.volume_shape;
        let values = &self.volume;
        let texel = |p: Vec3| -> f32 {
            // GL_LINEAR with clamp-to-edge at texture coordinates (p (n-1) + .5) / n.
            let mut i0 = [0usize; 3];
            let mut i1 = [0usize; 3];
            let mut t = [0f32; 3];
            for a in 0..3 {
                let n = shape[a];
                let x = (p[a] * (n - 1) as f32).clamp(0.0, (n - 1) as f32);
                let f = x.floor();
                i0[a] = f as usize;
                i1[a] = (i0[a] + 1).min(n - 1);
                t[a] = x - f;
            }
            let at = |x: usize, y: usize, z: usize| values[(z * shape[1] + y) * shape[0] + x];
            let mut v = 0.0;
            for (z, wz) in [(i0[2], 1.0 - t[2]), (i1[2], t[2])] {
                for (y, wy) in [(i0[1], 1.0 - t[1]), (i1[1], t[1])] {
                    for (x, wx) in [(i0[0], 1.0 - t[0]), (i1[0], t[0])] {
                        v += at(x, y, z) * wx * wy * wz;
                    }
                }
            }
            v
        };
        let background = Vec3::from_array(BACKGROUND);
        let pixel = |px: usize, py: usize| -> [u8; 4] {
            let sx = (px as f32 + 0.5) / width as f32 * 2.0 - 1.0;
            let sy = 1.0 - (py as f32 + 0.5) / height as f32 * 2.0;
            let origin = self.inverse_cell * (inverse_rotation * Vec3::new(sx * extent * aspect - camera.pan[0], sy * extent - camera.pan[1], 2.0)) + Vec3::splat(0.5);
            let (mut entry, mut leave) = (0.0f32, 1e20f32);
            for a in 0..3 {
                if direction[a].abs() < 1e-8 {
                    if origin[a] < 0.0 || origin[a] > 1.0 {
                        return rgba(background, 1.0);
                    }
                } else {
                    let (x, y) = (-origin[a] / direction[a], (1.0 - origin[a]) / direction[a]);
                    entry = entry.max(x.min(y));
                    leave = leave.min(x.max(y));
                }
            }
            if leave <= entry {
                return rgba(background, 1.0);
            }
            let step = (leave - entry) / steps as f32;
            let (mut colour, mut alpha) = (Vec3::ZERO, 0.0f32);
            for i in 0..steps {
                let p = origin + direction * (entry + (i as f32 + 0.5) * step);
                let value = texel(p);
                if value < threshold {
                    continue;
                }
                let t = if high > low { ((value - low) / (high - low)).clamp(0.0, 1.0) } else { 0.5 };
                let depth = 0.65 * opacity / (1.0 - opacity).max(1e-6);
                let a = if opacity >= 1.0 { 1.0 } else { 1.0 - (-depth * 12.0 * t.max(0.05) * step).exp() };
                colour += (1.0 - alpha) * a * ramp(t, palette);
                alpha += (1.0 - alpha) * a;
                if alpha > 0.995 {
                    break;
                }
            }
            rgba(colour + (1.0 - alpha) * background, 1.0)
        };
        let mut out = vec![0u8; width * height * 4];
        let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).min(16);
        let rows_per = height.div_ceil(threads).max(1);
        std::thread::scope(|scope| {
            for (chunk, rows) in out.chunks_mut(rows_per * width * 4).enumerate() {
                let pixel = &pixel;
                scope.spawn(move || {
                    for (i, px) in rows.chunks_exact_mut(4).enumerate() {
                        let (x, y) = (i % width, chunk * rows_per + i / width);
                        px.copy_from_slice(&pixel(x, y));
                    }
                });
            }
        });
        out
    }
}

fn rgba(c: Vec3, a: f32) -> [u8; 4] {
    let b = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
    [b(c.x), b(c.y), b(c.z), b(a)]
}

/// Software rasterisation of sorted, alpha-blended triangles (PNG export of a surface).
fn rasterize(width: usize, height: usize, triangles: &[[(egui::Vec2, [f32; 4]); 3]]) -> Vec<u8> {
    let mut image: Vec<[f32; 3]> = vec![BACKGROUND; width * height];
    for tri in triangles {
        let p: Vec<(f32, f32)> = tri.iter().map(|(v, _)| ((v.x + 1.0) * 0.5 * width as f32, (1.0 - v.y) * 0.5 * height as f32)).collect();
        let area = (p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1);
        if area.abs() < 1e-12 {
            continue;
        }
        let x0 = p.iter().map(|q| q.0).fold(f32::MAX, f32::min).floor().max(0.0) as usize;
        let x1 = (p.iter().map(|q| q.0).fold(f32::MIN, f32::max).ceil() as isize).min(width as isize - 1);
        let y0 = p.iter().map(|q| q.1).fold(f32::MAX, f32::min).floor().max(0.0) as usize;
        let y1 = (p.iter().map(|q| q.1).fold(f32::MIN, f32::max).ceil() as isize).min(height as isize - 1);
        if x1 < 0 || y1 < 0 {
            continue;
        }
        for y in y0..=y1 as usize {
            for x in x0..=x1 as usize {
                let (cx, cy) = (x as f32 + 0.5, y as f32 + 0.5);
                let w0 = ((p[1].0 - cx) * (p[2].1 - cy) - (p[2].0 - cx) * (p[1].1 - cy)) / area;
                let w1 = ((p[2].0 - cx) * (p[0].1 - cy) - (p[0].0 - cx) * (p[2].1 - cy)) / area;
                let w2 = 1.0 - w0 - w1;
                if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 {
                    continue;
                }
                let c: [f32; 4] = std::array::from_fn(|k| tri[0].1[k] * w0 + tri[1].1[k] * w1 + tri[2].1[k] * w2);
                let dst = &mut image[y * width + x];
                for k in 0..3 {
                    dst[k] = c[k] * c[3] + dst[k] * (1.0 - c[3]);
                }
            }
        }
    }
    image.into_iter().flat_map(|c| rgba(Vec3::from_array(c), 1.0)).collect()
}

fn write_png(path: &str, width: usize, height: usize, rgba: &[u8]) -> Result<(), String> {
    let file = std::fs::File::create(path).map_err(|_| "Failed to write viewport PNG".to_string())?;
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), width as u32, height as u32);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().map_err(|e| e.to_string())?;
    writer.write_image_data(rgba).map_err(|_| "Failed to write viewport PNG".to_string())
}

pub struct ElectronicDialog {
    open: bool,
    /// Files dropped on the window load as the reference volume.
    drop_reference: bool,
    output: String,
    reflections: String,
    charges: String,
    thresholds: String,
    quantity: i32,
    cube_units: i32,
    selected: i32,
    reference: i32,
    local_reference: bool,
    operation: i32,
    tool_group: i32,
    axis: i32,
    count: i32,
    window: i32,
    radius: i32,
    export_format: i32,
    scalar: f32,
    mask_low: f32,
    mask_high: f32,
    reference_weight: f32,
    boolean_operation: i32,
    surface_level: f32,
    suggested_level: f32,
    render_mode: i32,
    view_layout: i32,
    show_slice_plane: bool,
    sigma: f32,
    bonding_floor: f32,
    start: [f32; 3],
    end: [f32; 3],
    u: [f32; 3],
    v: [f32; 3],
    alpha: f32,
    real_cutoff: f32,
    reciprocal_cutoff: f32,
    yaw: f32,
    pitch: f32,
    opacity: f32,
    zoom: f32,
    pan: [f32; 2],
    specular: f32,
    shininess: f32,
    colour_low: f32,
    colour_high: f32,
    auto_low: f32,
    auto_high: f32,
    palette: i32,
    quality: i32,
    auto_range: bool,
    colour_unit: String,
    loaded_path: String,
    reference_path: String,
    slice_field: i32,
    slice_low: f32,
    slice_high: f32,
    colour_surface: bool,
    append_surface: bool,
    error: String,
    volume: Volume,
    reference_volume: Volume,
    result: Output,
    surface: Mesh,
    task: Option<Receiver<Result<Output, String>>>,
    pending_operation: String,
    undo: VecDeque<Workspace>,
    redo: VecDeque<Workspace>,
    history: Vec<String>,
    slice: SliceView,
    view: View3D,
    volume_dirty: bool,
    generate_surface: bool,
    picker: FileBrowser,
    picker_action: i32,
}

impl Default for ElectronicDialog {
    fn default() -> Self {
        Self {
            open: false,
            drop_reference: false,
            output: "electronic-result.xsf".into(),
            reflections: "0 0 0\n1 0 0\n-1 0 0".into(),
            charges: String::new(),
            thresholds: "0.1\n0.5\n1.0".into(),
            quantity: 0,
            cube_units: 0,
            selected: 0,
            reference: 0,
            local_reference: false,
            operation: 25,
            tool_group: 1,
            axis: 2,
            count: 100,
            window: 3,
            radius: 2,
            export_format: 0,
            scalar: 0.1,
            mask_low: 0.0,
            mask_high: 1.0,
            reference_weight: 1.0,
            boolean_operation: 0,
            surface_level: 0.1,
            suggested_level: 0.1,
            render_mode: 0,
            view_layout: 0,
            show_slice_plane: true,
            sigma: 0.5,
            bonding_floor: 1e-12,
            start: [0.0; 3],
            end: [1.0, 0.0, 0.0],
            u: [1.0, 0.0, 0.0],
            v: [0.0, 1.0, 0.0],
            alpha: 1.0,
            real_cutoff: 8.0,
            reciprocal_cutoff: 12.0,
            yaw: 0.6,
            pitch: 0.4,
            opacity: 0.35,
            zoom: 1.0,
            pan: [0.0; 2],
            specular: 0.28,
            shininess: 40.0,
            colour_low: 0.0,
            colour_high: 1.0,
            auto_low: 0.0,
            auto_high: 1.0,
            palette: 0,
            quality: 1,
            auto_range: true,
            colour_unit: "raw".into(),
            loaded_path: String::new(),
            reference_path: String::new(),
            slice_field: -1,
            slice_low: 0.0,
            slice_high: 1.0,
            colour_surface: false,
            append_surface: false,
            error: String::new(),
            volume: Volume::default(),
            reference_volume: Volume::default(),
            result: Output::default(),
            surface: Mesh::default(),
            task: None,
            pending_operation: String::new(),
            undo: VecDeque::new(),
            redo: VecDeque::new(),
            history: Vec::new(),
            slice: SliceView::default(),
            view: View3D { quality: 1, ..Default::default() },
            volume_dirty: true,
            generate_surface: false,
            picker: FileBrowser::default(),
            picker_action: 0,
        }
    }
}

fn v3(p: [f32; 3]) -> DVec3 {
    DVec3::new(p[0] as f64, p[1] as f64, p[2] as f64)
}

fn combo(ui: &mut egui::Ui, label: &str, selected: &mut i32, items: &[&str]) -> bool {
    let mut changed = false;
    field(ui, label, |ui| {
        egui::ComboBox::from_id_salt(label).width(ui.available_width()).selected_text(items.get(*selected as usize).copied().unwrap_or("")).show_ui(ui, |ui| {
            for (i, item) in items.iter().enumerate() {
                changed |= ui.selectable_value(selected, i as i32, *item).changed();
            }
        });
    });
    changed
}

fn input_float(ui: &mut egui::Ui, label: &str, value: &mut f32, decimals: Option<usize>) {
    field(ui, label, |ui| {
        let mut d = egui::DragValue::new(value).speed(0.01);
        d = match decimals {
            Some(n) => d.fixed_decimals(n),
            None => d.custom_formatter(|v, _| format_g(v, 6)),
        };
        ui.add_sized([ui.available_width(), 18.0], d);
    });
}

fn input_int(ui: &mut egui::Ui, label: &str, value: &mut i32) {
    field(ui, label, |ui| {
        ui.add_sized([ui.available_width(), 18.0], egui::DragValue::new(value).speed(0.2));
    });
}

fn input_vector(ui: &mut egui::Ui, label: &str, value: &mut [f32; 3]) {
    field(ui, label, |ui| {
        let w = (ui.available_width() - 8.0) / 3.0;
        ui.horizontal(|ui| {
            for c in value.iter_mut() {
                ui.add_sized([w, 18.0], egui::DragValue::new(c).speed(0.01).fixed_decimals(3));
            }
        });
    });
}

fn colour32(c: Vec3, a: f32) -> egui::Color32 {
    egui::Color32::from_rgba_unmultiplied((c.x.clamp(0.0, 1.0) * 255.0) as u8, (c.y.clamp(0.0, 1.0) * 255.0) as u8, (c.z.clamp(0.0, 1.0) * 255.0) as u8, (a.clamp(0.0, 1.0) * 255.0) as u8)
}

impl ElectronicDialog {
    /// The settings stored with a workspace (the C++ names, for project compatibility).
    fn settings(&self) -> BTreeMap<String, f64> {
        let f = |v: f32| v as f64;
        let i = |v: i32| v as f64;
        let b = |v: bool| v as i32 as f64;
        let entries: Vec<(&str, f64)> = vec![
            ("m_quantity", i(self.quantity)),
            ("m_cubeUnits", i(self.cube_units)),
            ("m_selected", i(self.selected)),
            ("m_reference", i(self.reference)),
            ("m_localReference", b(self.local_reference)),
            ("m_operation", i(self.operation)),
            ("m_toolGroup", i(self.tool_group)),
            ("m_axis", i(self.axis)),
            ("m_count", i(self.count)),
            ("m_window", i(self.window)),
            ("m_radius", i(self.radius)),
            ("m_exportFormat", i(self.export_format)),
            ("m_scalar", f(self.scalar)),
            ("m_maskLow", f(self.mask_low)),
            ("m_maskHigh", f(self.mask_high)),
            ("m_referenceWeight", f(self.reference_weight)),
            ("m_booleanOperation", i(self.boolean_operation)),
            ("m_surfaceLevel", f(self.surface_level)),
            ("m_suggestedLevel", f(self.suggested_level)),
            ("m_renderMode", i(self.render_mode)),
            ("m_viewLayout", i(self.view_layout)),
            ("m_showSlicePlane", b(self.show_slice_plane)),
            ("m_sigma", f(self.sigma)),
            ("m_bondingFloor", f(self.bonding_floor)),
            ("m_alpha", f(self.alpha)),
            ("m_realCutoff", f(self.real_cutoff)),
            ("m_reciprocalCutoff", f(self.reciprocal_cutoff)),
            ("m_yaw", f(self.yaw)),
            ("m_pitch", f(self.pitch)),
            ("m_opacity", f(self.opacity)),
            ("m_zoom", f(self.zoom)),
            ("m_specular", f(self.specular)),
            ("m_shininess", f(self.shininess)),
            ("m_colorLow", f(self.colour_low)),
            ("m_colorHigh", f(self.colour_high)),
            ("m_autoLow", f(self.auto_low)),
            ("m_autoHigh", f(self.auto_high)),
            ("m_palette", i(self.palette)),
            ("m_quality", i(self.quality)),
            ("m_autoRange", b(self.auto_range)),
            ("m_sliceLow", f(self.slice_low)),
            ("m_sliceHigh", f(self.slice_high)),
            ("m_colorSurface", b(self.colour_surface)),
            ("m_appendSurface", b(self.append_surface)),
            ("m_pan.x", f(self.pan[0])),
            ("m_pan.y", f(self.pan[1])),
            ("m_start[0]", f(self.start[0])),
            ("m_start[1]", f(self.start[1])),
            ("m_start[2]", f(self.start[2])),
            ("m_end[0]", f(self.end[0])),
            ("m_end[1]", f(self.end[1])),
            ("m_end[2]", f(self.end[2])),
            ("m_u[0]", f(self.u[0])),
            ("m_u[1]", f(self.u[1])),
            ("m_u[2]", f(self.u[2])),
            ("m_v[0]", f(self.v[0])),
            ("m_v[1]", f(self.v[1])),
            ("m_v[2]", f(self.v[2])),
        ];
        entries.into_iter().map(|(k, v)| (k.to_string(), v)).collect()
    }

    /// The dialog as a workspace (ElectronicPostProcessingDialog::snapshot).
    pub fn snapshot(&self, include_data: bool) -> Workspace {
        let mut out = Workspace::default();
        if include_data {
            out.volume = self.volume.clone();
            out.reference = self.reference_volume.clone();
            out.surface = self.surface.clone();
            out.table = self.result.table.clone();
            out.columns = self.result.columns;
            out.heading = self.result.heading.clone();
        }
        out.source_path = self.loaded_path.clone();
        out.reference_path = self.reference_path.clone();
        out.history = self.history.clone();
        out.slice_settings = self.slice.settings();
        out.settings = self.settings();
        out
    }

    /// Restores a workspace (ElectronicPostProcessingDialog::restore).
    pub fn restore(&mut self, saved: &Workspace) -> Result<(), String> {
        if self.task.is_some() {
            return Err("Wait for or cancel the calculation before restoring a workspace".into());
        }
        for settings in [&saved.settings, &saved.slice_settings] {
            for (name, &value) in settings {
                if !value.is_finite() || value.abs() > i32::MAX as f64 {
                    return Err(format!("Workspace setting outside the supported range: {name}"));
                }
            }
        }
        self.volume = saved.volume.clone();
        self.reference_volume = saved.reference.clone();
        self.surface = saved.surface.clone();
        self.result = Output { table: saved.table.clone(), columns: saved.columns, heading: saved.heading.clone(), ..Default::default() };
        self.loaded_path = saved.source_path.clone();
        self.reference_path = saved.reference_path.clone();
        self.history = saved.history.clone();
        let s = &saved.settings;
        let int = |name: &str, target: &mut i32| {
            if let Some(&v) = s.get(name) {
                *target = v as i32;
            }
        };
        let float = |name: &str, target: &mut f32| {
            if let Some(&v) = s.get(name) {
                *target = v as f32;
            }
        };
        let boolean = |name: &str, target: &mut bool| {
            if let Some(&v) = s.get(name) {
                *target = v != 0.0;
            }
        };
        int("m_quantity", &mut self.quantity);
        int("m_cubeUnits", &mut self.cube_units);
        int("m_selected", &mut self.selected);
        int("m_reference", &mut self.reference);
        boolean("m_localReference", &mut self.local_reference);
        int("m_operation", &mut self.operation);
        int("m_toolGroup", &mut self.tool_group);
        int("m_axis", &mut self.axis);
        int("m_count", &mut self.count);
        int("m_window", &mut self.window);
        int("m_radius", &mut self.radius);
        int("m_exportFormat", &mut self.export_format);
        float("m_scalar", &mut self.scalar);
        float("m_maskLow", &mut self.mask_low);
        float("m_maskHigh", &mut self.mask_high);
        float("m_referenceWeight", &mut self.reference_weight);
        int("m_booleanOperation", &mut self.boolean_operation);
        float("m_surfaceLevel", &mut self.surface_level);
        float("m_suggestedLevel", &mut self.suggested_level);
        int("m_renderMode", &mut self.render_mode);
        int("m_viewLayout", &mut self.view_layout);
        boolean("m_showSlicePlane", &mut self.show_slice_plane);
        float("m_sigma", &mut self.sigma);
        float("m_bondingFloor", &mut self.bonding_floor);
        float("m_alpha", &mut self.alpha);
        float("m_realCutoff", &mut self.real_cutoff);
        float("m_reciprocalCutoff", &mut self.reciprocal_cutoff);
        float("m_yaw", &mut self.yaw);
        float("m_pitch", &mut self.pitch);
        float("m_opacity", &mut self.opacity);
        float("m_zoom", &mut self.zoom);
        float("m_specular", &mut self.specular);
        float("m_shininess", &mut self.shininess);
        float("m_colorLow", &mut self.colour_low);
        float("m_colorHigh", &mut self.colour_high);
        float("m_autoLow", &mut self.auto_low);
        float("m_autoHigh", &mut self.auto_high);
        int("m_palette", &mut self.palette);
        int("m_quality", &mut self.quality);
        boolean("m_autoRange", &mut self.auto_range);
        float("m_sliceLow", &mut self.slice_low);
        float("m_sliceHigh", &mut self.slice_high);
        boolean("m_colorSurface", &mut self.colour_surface);
        boolean("m_appendSurface", &mut self.append_surface);
        float("m_pan.x", &mut self.pan[0]);
        float("m_pan.y", &mut self.pan[1]);
        for (k, target) in ["m_start[0]", "m_start[1]", "m_start[2]"].iter().zip(self.start.iter_mut()) {
            float(k, target);
        }
        for (k, target) in ["m_end[0]", "m_end[1]", "m_end[2]"].iter().zip(self.end.iter_mut()) {
            float(k, target);
        }
        for (k, target) in ["m_u[0]", "m_u[1]", "m_u[2]"].iter().zip(self.u.iter_mut()) {
            float(k, target);
        }
        for (k, target) in ["m_v[0]", "m_v[1]", "m_v[2]"].iter().zip(self.v.iter_mut()) {
            float(k, target);
        }
        self.selected = self.selected.clamp(0, 0.max(self.volume.fields.len() as i32 - 1));
        let references = if self.local_reference || self.reference_volume.fields.is_empty() { &self.volume } else { &self.reference_volume };
        self.reference = self.reference.clamp(0, 0.max(references.fields.len() as i32 - 1));
        self.operation = self.operation.clamp(0, 37);
        self.axis = self.axis.clamp(0, 2);
        self.quantity = self.quantity.clamp(0, 4);
        self.cube_units = self.cube_units.clamp(0, 1);
        self.export_format = self.export_format.clamp(0, 7);
        self.boolean_operation = self.boolean_operation.clamp(0, 3);
        self.quality = self.quality.clamp(0, 3);
        self.palette = self.palette.clamp(0, 2);
        self.render_mode = self.render_mode.clamp(0, 1);
        self.view_layout = self.view_layout.clamp(0, 2);
        self.opacity = self.opacity.clamp(0.0, 1.0);
        self.zoom = self.zoom.clamp(0.1, 100.0);
        self.view.quality = self.quality;
        let surface = self.surface.clone();
        self.view.set_mesh(&surface);
        self.volume_dirty = true;
        self.slice.set_resolution(SLICE_RESOLUTIONS[self.quality as usize]);
        self.slice.restore(&saved.slice_settings);
        self.slice_field = self.selected;
        self.open = !self.volume.fields.is_empty() || !self.reference_volume.fields.is_empty();
        Ok(())
    }

    /// Keeps the state before a calculation (four states, 32 MiB each).
    fn remember(&mut self) {
        let mut bytes = self.surface.vertices.len() * 24 + self.surface.colors.len() * 8;
        for g in self.volume.fields.iter().chain(&self.reference_volume.fields) {
            bytes += g.values.len() * 8;
        }
        self.redo.clear();
        if bytes > 32 * 1024 * 1024 {
            self.undo.clear();
            return;
        }
        self.undo.push_back(self.snapshot(true));
        while self.undo.len() > 4 {
            self.undo.pop_front();
        }
    }

    fn undo_calculation(&mut self) -> Result<(), String> {
        let previous = self.undo.pop_back().ok_or("Nothing to undo")?;
        self.redo.push_back(self.snapshot(true));
        self.restore(&previous)
    }

    fn redo_calculation(&mut self) -> Result<(), String> {
        let next = self.redo.pop_back().ok_or("Nothing to redo")?;
        self.undo.push_back(self.snapshot(true));
        self.restore(&next)
    }

    fn reset_camera(&mut self) {
        self.yaw = 0.6;
        self.pitch = 0.4;
        self.zoom = 1.0;
        self.pan = [0.0; 2];
    }

    fn camera(&self) -> Camera {
        Camera { yaw: self.yaw, pitch: self.pitch, zoom: self.zoom, pan: self.pan }
    }

    fn clear(&mut self) {
        self.volume = Volume::default();
        self.reference_volume = Volume::default();
        self.result = Output::default();
        self.surface = Mesh::default();
        self.view.set_mesh(&Mesh::default());
        self.selected = 0;
        self.reference = 0;
        self.slice_field = -1;
        self.loaded_path.clear();
        self.reference_path.clear();
        self.reset_camera();
    }

    fn reference_volume(&self) -> &Volume {
        if self.reference_volume.fields.is_empty() || self.local_reference {
            &self.volume
        } else {
            &self.reference_volume
        }
    }

    /// A loaded volume, a result or an error (from the background or a script).
    fn apply(&mut self, result: Result<Output, String>) {
        let out = match result {
            Ok(o) => o,
            Err(e) => {
                self.error = e;
                return;
            }
        };
        self.error.clear();
        self.remember();
        if !self.pending_operation.is_empty() {
            self.history.push(std::mem::take(&mut self.pending_operation));
        }
        self.result = out;
        if self.result.surface_ready {
            if !self.result.append_mesh {
                self.surface = Mesh::default();
            }
            if self.surface.vertices.len() + self.result.mesh.vertices.len() <= MESH_LIMIT {
                self.surface.vertices.extend(self.result.mesh.vertices.iter().copied());
                self.surface.colors.extend(self.result.mesh.colors.iter().copied());
                let surface = self.surface.clone();
                self.view.set_mesh(&surface);
                let (low, high) = (self.result.colour_low as f32, self.result.colour_high as f32);
                self.auto_low = if self.result.append_mesh { self.auto_low.min(low) } else { low };
                self.auto_high = if self.result.append_mesh { self.auto_high.max(high) } else { high };
                self.colour_unit = self.result.colour_unit.clone();
                if self.auto_range {
                    self.colour_low = self.auto_low;
                    self.colour_high = self.auto_high;
                }
            } else {
                self.error = "Combined surfaces exceed the mesh limit; clear or replace the surface.".into();
            }
            self.result.mesh = Mesh::default();
        }
        if self.result.reference_loaded {
            self.reference_volume = std::mem::take(&mut self.result.volume);
            self.reference = 0;
            self.reference_path = self.result.source_path.clone();
        } else if self.result.loaded {
            self.volume = std::mem::take(&mut self.result.volume);
            self.selected = 0;
            self.reference = 0;
            self.loaded_path = self.result.source_path.clone();
            self.slice_field = -1;
            self.slice.invalidate();
            self.slice.reset();
            self.surface = Mesh::default();
            self.view.set_mesh(&Mesh::default());
            self.reset_camera();
        } else if !self.result.volume.fields.is_empty() {
            self.selected = self.volume.fields.len() as i32;
            let fields = std::mem::take(&mut self.result.volume.fields);
            self.volume.fields.extend(fields);
        }
    }

    fn start(&mut self, f: impl FnOnce() -> Result<Output, String> + Send + 'static, sync: bool) {
        self.error.clear();
        if sync {
            let r = f();
            self.apply(r);
            return;
        }
        let (tx, rx) = channel();
        std::thread::spawn(move || {
            let _ = tx.send(f());
        });
        self.task = Some(rx);
    }

    pub fn load(&mut self, path: &str, reference: bool, sync: bool) {
        let q = QUANTITY_NAMES[self.quantity.clamp(0, 4) as usize].to_string();
        let units = if self.cube_units == 1 { "angstrom" } else { "bohr" }.to_string();
        self.pending_operation = format!("{}{path}", if reference { "Load reference: " } else { "Load: " });
        let path = path.to_string();
        self.start(
            move || {
                let volume = load_volume(&path, &q, &units)?;
                Ok(Output { volume, loaded: !reference, reference_loaded: reference, source_path: path, ..Default::default() })
            },
            sync,
        );
    }

    fn job(&self, op: i32) -> Result<Job, String> {
        let source = self.volume.fields.get(self.selected as usize).cloned().ok_or("Load a volume first.")?;
        let reference = self.reference_volume().fields.get(self.reference as usize).cloned().ok_or("Load a volume first.")?;
        Ok(Job {
            op,
            source,
            reference,
            sites: self.volume.sites.clone(),
            count: self.count,
            axis: self.axis,
            window: self.window,
            radius: self.radius,
            mask_low: self.mask_low as f64,
            mask_high: self.mask_high as f64,
            reference_weight: self.reference_weight as f64,
            boolean_operation: self.boolean_operation.clamp(0, 3) as usize,
            scalar: if op == 22 { self.surface_level as f64 } else { self.scalar as f64 },
            sigma: self.sigma as f64,
            alpha: self.alpha as f64,
            real_cutoff: self.real_cutoff as f64,
            reciprocal_cutoff: self.reciprocal_cutoff as f64,
            start: v3(self.start),
            end: v3(self.end),
            u: v3(self.u),
            v: v3(self.v),
            reflections: self.reflections.clone(),
            charges: self.charges.clone(),
            thresholds: self.thresholds.clone(),
            bonding_floor: self.bonding_floor as f64,
            coloured: self.colour_surface,
            append_surface: self.append_surface,
            existing_surface: !self.surface.vertices.is_empty(),
            existing_unit: self.colour_unit.clone(),
        })
    }

    /// Calculate (or Update Surface when `op` is 22 from the 3D view).
    pub fn calculate(&mut self, op: i32, sync: bool) {
        let job = match self.job(op) {
            Ok(j) => j,
            Err(e) => {
                self.error = e;
                return;
            }
        };
        let mut provenance = format!("Operation {op}, input {}, reference {}", job.source.name, job.reference.name);
        for (k, v) in self.settings() {
            provenance += &format!("; {k}={}", format_g(v, 6));
        }
        self.pending_operation = provenance;
        self.start(move || run(job), sync);
    }

    /// Export by the selected format (ElectronicPostProcessingDialog::save).
    pub fn save(&mut self, path: &str) {
        let r = self.save_inner(path);
        self.error = match r {
            Ok(()) => {
                self.output = path.to_string();
                format!("Saved {}", file_name(path))
            }
            Err(e) => e,
        };
    }

    fn save_inner(&mut self, path: &str) -> Result<(), String> {
        if self.volume.fields.is_empty() {
            return Err("Load a volume first.".into());
        }
        let grid = self.volume.fields[self.selected as usize].clone();
        match self.export_format {
            0..=2 => save_volume(&Volume { sites: self.volume.sites.clone(), fields: vec![grid] }, path, ["xsf", "cube", "vasp"][self.export_format as usize]),
            6 => {
                let (width, height) = (1600usize, 1200usize);
                let pixels = if self.render_mode == 0 {
                    self.view.quality = self.quality;
                    self.view.set_volume(&grid)?;
                    self.volume_dirty = true;
                    self.view.render_volume(width, height, self.camera(), self.opacity, self.surface_level, self.slice_low, self.slice_high, self.palette)
                } else {
                    if self.surface.vertices.is_empty() {
                        return Err("Calculate an isosurface first".into());
                    }
                    let tris = self.view.surface_triangles(width as f32 / height as f32, self.camera(), self.opacity, self.colour_low, self.colour_high, self.palette, self.specular, self.shininess);
                    rasterize(width, height, &tris)
                };
                write_png(path, width, height, &pixels)
            }
            7 => {
                let svg = self.slice.svg(&grid, self.slice_low, self.slice_high, self.palette)?;
                std::fs::write(path, svg).map_err(|_| "Failed to write section SVG".to_string())
            }
            4 | 5 => {
                if self.surface.vertices.is_empty() {
                    return Err("Calculate an isosurface first.".into());
                }
                save_mesh(&self.surface, std::path::Path::new(path))
            }
            _ => {
                if self.result.columns == 0 {
                    return Err("No table to export".into());
                }
                let mut text = format!("# {}\n", self.result.heading);
                let columns = self.result.columns as usize;
                for (i, v) in self.result.table.iter().enumerate() {
                    text += &format_g(*v, 17);
                    text.push(if (i + 1) % columns != 0 { ',' } else { '\n' });
                }
                std::fs::write(path, text).map_err(|_| "Cannot open CSV output".to_string())
            }
        }
    }

    pub fn save_workspace(&self, path: &str) -> Result<(), String> {
        save_project(&[self.snapshot(true)], path, FORMAT)
    }

    pub fn open_workspace(&mut self, path: &str) -> Result<(), String> {
        let saved = load_project(path)?;
        let first = saved.into_iter().next().ok_or("The workspace has no tabs")?;
        self.remember();
        self.restore(&first)
    }

    /// The field's display range, when the selected field changed (drawPreview).
    fn refresh_field(&mut self) {
        if self.volume.fields.is_empty() || self.slice_field == self.selected {
            return;
        }
        let grid = &self.volume.fields[self.selected as usize];
        if let Ok(r) = display_range(&grid.values) {
            self.slice_low = r.low as f32;
            self.slice_high = r.high as f32;
            self.surface_level = r.suggested as f32;
            self.suggested_level = r.suggested as f32;
        }
        let unit = grid.unit.clone();
        self.volume_dirty = true;
        self.slice_field = self.selected;
        self.slice.invalidate();
        if self.surface.vertices.is_empty() {
            self.auto_low = self.slice_low;
            self.auto_high = self.slice_high;
            self.colour_unit = unit;
            if self.auto_range {
                self.colour_low = self.auto_low;
                self.colour_high = self.auto_high;
            }
        }
    }

    fn set_export_format(&mut self, format: i32) {
        self.export_format = format;
        let mut p = std::path::PathBuf::from(&self.output);
        p.set_extension(&EXPORT_EXTENSIONS[format as usize][1..]);
        self.output = p.to_string_lossy().into_owned();
    }

    fn select_field(ui: &mut egui::Ui, label: &str, selected: &mut i32, volume: &Volume) {
        *selected = (*selected).clamp(0, volume.fields.len() as i32 - 1);
        field(ui, label, |ui| {
            egui::ComboBox::from_id_salt(label).width(ui.available_width()).selected_text(volume.fields[*selected as usize].name.clone()).show_ui(ui, |ui| {
                for (i, f) in volume.fields.iter().enumerate() {
                    ui.selectable_value(selected, i as i32, &f.name);
                }
            });
        });
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        let running = self.task.is_some();
        ui.add_enabled_ui(!running, |ui| {
            ui.horizontal(|ui| {
                if ui.button("Open Workspace").clicked() {
                    self.picker_action = 3;
                    self.picker.show_for(Purpose::OpenProject, "Open electronic workspace", "workspace.afproject", &["afproject"]);
                }
                if ui.button("Save Workspace").clicked() {
                    self.picker_action = 4;
                    self.picker.show_for(Purpose::SaveProject, "Save electronic workspace", "workspace.afproject", &["afproject"]);
                }
            });
            ui.horizontal(|ui| {
                if ui.add_enabled(!self.undo.is_empty(), egui::Button::new("Undo Calculation")).clicked() {
                    if let Err(e) = self.undo_calculation() {
                        self.error = e;
                    }
                }
                if ui.add_enabled(!self.redo.is_empty(), egui::Button::new("Redo")).clicked() {
                    if let Err(e) = self.redo_calculation() {
                        self.error = e;
                    }
                }
            });
            ui.label("Undo retains 4 states up to 32 MiB each. Save larger workspaces to disk.");
            crate::dialog_widgets::CollapsingHeader::new("Operation history").show(ui, |ui| {
                for entry in &self.history {
                    ui.label(entry);
                }
            });
            heading(ui, "Tools");
            ui.horizontal(|ui| {
                if crate::dialog_widgets::radio(ui, self.tool_group == 1, "Charge transfer").clicked() {
                    self.tool_group = 1;
                    self.operation = 25;
                }
                if crate::dialog_widgets::radio(ui, self.tool_group == 0, "General analysis").clicked() {
                    self.tool_group = 0;
                    self.operation = 0;
                }
                if crate::dialog_widgets::radio(ui, self.tool_group == 2, "Bonding & topology").clicked() {
                    self.tool_group = 2;
                    self.operation = 33;
                }
            });
            let grid_tools = |ui: &mut egui::Ui, labels: &[&str], first: i32, operation: &mut i32| {
                // Two columns of ImGui selectables (full-cell highlight).
                let column = (ui.available_width() - ui.spacing().item_spacing.x) * 0.5;
                for pair in labels.chunks(2).enumerate() {
                    let (row, items) = pair;
                    ui.horizontal(|ui| {
                        ui.spacing_mut().item_spacing.y = 0.0;
                        for (k, label) in items.iter().enumerate() {
                            let index = first + (row * 2 + k) as i32;
                            if crate::dialog_widgets::selectable(ui, *operation == index, label, column).clicked() {
                                *operation = index;
                            }
                        }
                    });
                }
            };
            match self.tool_group {
                1 => grid_tools(ui, &CHARGE_TOOLS, 25, &mut self.operation),
                2 => grid_tools(ui, &BONDING_TOOLS, 33, &mut self.operation),
                _ => {
                    combo(ui, "Tool", &mut self.operation, &GENERAL_TOOLS);
                }
            }
            if self.tool_group != 0 && self.volume.fields.is_empty() {
                ui.label("Load a density file below to use the selected tool.");
            }
            heading(ui, "Input data");
            let half = (ui.available_width() - ui.spacing().item_spacing.x) * 0.5;
            ui.horizontal(|ui| {
                if crate::dialog_widgets::primary_button(ui, "Open Volume", egui::vec2(half, crate::theme::FRAME_HEIGHT)).clicked() {
                    self.picker_action = 0;
                    let name = file_name(&self.loaded_path);
                    self.picker.show_for(Purpose::Open, "Open electronic volume", &name, &[]);
                }
                if ui.add_sized([half, crate::theme::FRAME_HEIGHT], egui::Button::new("Reference")).clicked() {
                    self.picker_action = 1;
                    let name = file_name(if self.reference_path.is_empty() { &self.loaded_path } else { &self.reference_path });
                    self.picker.show_for(Purpose::Open, "Open reference volume", &name, &[]);
                }
            });
            ui.horizontal(|ui| {
                if ui.add_enabled(!self.loaded_path.is_empty(), egui::Button::new("Reload Volume").min_size(egui::vec2(half, crate::theme::FRAME_HEIGHT))).clicked() {
                    let path = self.loaded_path.clone();
                    self.load(&path, false, false);
                }
                if ui.add_sized([half, crate::theme::FRAME_HEIGHT], egui::Button::new("Clear")).clicked() {
                    self.clear();
                }
            });
            if !self.loaded_path.is_empty() {
                ui.label(format!("Volume: {}", file_name(&self.loaded_path))).on_hover_text(&self.loaded_path);
            }
            if !self.reference_path.is_empty() {
                ui.label(format!("Reference: {}", file_name(&self.reference_path))).on_hover_text(&self.reference_path);
            }
            crate::dialog_widgets::checkbox(ui, &mut self.drop_reference, "Drop files as reference");
            ui.add(egui::Label::new("Drop a VASP, Cube or XSF file onto the app to load it here.").wrap());
            crate::dialog_widgets::CollapsingHeader::new("Import options").show(ui, |ui| {
                combo(ui, "Quantity", &mut self.quantity, &QUANTITIES);
                combo(ui, "Cube units", &mut self.cube_units, &["Bohr", "Angstrom"]);
                ui.label("Choose import settings before browsing, or reload the selected volume to apply changes.");
            });
            if !self.volume.fields.is_empty() {
                self.analysis_controls(ui);
            }
        });
        if running {
            ui.label("Calculating...");
            if ui.button("Cancel Calculation").clicked() {
                // The calculation has no cancellation point here: its result is discarded.
                self.task = None;
                self.pending_operation.clear();
                self.error = "Calculation cancelled.".into();
            }
        }
        if !self.error.is_empty() {
            ui.label(&self.error);
        }
        crate::dialog_widgets::CollapsingHeader::new("Appearance").show(ui, |ui| {
            if combo(ui, "Preview quality", &mut self.quality, &QUALITIES) {
                self.view.quality = self.quality;
                self.slice.set_resolution(SLICE_RESOLUTIONS[self.quality as usize]);
                self.volume_dirty = true;
            }
            if combo(ui, "Colors", &mut self.palette, &PALETTES) {
                self.view.texture = None;
            }
            if crate::dialog_widgets::checkbox(ui, &mut self.auto_range, "Automatic range").changed() && self.auto_range {
                self.colour_low = self.auto_low;
                self.colour_high = self.auto_high;
            }
            if !self.auto_range {
                input_float(ui, "Minimum", &mut self.colour_low, None);
                input_float(ui, "Maximum", &mut self.colour_high, None);
            }
            crate::dialog_widgets::CollapsingHeader::new("Lighting").show(ui, |ui| {
                field(ui, "Specular", |ui| ui.add(crate::dialog_widgets::Slider::new(&mut self.specular, 0.0..=1.0).fixed_decimals(2)));
                field(ui, "Shininess", |ui| ui.add(crate::dialog_widgets::Slider::new(&mut self.shininess, 4.0..=128.0).fixed_decimals(0)));
            });
        });
    }

    fn analysis_controls(&mut self, ui: &mut egui::Ui) {
        heading(ui, "Analysis");
        let mut selected = self.selected;
        Self::select_field(ui, "Field", &mut selected, &self.volume);
        self.selected = selected;
        if !self.reference_volume.fields.is_empty() {
            crate::dialog_widgets::checkbox(ui, &mut self.local_reference, "Use loaded/result field as reference");
        }
        let reference_count = self.reference_volume().fields.len() as i32;
        self.reference = self.reference.clamp(0, reference_count - 1);
        let g = &self.volume.fields[self.selected as usize];
        ui.weak(format!("{} × {} × {} | {}", g.shape[0], g.shape[1], g.shape[2], unit_display(&g.unit)));
        let op = self.operation;
        let select_reference = |s: &mut Self, ui: &mut egui::Ui, label: &str| {
            let mut r = s.reference;
            let volume = if s.reference_volume.fields.is_empty() || s.local_reference { &s.volume } else { &s.reference_volume };
            Self::select_field(ui, label, &mut r, volume);
            s.reference = r;
        };
        if (1..=4).contains(&op) || op == 23 || op == 22 {
            select_reference(self, ui, "Reference");
        }
        if op == 5 {
            input_float(ui, "Scale factor", &mut self.scalar, Some(3));
        }
        if op == 6 {
            input_float(ui, "Sigma (A)", &mut self.sigma, Some(3));
            input_int(ui, "Radius (steps)", &mut self.radius);
        }
        if op == 9 {
            ui.label("Gradient-expansion kinetic and local-virial energy densities (eV/Å³). Requires nonnegative e/Å³. Samples <= 1e-12 e/Å³ are masked to zero.");
        }
        if op == 10 || op == 13 || op == 17 {
            input_vector(ui, "Origin (A)", &mut self.start);
        }
        if op == 10 {
            input_vector(ui, "End (A)", &mut self.end);
        }
        if op == 10 || op == 13 {
            input_int(ui, "Sample count", &mut self.count);
        }
        if op == 11 || op == 12 {
            combo(ui, "Plane normal", &mut self.axis, &AXES);
        }
        if op == 12 {
            input_int(ui, "Window (odd)", &mut self.window);
        }
        if op == 13 {
            input_vector(ui, "Plane span u (A)", &mut self.u);
            input_vector(ui, "Plane span v (A)", &mut self.v);
        }
        if op == 14 {
            input_float(ui, "Isovalue", &mut self.scalar, Some(3));
        }
        if op == 17 {
            input_float(ui, "Sphere radius (A)", &mut self.scalar, Some(3));
        }
        if op == 15 {
            ui.label("26-neighbor grid maxima; adjacent equal values use a deterministic tie break. Refine the grid to check peak positions.");
        }
        if op == 16 {
            ui.label("Integrates nearest-site Voronoi cells around imported atoms. Boundary samples share weight. These are geometric partitions, not Bader basins.");
        }
        if op == 18 || op == 19 {
            ui.label(if op == 18 { "One h k l per line" } else { "One h k l real imag per line; include conjugate pairs" });
            ui.add(egui::TextEdit::multiline(&mut self.reflections).desired_rows(5).desired_width(f32::INFINITY).code_editor());
        }
        if op == 21 {
            ui.label("Explicit charges in electron units, one per imported atom in file order. Requires a neutral cell; conducting boundary. Increase both cutoffs to check convergence.");
            ui.label("Charges");
            ui.add(egui::TextEdit::multiline(&mut self.charges).desired_rows(4).desired_width(f32::INFINITY).code_editor());
            input_float(ui, "Ewald alpha (1/A)", &mut self.alpha, Some(3));
            input_float(ui, "Real cutoff (A)", &mut self.real_cutoff, Some(3));
            input_float(ui, "Reciprocal cutoff (1/A, includes 2*pi)", &mut self.reciprocal_cutoff, Some(3));
        }
        if op == 22 {
            crate::dialog_widgets::checkbox(ui, &mut self.colour_surface, "Color from reference");
            crate::dialog_widgets::checkbox(ui, &mut self.append_surface, "Append surface level");
        }
        if op == 24 {
            ui.label("For Cube/XSF covering a complete periodic cell: checks all opposite endpoint values, then removes duplicate planes. Matching values alone do not establish physical periodicity.");
        }
        if op == 25 || op == 27 || op == 28 {
            select_reference(self, ui, "Operand");
        }
        if op == 25 {
            input_float(ui, "Reference weight", &mut self.reference_weight, Some(3));
            ui.label("Result = field - weight * reference. Use aligned fragment densities from the same cell and geometry. Select the result and repeat for additional fragments.");
        }
        if op == 26 {
            input_float(ui, "Lower bound", &mut self.mask_low, None);
            input_float(ui, "Upper bound", &mut self.mask_high, None);
            ui.label("Creates a binary mask for the inclusive value range.");
        }
        if op == 27 {
            combo(ui, "Boolean operation", &mut self.boolean_operation, &BOOLEAN_OPERATIONS);
        }
        if op == 28 {
            ui.label("The operand must be a binary mask. Keeps field values inside the mask; integrate the result to measure regional charge.");
        }
        if op == 29 || op == 30 {
            ui.label("Apply to a difference density in e/Å³. Accumulation and depletion are positive electron counts, not atom-assigned transfer charges.");
        }
        if op == 31 {
            combo(ui, "Plane normal", &mut self.axis, &AXES);
            ui.label("Cumulative electrons from the lower cell face. The final row includes the whole cell; periodic profiles depend on the cell origin.");
        }
        if (33..=35).contains(&op) {
            input_float(ui, "Density floor (e/Å³)", &mut self.bonding_floor, None);
            ui.label("Requires nonnegative electron density in e/Å³. Samples at or below the floor are masked to zero, avoiding the divergent gradient/density ratio near vacuum.");
            if op == 34 {
                ui.label("Sign(lambda_2)*rho: negative values mark bonding/hydrogen-bonding accumulation, positive values mark steric repulsion.");
            }
            if op == 35 {
                ui.label("Density Overlap Regions Indicator, in [0,1); values approaching 1 mark boundaries between density basins.");
            }
        }
        if op == 36 {
            ui.label("One density threshold per line (e/Å³). Finite (non-periodic) grids only.");
            ui.add(egui::TextEdit::multiline(&mut self.thresholds).desired_rows(4).desired_width(f32::INFINITY).code_editor());
        }
        if op == 37 {
            ui.label("Native QTAIM/Bader on-grid steepest-ascent partitioning (requires nonnegative e/Å³). Adds a basin-index field and a per-basin charge/volume/maximum table.");
        }
        ui.add_space(4.0);
        if ui.add(egui::Button::new(egui::RichText::new("Calculate").strong()).min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
            self.calculate(self.operation, false);
        }
        ui.add_space(4.0);
        crate::dialog_widgets::CollapsingHeader::new("Export").default_open(true).show(ui, |ui| {
            let mut format = self.export_format;
            if combo(ui, "Format", &mut format, &EXPORT_FORMATS) {
                self.set_export_format(format);
            }
            if ui.add(egui::Button::new("Save As").min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
                let mut output = std::path::PathBuf::from(&self.output);
                if output.parent().is_none_or(|p| p.as_os_str().is_empty()) && !self.loaded_path.is_empty() {
                    if let Some(parent) = std::path::Path::new(&self.loaded_path).parent() {
                        self.picker.folder = parent.to_path_buf();
                    }
                    output = output.file_name().map(std::path::PathBuf::from).unwrap_or(output);
                }
                self.picker_action = 2;
                self.picker.show_for(Purpose::SaveAs, "Save electronic result", &output.to_string_lossy(), &[]);
            }
        });
    }

    fn results(&mut self, ui: &mut egui::Ui) {
        if !self.result.table.is_empty() {
            crate::dialog_widgets::CollapsingHeader::new("Results table").default_open(true).show(ui, |ui| {
                ui.label(&self.result.heading);
                let t = &self.result.table;
                if self.result.heading == "distance_A,value" || self.result.heading == "distance_A,cumulative_e" {
                    let x: Vec<f64> = t.iter().step_by(2).copied().collect();
                    let y: Vec<f64> = t.iter().skip(1).step_by(2).copied().collect();
                    line_plot(ui, "Profile", "distance (Å)", "value", &x, vec![("Profile".into(), y)], 150.0);
                } else if self.result.heading == "threshold,betti0,betti1,betti2" {
                    let x: Vec<f64> = t.iter().step_by(4).copied().collect();
                    for (k, name) in ["Connected components (betti0)", "Loops (betti1)", "Enclosed cavities (betti2)"].iter().enumerate() {
                        let y: Vec<f64> = t.iter().skip(k + 1).step_by(4).copied().collect();
                        line_plot(ui, name, "threshold", "count", &x, vec![(name.to_string(), y)], 100.0);
                    }
                }
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    egui::ScrollArea::vertical().id_salt("electronic table").max_height(100.0).show(ui, |ui| {
                        ui.set_width(ui.available_width());
                        let columns = self.result.columns.max(1) as usize;
                        for row in t.chunks(columns).take(200) {
                            let text: String = row.iter().map(|v| format!("{}   ", format_g(*v, 9))).collect();
                            ui.label(egui::RichText::new(text).monospace());
                        }
                    });
                });
            });
        }
        self.preview(ui);
    }

    fn preview(&mut self, ui: &mut egui::Ui) {
        egui::ComboBox::from_id_salt("view layout").width(ui.available_width()).selected_text(["3D and 2D", "Only 3D", "Only 2D"][self.view_layout as usize]).show_ui(ui, |ui| {
            for (i, name) in ["3D and 2D", "Only 3D", "Only 2D"].iter().enumerate() {
                ui.selectable_value(&mut self.view_layout, i as i32, *name);
            }
        });
        if self.volume.fields.is_empty() {
            ui.label("Open a VASP, Cube or XSF file to display both views.");
            return;
        }
        self.refresh_field();
        let width = ui.available_width();
        let view_width = if self.view_layout == 0 { (width - ui.spacing().item_spacing.x) * 0.5 } else { width };
        let height = ui.available_height().max(420.0);
        let top_down = egui::Layout::top_down(egui::Align::Min);
        ui.horizontal_top(|ui| {
            if self.view_layout != 2 {
                ui.allocate_ui_with_layout(egui::vec2(view_width, height), top_down, |ui| {
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        ui.set_width(view_width - 14.0);
                        self.view_3d(ui, height - 14.0);
                    });
                });
            }
            if self.view_layout != 1 {
                ui.allocate_ui_with_layout(egui::vec2(view_width, height), top_down, |ui| {
                    egui::Frame::group(ui.style()).show(ui, |ui| {
                        ui.set_width(view_width - 14.0);
                        self.view_2d(ui, height - 14.0);
                    });
                });
            }
        });
    }

    fn view_3d(&mut self, ui: &mut egui::Ui, height: f32) {
        let top = ui.cursor().top();
        let grid_unit = self.volume.fields[self.selected as usize].unit.clone();
        ui.label("3D view");
        egui::ComboBox::from_id_salt("render mode").width(ui.available_width()).selected_text(["Volume (interior density)", "Isosurface (boundary)"][self.render_mode as usize]).show_ui(ui, |ui| {
            ui.selectable_value(&mut self.render_mode, 0, "Volume (interior density)");
            ui.selectable_value(&mut self.render_mode, 1, "Isosurface (boundary)");
        });
        let running = self.task.is_some();
        ui.add_enabled_ui(!running, |ui| {
            ui.label(if self.render_mode == 0 { "Density threshold" } else { "3D isovalue" });
            ui.add_sized([ui.available_width(), 18.0], egui::DragValue::new(&mut self.surface_level).speed(0.001).custom_formatter(|v, _| format_g(v, 5)));
            if ui.button("Estimate Level").clicked() {
                self.surface_level = self.suggested_level;
            }
            ui.label(format!("Field range: {} to {} {}", format_g(self.slice_low as f64, 5), format_g(self.slice_high as f64, 5), unit_display(&grid_unit)));
            let invalid = !self.surface_level.is_finite() || self.surface_level < self.slice_low || self.surface_level > self.slice_high;
            if invalid {
                ui.colored_label(egui::Color32::from_rgb(217, 82, 20), "Warning: value is outside the field range. Choose a value within the range above.");
            } else if self.slice_low == self.slice_high {
                ui.label("Constant field: no distinct isosurface exists.");
            }
            if self.render_mode == 1 && ui.add_enabled(!invalid && self.slice_low != self.slice_high, egui::Button::new("Update Surface").min_size(egui::vec2(ui.available_width(), 0.0))).clicked() {
                self.generate_surface = true;
            }
        });
        ui.label("Transparency");
        let mut transparency = 1.0 - self.opacity;
        if ui.add(crate::dialog_widgets::Slider::new(&mut transparency, 0.0..=1.0).fixed_decimals(2)).changed() {
            self.opacity = 1.0 - transparency;
        }
        crate::dialog_widgets::checkbox(ui, &mut self.show_slice_plane, "Show slice plane");
        if ui.button("Fit View").clicked() {
            self.reset_camera();
        }
        ui.label("Drag: orbit | Right drag: pan | Wheel: zoom");
        let used = ui.cursor().top() - top;
        let size = egui::vec2(ui.available_width().max(100.0), (height - used - 52.0).max(150.0));
        let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click_and_drag());
        if response.hovered() {
            let wheel = ui.input(|i| i.smooth_scroll_delta.y) / 50.0;
            if wheel != 0.0 {
                self.zoom = (self.zoom * (wheel * 0.12).exp()).clamp(0.15, 15.0);
            }
            if response.double_clicked() {
                self.reset_camera();
            }
        }
        if response.dragged_by(egui::PointerButton::Primary) {
            let d = response.drag_delta();
            self.yaw = (self.yaw + d.x * 0.008 + std::f32::consts::PI).rem_euclid(std::f32::consts::TAU) - std::f32::consts::PI;
            self.pitch = (self.pitch + d.y * 0.008).clamp(-1.55, 1.55);
        }
        if response.dragged_by(egui::PointerButton::Secondary) {
            let d = response.drag_delta();
            self.pan[0] += d.x * 2.5 / (size.y * self.zoom);
            self.pan[1] -= d.y * 2.5 / (size.y * self.zoom);
        }
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(245, 247, 251));
        let aspect = size.x / size.y;
        let camera = self.camera();
        let to_screen = |ndc: egui::Vec2| egui::pos2(rect.left() + (ndc.x + 1.0) * size.x * 0.5, rect.top() + (1.0 - ndc.y) * size.y * 0.5);
        if self.render_mode == 0 || !self.surface.vertices.is_empty() {
            if !self.colour_low.is_finite() || !self.colour_high.is_finite() || self.colour_high < self.colour_low {
                self.error = "Color range must be finite with maximum >= minimum.".into();
            } else if self.render_mode == 0 {
                if self.volume_dirty {
                    self.view.quality = self.quality;
                    let grid = self.volume.fields[self.selected as usize].clone();
                    if let Err(e) = self.view.set_volume(&grid) {
                        self.error = e;
                    }
                    self.volume_dirty = false;
                }
                if !self.view.volume.is_empty() {
                    // The image is ray marched at up to 480 pixels across and scaled to the view.
                    let scale = (480.0 / size.x).min(1.0);
                    let (w, h) = (((size.x * scale) as usize).max(16), ((size.y * scale) as usize).max(16));
                    let threshold = if self.surface_level.is_finite() { self.surface_level } else { self.slice_high };
                    let key = format!("{w} {h} {:?} {} {threshold} {} {} {} {} {:p}", (camera.yaw, camera.pitch, camera.zoom, camera.pan), self.opacity, self.slice_low, self.slice_high, self.palette, self.quality, self.view.volume.as_ptr());
                    if self.view.texture.as_ref().is_none_or(|t| t.1 != key) {
                        let pixels = self.view.render_volume(w, h, camera, self.opacity, threshold, self.slice_low, self.slice_high, self.palette);
                        let image = egui::ColorImage::from_rgba_unmultiplied([w, h], &pixels);
                        let handle = ui.ctx().load_texture("electronic volume", image, egui::TextureOptions::LINEAR);
                        self.view.texture = Some((handle, key));
                    }
                    if let Some((t, _)) = &self.view.texture {
                        painter.image(t.id(), rect, egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)), egui::Color32::WHITE);
                    }
                }
            } else {
                let tris = self.view.surface_triangles(aspect, camera, self.opacity, self.colour_low, self.colour_high, self.palette, self.specular, self.shininess);
                let mut mesh = egui::Mesh::default();
                for t in &tris {
                    for (p, c) in t {
                        let base = mesh.vertices.len() as u32;
                        mesh.colored_vertex(to_screen(*p), colour32(Vec3::new(c[0], c[1], c[2]), c[3]));
                        mesh.indices.push(base);
                    }
                }
                painter.add(egui::Shape::mesh(mesh));
            }
            if self.show_slice_plane {
                let grid = &self.volume.fields[self.selected as usize];
                if let Ok(plane) = self.slice.plane(grid) {
                    let points: Vec<egui::Pos2> = plane.boundary.iter().map(|&p| to_screen(self.view.project(p, self.render_mode == 0, aspect, camera))).collect();
                    if points.len() >= 3 {
                        painter.add(egui::Shape::convex_polygon(points.clone(), egui::Color32::from_rgba_unmultiplied(65, 135, 200, 48), egui::Stroke::NONE));
                        painter.add(egui::Shape::closed_line(points, egui::Stroke::new(1.5, egui::Color32::from_rgba_unmultiplied(35, 95, 155, 220))));
                    }
                }
            }
            let caption = if self.render_mode == 0 { "Volume density".to_string() } else { format!("{} triangles", self.surface.vertices.len() / 3) };
            painter.text(rect.min + egui::vec2(12.0, 12.0), egui::Align2::LEFT_TOP, caption, egui::FontId::proportional(13.0), egui::Color32::from_rgb(70, 80, 95));
        } else {
            painter.text(rect.min + egui::vec2(20.0, 25.0), egui::Align2::LEFT_TOP, "Set the 3D isovalue and click Update surface.", egui::FontId::proportional(13.0), egui::Color32::from_rgb(90, 100, 115));
        }
        // The colour bar.
        let (bar, _) = ui.allocate_exact_size(egui::vec2(size.x, 12.0), egui::Sense::hover());
        for i in 0..128 {
            let c = palette_colour(i as f32 / 127.0, self.palette);
            let x0 = bar.left() + bar.width() * i as f32 / 128.0;
            let x1 = bar.left() + bar.width() * (i + 1) as f32 / 128.0;
            ui.painter().rect_filled(egui::Rect::from_min_max(egui::pos2(x0, bar.top()), egui::pos2(x1, bar.top() + 10.0)), 0.0, colour32(c, 1.0));
        }
        let (low, high, unit) = if self.render_mode == 0 { (self.slice_low, self.slice_high, grid_unit) } else { (self.colour_low, self.colour_high, self.colour_unit.clone()) };
        ui.horizontal(|ui| {
            ui.label(format_g(low as f64, 5));
            ui.add_space((size.x * 0.4 - 60.0).max(0.0));
            ui.label(unit_display(&unit));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(format_g(high as f64, 5));
            });
        });
    }

    fn view_2d(&mut self, ui: &mut egui::Ui, height: f32) {
        let top = ui.cursor().top();
        let grid = self.volume.fields[self.selected as usize].clone();
        ui.label("2D section");
        let s = &mut self.slice;
        let planes = ["bc plane (normal a*)", "ac plane (normal b*)", "ab plane (normal c*)", "Arbitrary plane"];
        egui::ComboBox::from_id_salt("section plane").width(ui.available_width()).selected_text(planes[s.axis as usize]).show_ui(ui, |ui| {
            for (i, name) in planes.iter().enumerate() {
                s.dirty |= ui.selectable_value(&mut s.axis, i as i32, *name).changed();
            }
        });
        if s.axis == 3 {
            ui.label("Point (fractional a, b, c)");
            ui.horizontal(|ui| {
                for c in &mut s.point {
                    s.dirty |= ui.add(egui::DragValue::new(c).speed(0.01).custom_formatter(|v, _| format_g(v, 4))).changed();
                }
            });
            ui.label("Normal (Cartesian x, y, z)");
            ui.horizontal(|ui| {
                for c in &mut s.normal {
                    s.dirty |= ui.add(egui::DragValue::new(c).speed(0.01).custom_formatter(|v, _| format_g(v, 4))).changed();
                }
            });
        } else {
            ui.label("Slice position (fraction of cell)");
            s.dirty |= ui.add(crate::dialog_widgets::Slider::new(&mut s.position, 0.0..=1.0).fixed_decimals(3)).changed();
        }
        let geometry = match s.update(&grid) {
            Ok(g) => g,
            Err(e) => {
                ui.label(e);
                return;
            }
        };
        if geometry.boundary.len() < 3 {
            ui.label("Warning: this plane does not intersect the cell. Change the point or normal.");
            return;
        }
        crate::dialog_widgets::checkbox(ui, &mut s.show_contour, "Show contour line");
        ui.label("2D isovalue (contour)");
        let mut level_changed = ui.add_sized([ui.available_width(), 18.0], egui::DragValue::new(&mut s.level).speed(0.001).custom_formatter(|v, _| format_g(v, 5))).changed();
        if ui.button("Estimate 2D level").clicked() {
            if let Ok(r) = display_range(&s.visible) {
                s.level = r.suggested as f32;
                level_changed = true;
            }
        }
        ui.horizontal(|ui| {
            if ui.button("Reset 2D").clicked() {
                s.reset();
            }
            ui.add(crate::dialog_widgets::Slider::new(&mut s.angle, -180.0..=180.0).fixed_decimals(0).suffix("°"));
        });
        ui.label("Drag: rotate in plane | Right drag: pan | Wheel: zoom");
        if !s.level.is_finite() {
            ui.label("Enter a finite 2D isovalue.");
            return;
        }
        if s.dirty || level_changed {
            if let Err(e) = s.recontour() {
                ui.label(e);
                return;
            }
        }
        let low = s.visible.iter().copied().fold(f64::INFINITY, f64::min);
        let high = s.visible.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if (s.level as f64) < low || (s.level as f64) > high {
            ui.colored_label(
                egui::Color32::from_rgb(217, 82, 20),
                format!("Warning: contour value is outside this slice's range ({} to {} {}).", format_g(low, 5), format_g(high, 5), unit_display(&grid.unit)),
            );
        } else if low == high {
            ui.label("Constant slice: no distinct contour exists.");
        }
        ui.label(format!("Range: {} to {} {}", format_g(self.slice_low as f64, 5), format_g(self.slice_high as f64, 5), unit_display(&grid.unit)));
        let used = ui.cursor().top() - top;
        let size = egui::vec2(ui.available_width().max(80.0), (height - used).max(80.0));
        let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click_and_drag());
        let s = &mut self.slice;
        if response.hovered() {
            let wheel = ui.input(|i| i.smooth_scroll_delta.y) / 50.0;
            if wheel != 0.0 {
                s.zoom = (s.zoom * (wheel * 0.12).exp()).clamp(0.1, 20.0);
            }
            if response.double_clicked() {
                s.reset();
            }
        }
        if response.dragged_by(egui::PointerButton::Primary) {
            s.angle = (s.angle + response.drag_delta().x * 0.5 + 180.0).rem_euclid(360.0) - 180.0;
        }
        if response.dragged_by(egui::PointerButton::Secondary) {
            let d = response.drag_delta();
            s.pan[0] += d.x / size.x;
            s.pan[1] += d.y / size.y;
        }
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(245, 247, 251));
        let Some(slice) = s.slice.clone() else { return };
        let screen = s.screen(&slice, size);
        let (lo, hi) = (self.slice_low, self.slice_high);
        let mut mesh = egui::Mesh::default();
        for (polygon, value) in s.cells(&grid) {
            let t = if hi > lo { ((value - lo as f64) / (hi as f64 - lo as f64)) as f32 } else { 0.5 };
            let c = colour32(palette_colour(t, self.palette), 1.0);
            let base = mesh.vertices.len() as u32;
            for p in &polygon {
                mesh.colored_vertex(rect.min + screen(*p), c);
            }
            for k in 1..polygon.len() as u32 - 1 {
                mesh.add_triangle(base, base + k, base + k + 1);
            }
        }
        painter.add(egui::Shape::mesh(mesh));
        if s.show_contour {
            for (a, b) in s.segments(&grid) {
                painter.line_segment([rect.min + screen(a), rect.min + screen(b)], egui::Stroke::new(1.0, egui::Color32::from_rgb(30, 40, 55)));
            }
        }
    }

    fn poll(&mut self) {
        if let Some(rx) = &self.task {
            if let Ok(r) = rx.try_recv() {
                self.task = None;
                self.apply(r);
            }
        }
    }
}

fn floats(args: &str) -> Vec<f32> {
    args.split_whitespace().filter_map(|w| w.parse().ok()).collect()
}

fn on(args: &str) -> bool {
    matches!(args.trim(), "on" | "true" | "1")
}

/// The operation index of a tool name (any group) or a number.
fn operation_index(name: &str) -> Option<i32> {
    if let Ok(n) = name.parse::<i32>() {
        return (0..=37).contains(&n).then_some(n);
    }
    let lower = name.to_ascii_lowercase();
    let find = |list: &[&str], first: i32| list.iter().position(|l| l.to_ascii_lowercase() == lower).map(|i| first + i as i32);
    find(&GENERAL_TOOLS, 0).or_else(|| find(&CHARGE_TOOLS, 25)).or_else(|| find(&BONDING_TOOLS, 33))
}

impl StructureDialog for ElectronicDialog {
    fn close(&mut self) {
        self.open = false;
    }
    fn drop_file(&mut self, path: &str) -> bool {
        let reference = self.drop_reference;
        self.load(path, reference, false);
        true
    }
    fn workspace(&self) -> Option<af_electronic::project::Workspace> {
        Some(self.snapshot(true))
    }
    fn restore_workspace(&mut self, saved: &af_electronic::project::Workspace) -> Result<(), String> {
        self.restore(saved)
    }
    fn busy(&self) -> bool {
        self.task.is_some()
    }
    fn steps(&self) -> &'static [&'static str] {
        &["electronic"]
    }

    fn title(&self) -> &'static str {
        "Electronic Post-Processing"
    }

    fn open(&mut self, _step: &str) {
        self.open = true;
    }

    fn is_open(&self) -> bool {
        self.open
    }

    fn show(&mut self, ctx: &egui::Context, _input: &Input) -> Option<Outcome> {
        self.poll();
        if self.task.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
        }
        if !self.open {
            return None;
        }
        if self.generate_surface && self.task.is_none() {
            self.generate_surface = false;
            self.calculate(22, false);
        }
        let mut open = self.open;
        let screen = ctx.content_rect();
        crate::theme::window("Electronic Post-Processing")
            .id(egui::Id::unique("electronic post-processing"))
            .open(&mut open)
            .collapsible(false)
            .default_size([1380.0f32.min(screen.width() - 20.0), 860.0f32.min(screen.height() - 60.0)])
            .min_size([600.0, 480.0])
            .show(ctx, |ui| {
                crate::analysis_panels::two_panels(ui, "electronic", self, |s, ui| s.controls(ui), |s, ui| s.results(ui));
            });
        self.open = open;
        if let Some(path) = self.picker.show(ctx) {
            let path = path.to_string_lossy().into_owned();
            let r = match self.picker_action {
                4 => self.save_workspace(&path),
                3 => self.open_workspace(&path),
                2 => {
                    self.save(&path);
                    Ok(())
                }
                action => {
                    self.load(&path, action == 1, false);
                    Ok(())
                }
            };
            if let Err(e) = r {
                self.error = e;
            }
        }
        None
    }

    /// try COMMAND ARGS | load PATH | reference PATH | quantity auto|density|potential|elf|raw | units bohr|angstrom |
    /// clear | tool NAME|N | field N | operand N | local-reference on|off |
    /// set NAME VALUE (scalar, sigma, radius, count, axis, window, mask-low, mask-high, weight,
    /// boolean, floor, start/end/u/v X Y Z, alpha, real-cutoff, reciprocal-cutoff, colour-surface,
    /// append, level, opacity, quality, palette, render volume|surface, layout N, reflections A; B,
    /// charges ..., thresholds A; B) | calculate | surface | undo | redo | format N | export PATH |
    /// save-workspace PATH | open-workspace PATH | slice AXIS [POSITION] | slice-arbitrary PX PY PZ NX NY NZ |
    /// contour on|off [LEVEL] | expect-fields N | expect-field-name TEXT | expect-table ROW COL VALUE |
    /// expect-rows N | expect-heading TEXT | expect-triangles MIN | expect-range LOW HIGH |
    /// expect-error TEXT | expect-no-error | expect-history N | expect-contours MIN
    fn script(&mut self, command: &str, args: &str, _input: &Input) -> Result<Option<Outcome>, String> {
        let args = args.trim();
        let n = floats(args);
        let num = |i: usize| n.get(i).copied().ok_or(format!("{command} wants a number"));
        let fail_on_error = |s: &Self| if s.error.is_empty() || s.error.starts_with("Saved") { Ok(None) } else { Err(s.error.clone()) };
        match command {
            "try" => {
                // try COMMAND ARGS: runs a command whose failure is checked by expect-error.
                let (inner, rest) = args.split_once(' ').unwrap_or((args, ""));
                if let Err(e) = self.script(inner, rest, _input) {
                    self.error = e;
                }
            }
            "load" | "reference" => {
                self.load(args, command == "reference", true);
                return fail_on_error(self);
            }
            "quantity" => self.quantity = QUANTITY_NAMES.iter().position(|q| *q == args).ok_or("quantity auto|density|potential|elf|raw")? as i32,
            "units" => self.cube_units = (args == "angstrom") as i32,
            "clear" => self.clear(),
            "tool" => {
                let op = operation_index(args).ok_or(format!("no tool {args}"))?;
                self.operation = op;
                self.tool_group = if op >= 33 { 2 } else if op >= 25 { 1 } else { 0 };
            }
            "field" => {
                let i = num(0)? as i32;
                if i < 0 || i as usize >= self.volume.fields.len() {
                    return Err(format!("no field {i}: {} fields", self.volume.fields.len()));
                }
                self.selected = i;
            }
            "operand" => {
                let i = num(0)? as i32;
                if i < 0 || i as usize >= self.reference_volume().fields.len() {
                    return Err(format!("no reference field {i}"));
                }
                self.reference = i;
            }
            "local-reference" => self.local_reference = on(args),
            "set" => {
                let (name, value) = args.split_once(' ').unwrap_or((args, ""));
                let value = value.trim();
                let v = floats(value);
                let one = || v.first().copied().ok_or(format!("set {name} wants a number"));
                let three = || if v.len() == 3 { Ok([v[0], v[1], v[2]]) } else { Err(format!("set {name} wants X Y Z")) };
                match name {
                    "scalar" => self.scalar = one()?,
                    "sigma" => self.sigma = one()?,
                    "radius" => self.radius = one()? as i32,
                    "count" => self.count = one()? as i32,
                    "axis" => self.axis = one()? as i32,
                    "window" => self.window = one()? as i32,
                    "mask-low" => self.mask_low = one()?,
                    "mask-high" => self.mask_high = one()?,
                    "weight" => self.reference_weight = one()?,
                    "boolean" => self.boolean_operation = BOOLEAN_NAMES.iter().position(|b| *b == value).ok_or("boolean union|intersection|difference|xor")? as i32,
                    "floor" => self.bonding_floor = one()?,
                    "start" => self.start = three()?,
                    "end" => self.end = three()?,
                    "u" => self.u = three()?,
                    "v" => self.v = three()?,
                    "alpha" => self.alpha = one()?,
                    "real-cutoff" => self.real_cutoff = one()?,
                    "reciprocal-cutoff" => self.reciprocal_cutoff = one()?,
                    "colour-surface" => self.colour_surface = on(value),
                    "append" => self.append_surface = on(value),
                    "level" => self.surface_level = one()?,
                    "opacity" => self.opacity = one()?.clamp(0.0, 1.0),
                    "quality" => {
                        self.quality = (one()? as i32).clamp(0, 3);
                        self.view.quality = self.quality;
                        self.slice.set_resolution(SLICE_RESOLUTIONS[self.quality as usize]);
                        self.volume_dirty = true;
                    }
                    "palette" => self.palette = (one()? as i32).clamp(0, 2),
                    "render" => self.render_mode = (value == "surface") as i32,
                    "layout" => self.view_layout = (one()? as i32).clamp(0, 2),
                    "reflections" => self.reflections = value.split(';').map(str::trim).collect::<Vec<_>>().join("\n"),
                    "charges" => self.charges = value.to_string(),
                    "thresholds" => self.thresholds = value.split(';').map(str::trim).collect::<Vec<_>>().join("\n"),
                    _ => return Err(format!("unknown electronic setting {name}")),
                }
            }
            "calculate" => {
                self.calculate(self.operation, true);
                return fail_on_error(self);
            }
            "surface" => {
                self.refresh_field();
                if let Some(&level) = n.first() {
                    self.surface_level = level;
                }
                self.render_mode = 1;
                self.calculate(22, true);
                return fail_on_error(self);
            }
            "undo" => self.undo_calculation()?,
            "redo" => self.redo_calculation()?,
            "format" => self.set_export_format((num(0)? as i32).clamp(0, 7)),
            "export" => {
                let ext = std::path::Path::new(args).extension().map(|e| format!(".{}", e.to_string_lossy().to_ascii_lowercase())).unwrap_or_default();
                // The format from the extension unless the selected one already writes it (OBJ/PLY, PNG/SVG).
                if EXPORT_EXTENSIONS.get(self.export_format as usize) != Some(&ext.as_str()) {
                    if let Some(i) = EXPORT_EXTENSIONS.iter().position(|e| *e == ext) {
                        self.export_format = i as i32;
                    }
                }
                self.refresh_field();
                if self.export_format == 7 {
                    let grid = self.volume.fields.get(self.selected as usize).cloned().ok_or("Load a volume first.")?;
                    self.slice.update(&grid)?;
                }
                self.save(args);
                return fail_on_error(self);
            }
            "save-workspace" => self.save_workspace(args)?,
            "open-workspace" => self.open_workspace(args)?,
            "slice" => {
                self.slice.axis = (num(0)? as i32).clamp(0, 2);
                if let Some(&p) = n.get(1) {
                    self.slice.position = p.clamp(0.0, 1.0);
                }
                self.slice.dirty = true;
            }
            "slice-arbitrary" => {
                if n.len() != 6 {
                    return Err("slice-arbitrary PX PY PZ NX NY NZ".into());
                }
                self.slice.axis = 3;
                self.slice.point = [n[0], n[1], n[2]];
                self.slice.normal = [n[3], n[4], n[5]];
                self.slice.dirty = true;
            }
            "contour" => {
                let w: Vec<&str> = args.split_whitespace().collect();
                self.slice.show_contour = w.first().is_some_and(|w| on(w));
                if let Some(level) = w.get(1).and_then(|l| l.parse::<f32>().ok()) {
                    self.slice.level = level;
                    self.slice.estimate_level = false;
                }
                self.refresh_field();
                let grid = self.volume.fields.get(self.selected as usize).cloned().ok_or("Load a volume first.")?;
                self.slice.update(&grid)?;
                self.slice.recontour()?;
            }
            "expect-fields" => {
                if self.volume.fields.len() != num(0)? as usize {
                    return Err(format!("{} fields ({})", self.volume.fields.len(), self.error));
                }
            }
            "expect-field-name" => {
                let name = self.volume.fields.get(self.selected as usize).map(|f| f.name.as_str()).unwrap_or("");
                if name != args {
                    return Err(format!("the selected field is {name:?}"));
                }
            }
            "expect-rows" => {
                let rows = self.result.table.len() / self.result.columns.max(1) as usize;
                if rows != num(0)? as usize {
                    return Err(format!("{rows} table rows"));
                }
            }
            "expect-table" => {
                let w: Vec<f64> = args.split_whitespace().filter_map(|w| w.parse().ok()).collect();
                if w.len() != 3 {
                    return Err("expect-table ROW COL VALUE".into());
                }
                let i = w[0] as usize * self.result.columns.max(1) as usize + w[1] as usize;
                let got = self.result.table.get(i).copied().ok_or(format!("no table cell {} {}", w[0], w[1]))?;
                if (got - w[2]).abs() > 1e-6 * w[2].abs().max(1.0) {
                    return Err(format!("table cell {} {} is {got}", w[0], w[1]));
                }
            }
            "expect-heading" => {
                if self.result.heading != args {
                    return Err(format!("the table heading is {:?}", self.result.heading));
                }
            }
            "expect-triangles" => {
                let t = self.surface.vertices.len() / 3;
                if t < num(0)? as usize {
                    return Err(format!("{t} surface triangles ({})", self.error));
                }
            }
            "expect-range" => {
                self.refresh_field();
                if n.len() != 2 || (self.slice_low - n[0]).abs() > 1e-4 * n[0].abs().max(1.0) || (self.slice_high - n[1]).abs() > 1e-4 * n[1].abs().max(1.0) {
                    return Err(format!("field range {} to {}", self.slice_low, self.slice_high));
                }
            }
            "expect-contours" => {
                let segments = self.slice.contours.len() / 2;
                if segments < num(0)? as usize {
                    return Err(format!("{segments} contour segments"));
                }
            }
            "expect-history" => {
                if self.history.len() != num(0)? as usize {
                    return Err(format!("{} history entries: {:?}", self.history.len(), self.history));
                }
            }
            "expect-error" => {
                if self.error.is_empty() || !self.error.contains(args) {
                    return Err(format!("error is {:?}", self.error));
                }
            }
            "expect-no-error" => {
                if !self.error.is_empty() && !self.error.starts_with("Saved") {
                    return Err(self.error.clone());
                }
            }
            other => return Err(format!("unknown electronic command {other}")),
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reflections_parse_like_the_cpp_streams() {
        let (h, r) = parse_reflections("0 0 0\n\n 1 0 0\r\n-1 0 0", false).unwrap();
        assert_eq!(h, vec![IVec3::ZERO, IVec3::X, -IVec3::X]);
        assert_eq!(r.len(), 3);
        let (_, r) = parse_reflections("1 0 0 0.5 -0.25", true).unwrap();
        assert_eq!((r[0].value.re, r[0].value.im), (0.5, -0.25));
        assert_eq!(parse_reflections("1 0", false).unwrap_err(), "Malformed reflection row");
        assert_eq!(parse_reflections("1 0 0 2", false).unwrap_err(), "Extra reflection columns");
        assert_eq!(parse_reflections("1 0 0 2", true).unwrap_err(), "Malformed reflection row");
    }

    #[test]
    fn tools_are_found_by_name_or_number() {
        assert_eq!(operation_index("Total integral"), Some(0));
        assert_eq!(operation_index("bader partition"), Some(37));
        assert_eq!(operation_index("Density difference"), Some(25));
        assert_eq!(operation_index("12"), Some(12));
        assert_eq!(operation_index("38"), None);
    }

    #[test]
    fn rasterised_triangle_covers_its_area() {
        let red = [1.0, 0.0, 0.0, 1.0];
        let tri = [[(egui::vec2(-1.0, -1.0), red), (egui::vec2(1.0, -1.0), red), (egui::vec2(-1.0, 1.0), red)]];
        let image = rasterize(10, 10, &tri);
        let reds = image.chunks_exact(4).filter(|p| p[0] == 255 && p[1] == 0).count();
        assert!((40..=60).contains(&reds), "{reds}");
    }
}
