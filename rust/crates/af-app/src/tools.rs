//! The scientific tools of the Analysis menu (the C++ ui/ScientificToolsDialog,
//! ui/ScientificToolsLayouts and ui/SurfaceMeshOverlay): one window per tool
//! family (per-atom analyses beside the view, plot workspaces, trajectory
//! analyses, property calculators, simulation runs and input generators),
//! inputs from the active structure, files or typed values, background runs
//! whose results (report, summary, tensors, plots, per-atom values, generated
//! files, structures) live in a temporary run folder, kept runs for comparison,
//! saving, a result mesh drawn translucent over the view, and pipeline step
//! editing for relax, nvt-dynamics, npt-dynamics and standardize-cell.

use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use af_core::format::format_g;
use af_core::json::Json;
use af_core::structure::Structure;
use af_render::scene::FrameView;
use af_science::atom_properties::{per_atom_properties, AtomProperty};
use af_science::catalog::{find_tool, ParameterDef, ToolDef, ToolView};
use af_science::data::{read_frames, structure_json, write_extxyz};
use af_science::plots::{overlay_plots, plot_csv, result_plots, PlotSpec};
use af_science::tools::{result_document, result_report, run_tool};
use eframe::egui;
use glam::{Vec3, Vec4};

use crate::file_browser::{FileBrowser, Purpose};
use crate::plot_widget::unit_display;
use crate::structure_dialog::StepEdit;

/// What the tool window gives the app.
#[derive(Default)]
pub struct ToolEvents {
    /// Structures to open in new tabs.
    pub opened: Vec<(Structure, String)>,
    /// Colour the atoms by a per-atom result.
    pub colour: Option<AtomProperty>,
}

// ---------------------------------------------------------------------------
// Files of a run
// ---------------------------------------------------------------------------

fn write_text(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|_| format!("Cannot save {}", path.display()))
}

fn read_text(path: &Path) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|_| format!("Cannot read {}", path.display()))
}

fn stamp() -> u128 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0)
}

fn science_folder() -> PathBuf {
    std::env::temp_dir().join("AtomForge-science")
}

fn new_run_directory() -> Result<PathBuf, String> {
    let directory = science_folder().join(format!("run-{}", stamp()));
    std::fs::create_dir_all(&directory).map_err(|e| format!("Cannot create {}: {e}", directory.display()))?;
    Ok(directory)
}

fn generated_files(result: &Json) -> Vec<(String, String)> {
    match result.find("files") {
        Some(Json::Object(members)) => members.iter().map(|(name, text)| (name.clone(), text.string().unwrap_or("").to_string())).collect(),
        _ => Vec::new(),
    }
}

/// Writes tool input from the view to a new file in the scientific tools'
/// temporary folder (large inputs go through files).
fn save_input_file(kind: &str, text: &str) -> Result<String, String> {
    let root = science_folder();
    std::fs::create_dir_all(&root).map_err(|_| format!("Cannot save the {kind} for analysis"))?;
    let file = root.join(format!("{kind}-{}.json", stamp()));
    std::fs::write(&file, text).map_err(|_| format!("Cannot save the {kind} for analysis"))?;
    Ok(file.to_string_lossy().into_owned())
}

fn save_active_structure(structure: &Structure) -> Result<String, String> {
    save_input_file("structure", &structure_json(structure).dump(-1))
}

/// The view's per-atom property as a JSON list (NaN for atoms without a value).
fn save_active_property(structure: &Structure) -> Result<String, String> {
    let values: Vec<String> = structure.atom_property.iter().map(|&v| if v.is_finite() { format_g(v, 17) } else { "NaN".into() }).collect();
    save_input_file("property", &format!("[{}]\n", values.join(",")))
}

/// The atomic numbers of the view's atoms, for species inputs.
fn save_active_species(structure: &Structure) -> Result<String, String> {
    let species = Json::Array(structure.atoms.iter().map(|a| Json::Number(a.atomic_number as f64)).collect());
    save_input_file("species", &species.dump(-1))
}

// ---------------------------------------------------------------------------
// Potentials, steps and labels
// ---------------------------------------------------------------------------

/// Index order of the interatomic-potential selector.
const EMT: usize = 0;
const LENNARD_JONES: usize = 1;
const EAM: usize = 2;
const TERSOFF: usize = 3;
const STILLINGER_WEBER: usize = 4;
const BUCKINGHAM: usize = 5;
const POTENTIALS: [&str; 6] =
    ["EMT (effective-medium theory)", "Lennard-Jones 12-6", "EAM / Finnis-Sinclair file", "Tersoff (covalent)", "Stillinger-Weber (covalent)", "Buckingham + Coulomb (ionic)"];
const EAM_FORMATS: [&str; 4] = ["auto", "setfl", "fs", "funcfl"];
const EAM_FORMAT_LABELS: [&str; 4] = ["Detect from file name", "setfl (eam/alloy)", "Finnis-Sinclair (eam/fs)", "funcfl (single element eam)"];
const DEFAULT_BUCKINGHAM: &str = "{\"charges\": {\"Mg\": 2, \"O\": -2}, \"cutoff\": 10,\n \"pairs\": [{\"elements\": [\"Mg\", \"O\"], \"A\": 821.6, \"rho\": 0.3242, \"C\": 0},\n           {\"elements\": [\"O\", \"O\"], \"A\": 22764, \"rho\": 0.149, \"C\": 27.88}]}";

fn potential_json(potential: usize, epsilon: f64, sigma: f64, cutoff: f64, file: &str, format: usize, options: &str) -> Result<Json, String> {
    let mut result = Json::object();
    if potential == TERSOFF || potential == STILLINGER_WEBER {
        result.set("potential", if potential == TERSOFF { "Tersoff" } else { "StillingerWeber" });
        if !file.is_empty() {
            result.set("file", file); // empty: built-in Si
        }
        return Ok(result);
    }
    if potential == BUCKINGHAM {
        result = Json::parse(options).map_err(|e| format!("Buckingham parameters: {e}"))?;
        if !result.is_object() {
            return Err("Buckingham parameters must be a JSON object".into());
        }
        result.set("potential", "Buckingham");
        return Ok(result);
    }
    if potential == EMT {
        result.set("potential", "EMT");
        return Ok(result);
    }
    if potential == EAM {
        if file.is_empty() {
            return Err("Choose an EAM potential file".into());
        }
        result.set("potential", "EAM");
        result.set("file", file);
        result.set("format", EAM_FORMATS[format.min(3)]);
        return Ok(result);
    }
    result.set("potential", "LennardJones");
    result.set("epsilon", epsilon);
    result.set("sigma", sigma);
    result.set("cutoff", cutoff);
    Ok(result)
}

/// Pipeline steps that run a catalog tool.
fn step_tool(step: &str) -> Option<&'static str> {
    match step {
        "relax" => Some("relax"),
        "nvt-dynamics" => Some("nvt"),
        "npt-dynamics" => Some("npt"),
        "standardize-cell" => Some("symmetry"),
        _ => None,
    }
}

const STEP_CELLS: [&str; 3] = ["primitive", "conventional", "symmetrized"];

fn scalar_kind(kind: &str) -> bool {
    matches!(kind, "float" | "int" | "bool")
}

fn trimmed(text: &str) -> String {
    text.trim_matches([' ', '\t', '\r', '\n']).to_string()
}

/// A text field shows a JSON string literal ("conventional") without its quotes.
fn field_text(kind: &str, text: &str) -> String {
    if kind != "string" || text.len() < 2 || !text.starts_with('"') || !text.ends_with('"') {
        return text.to_string();
    }
    match Json::parse(text) {
        Ok(Json::String(s)) => s,
        _ => text.to_string(),
    }
}

/// Title case for window titles (the C++ uiText::titleCase).
pub fn title_case(text: &str) -> String {
    const SMALL: [&str; 18] = ["a", "an", "the", "and", "or", "to", "in", "of", "as", "for", "on", "at", "by", "from", "with", "into", "per", "vs"];
    let words: Vec<&str> = text.split(' ').collect();
    let n = words.len();
    words
        .iter()
        .enumerate()
        .map(|(k, word)| {
            let plain = word.as_bytes().first().is_some_and(|c| c.is_ascii_lowercase()) && word.bytes().all(|c| c.is_ascii_alphabetic() || c == b'-');
            let small = SMALL.contains(word);
            if plain && (!small || k == 0 || k + 1 == n) {
                let mut w = word.to_string();
                w[..1].make_ascii_uppercase();
                w
            } else {
                word.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// A result field name as a label: trailing unit words become a suffix
/// (max_force_eV_per_A -> "Max force (eV/A)").
pub fn pretty_label(key: &str) -> String {
    const UNITS: [(&str, &str); 17] = [
        ("eV", "eV"),
        ("A", "A"),
        ("A2", "A^2"),
        ("A3", "A^3"),
        ("K", "K"),
        ("GPa", "GPa"),
        ("fs", "fs"),
        ("ps", "ps"),
        ("deg", "deg"),
        ("THz", "THz"),
        ("amu", "amu"),
        ("kV", "kV"),
        ("mrad", "mrad"),
        ("atom", "atom"),
        ("inv", "1/"),
        ("cm", "cm"),
        ("meV", "meV"),
    ];
    let tokens: Vec<&str> = key.split('_').filter(|t| !t.is_empty()).collect();
    let unit_of = |w: &str| UNITS.iter().find(|(s, _)| *s == w).map(|(_, shown)| *shown);
    let mut first = tokens.len();
    while first > 1 && (unit_of(tokens[first - 1]).is_some() || tokens[first - 1] == "per") {
        first -= 1;
    }
    while first < tokens.len() && tokens[first] == "per" {
        first += 1;
    }
    let mut label = String::new();
    let mut unit = String::new();
    for (k, t) in tokens.iter().enumerate() {
        if k < first {
            if !label.is_empty() {
                label.push(' ');
            }
            label.push_str(t);
        } else if *t == "per" {
            unit.push('/');
        } else {
            unit.push_str(unit_of(t).unwrap_or(""));
        }
    }
    if let Some(c) = label.get(..1) {
        let upper = c.to_ascii_uppercase();
        label.replace_range(..1, &upper);
    }
    if unit.is_empty() {
        label
    } else {
        format!("{label} ({unit})")
    }
}

/// Scalars of a result as (label, text) rows.
fn summary_of(result: &Json) -> Vec<(String, String)> {
    let Json::Object(members) = result else { return Vec::new() };
    let mut rows = Vec::new();
    for (key, value) in members {
        let text = match value {
            Json::Number(n) => format_g(*n, 6),
            Json::Bool(b) => (if *b { "yes" } else { "no" }).to_string(),
            Json::String(s) if s.len() <= 80 && !s.contains('\n') => s.clone(),
            _ => continue,
        };
        rows.push((pretty_label(key), text));
    }
    rows
}

type Matrix = (String, Vec<Vec<f64>>);

/// Small numeric tables (tensors, cells), shown as grids.
fn matrices_of(result: &Json) -> Vec<Matrix> {
    let Json::Object(members) = result else { return Vec::new() };
    let mut matrices = Vec::new();
    'next: for (key, value) in members {
        let Json::Array(items) = value else { continue };
        if items.len() < 2 || items.len() > 9 {
            continue;
        }
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for row in items {
            let Json::Array(cells) = row else { continue 'next };
            if cells.len() < 2 || cells.len() > 9 || rows.first().is_some_and(|r| r.len() != cells.len()) {
                continue 'next;
            }
            let mut numbers = Vec::new();
            for item in cells {
                let Json::Number(n) = item else { continue 'next };
                numbers.push(*n);
            }
            rows.push(numbers);
        }
        // Round-off below 1e-9 of the largest entry is shown as zero.
        let largest = rows.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
        for v in rows.iter_mut().flatten() {
            if v.abs() <= 1e-9 * largest {
                *v = 0.0;
            }
        }
        matrices.push((pretty_label(key), rows));
    }
    matrices
}

fn companion_cell(parameter: &str) -> Option<&'static str> {
    match parameter {
        "positions" | "reference_sites" => Some("cell"),
        "reference" => Some("reference_cell"),
        "current" => Some("current_cell"),
        _ => None,
    }
}

fn accepts_trajectory(tool: &str, parameter: &str) -> bool {
    match tool {
        "msd" | "structure-factor" => parameter == "positions",
        "vacf" | "vibrational-spectrum" => parameter == "velocities",
        "trajectory-structure" | "trajectory-average" | "unwrap-trajectory" | "time-series" | "frame-difference" => parameter == "trajectory_file",
        _ => false,
    }
}

fn cell_json(structure: &Structure) -> String {
    Json::Array(structure.cell_vectors.iter().map(|r| Json::Array(r.iter().map(|&v| Json::Number(v)).collect())).collect()).dump(-1)
}

fn optional(p: &ParameterDef) -> bool {
    !p.required && p.value.is_empty()
}

// ---------------------------------------------------------------------------
// The result mesh over the view
// ---------------------------------------------------------------------------

/// The leading integer of an OBJ face item (std::stol of the part before '/').
fn leading_integer(text: &str) -> Option<i64> {
    let t = text.trim_start();
    let b = t.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let s = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == s {
        return None;
    }
    t[..i].parse().ok()
}

/// Triangles of a Wavefront OBJ text ("v x y z" and "f a b c ..." lines;
/// polygons are fanned, negative indices count back from the last vertex):
/// three points per triangle.
pub fn triangles_from_obj(text: &str) -> Result<Vec<[f32; 3]>, String> {
    let mut points: Vec<[f32; 3]> = Vec::new();
    let mut triangles = Vec::new();
    for line in text.lines() {
        let mut words = line.split_whitespace();
        match words.next() {
            Some("v") => {
                let v: Vec<f32> = words.take(3).map_while(|w| w.parse::<f32>().ok()).collect();
                if v.len() != 3 {
                    return Err(format!("Malformed OBJ vertex: {line}"));
                }
                points.push([v[0], v[1], v[2]]);
            }
            Some("f") => {
                let mut face = Vec::new();
                for item in words {
                    let head = item.split('/').next().unwrap_or("");
                    let index = leading_integer(head).ok_or_else(|| format!("Malformed OBJ face: {line}"))?;
                    face.push(if index < 0 { points.len() as i64 + index } else { index - 1 });
                }
                for k in 1..face.len().saturating_sub(1) {
                    for index in [face[0], face[k], face[k + 1]] {
                        if index < 0 || index as usize >= points.len() {
                            return Err("OBJ face refers to a missing vertex".into());
                        }
                        triangles.push(points[index as usize]);
                    }
                }
            }
            _ => {}
        }
    }
    Ok(triangles)
}

/// A translucent triangle mesh drawn over the structure view (the C++
/// SurfaceMeshOverlay), e.g. the surface of Surface Mesh.
struct SurfaceMesh {
    vertices: Vec<[f32; 3]>,
    visible: bool,
    opacity: f32,
    colour: [f32; 3],
}

impl Default for SurfaceMesh {
    fn default() -> Self {
        Self { vertices: Vec::new(), visible: false, opacity: 0.45, colour: [0.30, 0.62, 0.86] }
    }
}

// ---------------------------------------------------------------------------
// The window
// ---------------------------------------------------------------------------

/// One input of the tool.
#[derive(Clone, Debug)]
struct Field {
    value: String,
    enabled: bool,
    path: String,
    /// Shown instead of `path` when it holds the active structure.
    active_label: String,
    use_file: bool,
    select_column: bool,
    column: i32,
    field: String,
    potential: usize,
    epsilon: f64,
    sigma: f64,
    cutoff: f64,
    potential_file: String,
    eam_format: usize,
    potential_options: String,
}

impl Default for Field {
    fn default() -> Self {
        Self {
            value: String::new(),
            enabled: true,
            path: String::new(),
            active_label: String::new(),
            use_file: false,
            select_column: false,
            column: 0,
            field: String::new(),
            potential: 0,
            epsilon: 0.0104,
            sigma: 3.40,
            cutoff: 8.5,
            potential_file: String::new(),
            eam_format: 0,
            potential_options: DEFAULT_BUCKINGHAM.to_string(),
        }
    }
}

/// A finished run: its files (in a temporary run folder) and what is shown.
#[derive(Clone, Default)]
struct ToolResult {
    /// result.json ({"tool", "parameters", "result"}); empty before a run.
    output: PathBuf,
    report: String,
    /// frames.extxyz; empty when the tool made no structures.
    structures: PathBuf,
    plots: Vec<PlotSpec>,
    properties: Vec<AtomProperty>,
    files: Vec<(String, String)>,
    summary: Vec<(String, String)>,
    matrices: Vec<Matrix>,
}

impl ToolResult {
    fn from_result(id: &str, output: PathBuf, structures: PathBuf, result: &Json) -> Result<Self, String> {
        Ok(Self {
            output,
            report: result_report(id, result)?,
            structures,
            plots: result_plots(id, result)?,
            properties: per_atom_properties(id, result),
            files: generated_files(result),
            summary: summary_of(result),
            matrices: matrices_of(result),
        })
    }
}

/// Runs a tool in a run folder: result.json and frames.extxyz written there.
fn execute(id: &str, request: &Json, directory: &Path) -> Result<ToolResult, String> {
    let reader = |path: &Path| af_io::loader::load_structure_from_file(&path.to_string_lossy());
    let output = run_tool(id, request, directory, Some(&reader))?;
    let result_path = directory.join("result.json");
    let structures = directory.join("frames.extxyz");
    write_text(&result_path, &(result_document(id, request, &output).dump(2) + "\n"))?;
    if !output.frames.is_empty() {
        write_extxyz(&structures, &output.frames, &output.velocities, &output.times)?;
    }
    ToolResult::from_result(id, result_path, if output.frames.is_empty() { PathBuf::new() } else { structures }, &output.result)
}

struct KeptRun {
    label: String,
    plots: Vec<PlotSpec>,
    summary: Vec<(String, String)>,
}

type Task = Arc<Mutex<Option<Result<ToolResult, String>>>>;

pub struct ToolWindow {
    pub open: bool,
    /// Update Step on a pipeline step opened here: the step with its new parameters.
    pub updated: Option<StepEdit>,
    tool: &'static str,
    fields: Vec<Field>,
    picker: FileBrowser,
    /// What the file picker is for: a field (>= 0), its potential file
    /// (<= -100), results (-3), structures (-4), all files (-5), one file
    /// (-6) or a plot's data (-10 - plot).
    picker_target: i32,
    property: usize,
    plot: usize,
    kept: Vec<KeptRun>,
    overlay: bool,
    file: usize,
    use_active: bool,
    /// Apply the active structure to the inputs on the next draw.
    active_defaults_pending: bool,
    auto_colour: bool,
    /// Trajectory loaded in Trajectory Playback (empty when none).
    trajectory: String,
    task: Option<Task>,
    error: String,
    result: ToolResult,
    step: Option<StepEdit>,
    step_cell: usize,
    mesh: SurfaceMesh,
    /// The result whose mesh the view shows.
    mesh_from: PathBuf,
    /// Events from script commands, handed to the app with the next draw.
    pending: ToolEvents,
}

impl Default for ToolWindow {
    fn default() -> Self {
        Self {
            open: false,
            updated: None,
            tool: "",
            fields: Vec::new(),
            picker: FileBrowser::default(),
            picker_target: -1,
            property: 0,
            plot: 0,
            kept: Vec::new(),
            overlay: true,
            file: 0,
            use_active: true,
            active_defaults_pending: false,
            auto_colour: true,
            trajectory: String::new(),
            task: None,
            error: String::new(),
            result: ToolResult::default(),
            step: None,
            step_cell: 0,
            mesh: SurfaceMesh::default(),
            mesh_from: PathBuf::new(),
            pending: ToolEvents::default(),
        }
    }
}

fn section(ui: &mut egui::Ui, text: &str) {
    crate::dialog_layout::separator_text(ui, text);
}

/// One component of a range or vector field (a cell of ImGui InputScalarN,
/// "%.8g"): the text is kept while it is edited; the value follows it.
fn component_field(ui: &mut egui::Ui, salt: (&str, usize), value: &mut f64, width: f32) -> bool {
    let id = ui.make_persistent_id(("component", salt.0, salt.1));
    let shown = af_core::format::format_g(*value, 8);
    let mut text = ui.data_mut(|d| d.get_temp::<String>(id)).unwrap_or_else(|| shown.clone());
    let response = ui.add_sized(egui::vec2(width, crate::theme::FRAME_HEIGHT), egui::TextEdit::singleline(&mut text).id(id.with("edit")).margin(egui::vec2(8.0, 4.0)).align(egui::Align2::LEFT_CENTER));
    let mut changed = false;
    if response.has_focus() {
        if response.changed() {
            if let Ok(v) = text.trim().parse::<f64>() {
                if v.is_finite() {
                    *value = v;
                    changed = true;
                }
            }
        }
        ui.data_mut(|d| d.insert_temp(id, text));
    } else {
        ui.data_mut(|d| d.remove::<String>(id));
    }
    changed
}

fn error_colour() -> egui::Color32 {
    egui::Color32::from_rgb(220, 80, 60)
}

impl ToolWindow {
    fn def(&self) -> Option<&'static ToolDef> {
        find_tool(self.tool)
    }

    fn running(&self) -> bool {
        self.task.is_some()
    }

    fn select_tool(&mut self, id: &'static str) {
        self.tool = id;
        self.active_defaults_pending = true;
        self.fields = find_tool(id)
            .map(|t| {
                t.parameters
                    .iter()
                    .map(|p| Field {
                        value: field_text(p.kind, p.value),
                        enabled: p.required || !p.value.is_empty(),
                        use_file: ((p.kind == "data" || p.kind == "structure") && p.value.is_empty()) || p.kind == "file",
                        ..Default::default()
                    })
                    .collect()
            })
            .unwrap_or_default();
        self.error.clear();
        self.result = ToolResult::default();
        self.kept.clear();
    }

    /// Opens the window on a catalog tool (Analysis menu, palette).
    pub fn open_tool(&mut self, id: &'static str) {
        if self.running() || find_tool(id).is_none() {
            return;
        }
        self.step = None;
        self.select_tool(id);
        self.open = true;
    }

    /// Sets a parameter's text (UI scripts): a typed value.
    pub fn set_text(&mut self, name: &str, text: &str) {
        let Some(def) = self.def() else { return };
        if let Some(i) = def.parameters.iter().position(|p| p.name == name) {
            let f = &mut self.fields[i];
            f.value = text.to_string();
            f.enabled = true;
            f.use_file = false;
            f.active_label.clear();
        }
    }

    /// The plots of the last run (UI scripts).
    pub fn plot_titles(&self) -> Vec<String> {
        self.result.plots.iter().map(|p| p.title.clone()).collect()
    }

    /// The per-atom values of the last run (UI scripts).
    pub fn properties(&self) -> &[AtomProperty] {
        &self.result.properties
    }

    fn is_source(p: &ParameterDef) -> bool {
        matches!(p.kind, "data" | "structure" | "file")
    }

    /// 0 sources (data, structures, files), 1 settings, 2 optional, 3 potential.
    fn group_of(p: &ParameterDef) -> i32 {
        if p.kind == "calculator" {
            3
        } else if optional(p) {
            2
        } else if Self::is_source(p) {
            0
        } else {
            1
        }
    }

    /// First structure or position input: filled from the active structure in the atom layout.
    fn active_input(&self) -> Option<usize> {
        self.def()?.parameters.iter().position(|p| p.kind == "structure" || (p.kind == "data" && p.name == "positions"))
    }

    fn supplied_by_active(&self, index: usize) -> bool {
        let Some(def) = self.def() else { return false };
        let Some(input) = self.active_input() else { return false };
        if !self.use_active || def.view != ToolView::Atoms {
            return false;
        }
        if index == input {
            return true;
        }
        let name = def.parameters[index].name;
        companion_cell(def.parameters[input].name).is_some_and(|c| name == c || name == "pbc")
    }

    fn members(&self, group: i32, skip: Option<usize>, hide_supplied: bool) -> Vec<usize> {
        let Some(def) = self.def() else { return Vec::new() };
        (0..def.parameters.len()).filter(|&i| Some(i) != skip && !(hide_supplied && self.supplied_by_active(i)) && Self::group_of(&def.parameters[i]) == group).collect()
    }

    fn has_fields(&self, group: i32, skip: Option<usize>) -> bool {
        !self.members(group, skip, true).is_empty()
    }

    fn use_active_structure(&mut self, index: usize, structure: &Structure) -> Result<(), String> {
        let def = self.def().ok_or("No tool is open")?;
        let path = save_active_structure(structure)?;
        let f = &mut self.fields[index];
        f.path = path;
        f.active_label = format!("Active structure ({} atoms)", structure.atoms.len());
        f.use_file = true;
        f.enabled = true;
        f.field.clear();
        f.select_column = false;
        // Species inputs take the atomic numbers of the same atoms.
        for (j, p) in def.parameters.iter().enumerate() {
            if p.name == "types" {
                let path = save_active_species(structure)?;
                let other = &mut self.fields[j];
                other.path = path;
                other.active_label = "Atomic numbers of the active structure".into();
                other.use_file = true;
                other.enabled = true;
            }
        }
        // Periodic local analyses also need the matching cell and axes.
        let companion = if def.parameters[index].kind == "data" { companion_cell(def.parameters[index].name) } else { None };
        if let Some(companion) = companion.filter(|_| structure.has_unit_cell) {
            for (j, p) in def.parameters.iter().enumerate() {
                let other = &mut self.fields[j];
                if p.name == companion {
                    other.value = cell_json(structure);
                    other.use_file = false;
                    other.enabled = true;
                } else if p.name == "pbc" {
                    other.value = "[true,true,true]".into();
                    other.enabled = true;
                }
            }
        }
        Ok(())
    }

    fn use_active_property(&mut self, index: usize, structure: &Structure) -> Result<(), String> {
        if structure.atom_property.len() != structure.atoms.len() || structure.atoms.is_empty() {
            return Err("The view has no per-atom property; compute one first (Edit > Structure Operations > Compute Per-Atom Property)".into());
        }
        let path = save_active_property(structure)?;
        let name = if structure.atom_property_name.is_empty() { "values".to_string() } else { structure.atom_property_name.clone() };
        let f = &mut self.fields[index];
        f.path = path;
        f.active_label = format!("Active property: {name} ({} atoms)", structure.atoms.len());
        f.use_file = true;
        f.enabled = true;
        f.field.clear();
        f.select_column = false;
        Ok(())
    }

    /// The request of a run (startCalculation): the active structure for the
    /// atom layout's input, then each enabled field.
    fn request(&mut self, structure: &Structure) -> Result<Json, String> {
        let def = self.def().ok_or("No tool is open")?;
        if def.view == ToolView::Atoms && self.use_active {
            if let Some(input) = self.active_input() {
                if structure.atoms.is_empty() {
                    return Err("The active structure has no atoms; load a structure or choose a file instead".into());
                }
                self.use_active_structure(input, structure)?;
            }
        }
        let mut request = Json::object();
        for (p, f) in def.parameters.iter().zip(&self.fields) {
            if !f.enabled {
                continue;
            }
            if p.kind == "calculator" {
                request.set(p.name, potential_json(f.potential, f.epsilon, f.sigma, f.cutoff, &f.potential_file, f.eam_format, &f.potential_options)?);
            } else if f.use_file {
                if f.path.is_empty() {
                    return Err(format!("Choose a file for {}", p.label));
                }
                let mut reference = Json::object();
                reference.set("file", f.path.as_str());
                if p.kind == "data" && !f.field.is_empty() {
                    reference.set("field", f.field.as_str());
                }
                if p.kind == "data" && f.select_column {
                    reference.set("column", f.column);
                }
                request.set(p.name, reference);
            } else {
                let value = trimmed(&f.value);
                if value.is_empty() {
                    return Err(format!("Provide {}", p.label));
                }
                // Plain text (an element symbol) needs no JSON quotes; text
                // starting with a quote, { or [ is read as JSON.
                if p.kind == "string" && !value.starts_with(['"', '{', '[']) {
                    request.set(p.name, value.as_str());
                } else {
                    request.set(p.name, Json::parse(&value).map_err(|e| format!("{}: {e}", p.label))?);
                }
            }
        }
        Ok(request)
    }

    fn start_calculation(&mut self, structure: &Structure) {
        let started = (|| -> Result<Task, String> {
            let request = self.request(structure)?;
            let directory = new_run_directory()?;
            let id = self.tool;
            let slot: Task = Arc::new(Mutex::new(None));
            let writer = slot.clone();
            std::thread::spawn(move || {
                let result = execute(id, &request, &directory);
                if let Ok(mut s) = writer.lock() {
                    *s = Some(result);
                }
            });
            Ok(slot)
        })();
        match started {
            Ok(task) => {
                self.error.clear();
                self.result = ToolResult::default();
                self.plot = 0;
                self.file = 0;
                self.task = Some(task);
            }
            Err(e) => self.error = e,
        }
    }

    /// A finished run; atom analyses colour the view straight away.
    fn finish(&mut self, done: Result<ToolResult, String>, structure: &Structure) -> Option<AtomProperty> {
        match done {
            Ok(result) => {
                self.error.clear();
                self.result = result;
                let atoms = self.def().is_some_and(|d| d.view == ToolView::Atoms);
                if atoms && self.auto_colour && !self.result.properties.is_empty() {
                    self.property = self.property.min(self.result.properties.len() - 1);
                    let p = &self.result.properties[self.property];
                    if p.values.len() == structure.atoms.len() {
                        return Some(p.clone());
                    }
                }
                None
            }
            Err(e) => {
                self.error = e;
                None
            }
        }
    }

    /// Runs the tool on this thread (UI scripts); the result shows as after a run.
    pub fn run_blocking(&mut self, structure: &Structure) -> Result<(), String> {
        if self.step.is_some() {
            return Err("The tool is editing a pipeline step; use update-step".into());
        }
        let request = self.request(structure);
        let done = request.and_then(|r| execute(self.tool, &r, &new_run_directory()?));
        let error = done.as_ref().err().cloned();
        self.plot = 0;
        self.file = 0;
        if let Some(p) = self.finish(done, structure) {
            self.pending.colour = Some(p);
        }
        error.map_or(Ok(()), Err)
    }

    // -----------------------------------------------------------------------
    // Pipeline steps
    // -----------------------------------------------------------------------

    /// Opens the window on a relax, dynamics or standardize-cell pipeline step,
    /// filled from the step's parameters (missing ones: the catalog defaults).
    pub fn edit_step(&mut self, edit: StepEdit) -> bool {
        let Some(id) = step_tool(&edit.step) else { return false };
        if self.running() || find_tool(id).is_none() {
            return false;
        }
        self.select_tool(id);
        let def = self.def().expect("the step's tool");
        for (i, p) in def.parameters.iter().enumerate() {
            let Some(value) = edit.parameters.find(p.name) else { continue };
            if p.kind == "structure" || value.is_null() {
                continue;
            }
            let f = &mut self.fields[i];
            if scalar_kind(p.kind) {
                match value {
                    Json::Bool(b) => f.value = (if *b { "true" } else { "false" }).into(),
                    Json::Number(n) => f.value = format_g(*n, 15),
                    Json::String(s) => f.value = trimmed(s),
                    _ => {}
                }
                f.enabled = true;
                continue;
            }
            let text = trimmed(&match value {
                Json::String(s) => s.clone(),
                other => other.dump(-1),
            });
            if text.is_empty() {
                // Not given: optional inputs are switched off, required ones keep their default.
                if optional(p) {
                    f.enabled = false;
                }
                continue;
            }
            f.enabled = true;
            let parsed = Json::parse(&text).ok();
            if p.kind == "calculator" {
                let Some(parsed) = parsed.filter(|j| j.is_object()) else { continue };
                let potential = parsed.find("potential").and_then(|v| v.string().ok()).unwrap_or("EMT").to_string();
                let number = |key: &str, fallback: f64| parsed.find(key).and_then(|v| v.number().ok()).unwrap_or(fallback);
                f.potential_file = parsed.find("file").and_then(|v| v.string().ok()).unwrap_or("").to_string();
                match potential.as_str() {
                    "LennardJones" => {
                        f.potential = LENNARD_JONES;
                        f.epsilon = number("epsilon", f.epsilon);
                        f.sigma = number("sigma", f.sigma);
                        f.cutoff = number("cutoff", f.cutoff);
                    }
                    "EAM" => {
                        f.potential = EAM;
                        let format = parsed.find("format").and_then(|v| v.string().ok()).unwrap_or("");
                        f.eam_format = EAM_FORMATS.iter().position(|x| *x == format).unwrap_or(0);
                    }
                    "Tersoff" => f.potential = TERSOFF,
                    "StillingerWeber" => f.potential = STILLINGER_WEBER,
                    "Buckingham" => {
                        f.potential = BUCKINGHAM;
                        let options = Json::Object(parsed.members().map(|m| m.iter().filter(|(k, _)| k != "potential").cloned().collect()).unwrap_or_default());
                        f.potential_options = options.dump(1);
                    }
                    _ => f.potential = EMT,
                }
            } else if let Some(reference) = parsed.as_ref().filter(|j| (p.kind == "data" || p.kind == "file") && j.is_object() && j.contains("file")) {
                // A file reference: {"file": path, "field": ..., "column": ...}.
                f.use_file = true;
                f.path = reference.find("file").and_then(|v| v.string().ok()).unwrap_or("").to_string();
                f.field = reference.find("field").and_then(|v| v.string().ok()).unwrap_or("").to_string();
                let column = reference.find("column").and_then(|v| v.number().ok());
                f.select_column = column.is_some();
                f.column = column.unwrap_or(0.0) as i32;
            } else if let Some(Json::String(path)) = parsed.as_ref().filter(|_| p.kind == "file") {
                f.use_file = true;
                f.path = path.clone();
            } else {
                f.use_file = false;
                f.value = field_text(p.kind, &text);
            }
        }
        let cell = edit.parameters.find("cell").and_then(|v| v.string().ok()).unwrap_or("");
        self.step_cell = STEP_CELLS.iter().position(|c| *c == cell).unwrap_or(0);
        self.step = Some(edit);
        self.open = true;
        true
    }

    /// The tool's inputs as the step's parameters: scalars as JSON numbers and
    /// booleans, other values as their JSON text ("" when not given).
    fn step_parameters(&self) -> Result<Json, String> {
        let step = self.step.as_ref().ok_or("No pipeline step is being edited")?;
        let mut parameters = if step.parameters.is_object() { step.parameters.clone() } else { Json::object() };
        let Some(def) = self.def() else { return Ok(parameters) };
        for (p, f) in def.parameters.iter().zip(&self.fields) {
            if p.kind == "structure" {
                continue; // supplied by the pipeline
            }
            let value = trimmed(&f.value);
            if p.kind == "bool" {
                parameters.set(p.name, value == "true");
                continue;
            }
            if p.kind == "float" || p.kind == "int" {
                if !f.enabled {
                    continue;
                }
                let number = value.parse::<f64>().map_err(|_| format!("Provide a number for {}", p.label))?;
                parameters.set(p.name, if p.kind == "int" { number.round() } else { number });
                continue;
            }
            let text = if !f.enabled {
                String::new()
            } else if p.kind == "calculator" {
                potential_json(f.potential, f.epsilon, f.sigma, f.cutoff, &f.potential_file, f.eam_format, &f.potential_options)?.dump(-1)
            } else if f.use_file && (p.kind == "data" || p.kind == "file") {
                if f.path.is_empty() {
                    return Err(format!("Choose a file for {}", p.label));
                }
                let mut reference = Json::object();
                reference.set("file", f.path.as_str());
                if p.kind == "data" && !f.field.is_empty() {
                    reference.set("field", f.field.as_str());
                }
                if p.kind == "data" && f.select_column {
                    reference.set("column", f.column);
                }
                reference.dump(-1)
            } else if !value.is_empty() && p.kind == "string" && !value.starts_with(['"', '{', '[']) {
                Json::from(value.as_str()).dump(-1) // plain text, accepted by the tool without quotes
            } else {
                value
            };
            parameters.set(p.name, text);
        }
        if step.step == "standardize-cell" {
            parameters.set("cell", STEP_CELLS[self.step_cell.min(2)]);
        }
        Ok(parameters)
    }

    /// Update Step: the step with the inputs as parameters; ends the edit.
    fn update_step(&mut self) -> Result<(), String> {
        let parameters = self.step_parameters()?;
        let mut edit = self.step.take().ok_or("No pipeline step is being edited")?;
        edit.parameters = parameters;
        self.updated = Some(edit);
        self.error.clear();
        self.open = false;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Projects
    // -----------------------------------------------------------------------

    /// The selected tool, inputs and last result as JSON text (for project files).
    pub fn snapshot(&self) -> String {
        let mut state = Json::object();
        state.set("open", self.open);
        if self.def().is_none() {
            return state.dump(-1);
        }
        state.set("tool", self.tool);
        let fields: Vec<Json> = self
            .fields
            .iter()
            .map(|f| {
                let mut item = Json::object();
                item.set("value", f.value.as_str());
                item.set("enabled", f.enabled);
                item.set("path", f.path.as_str());
                item.set("use_file", f.use_file);
                item.set("select_column", f.select_column);
                item.set("column", f.column);
                item.set("field", f.field.as_str());
                item.set("potential", f.potential);
                item.set("epsilon", f.epsilon);
                item.set("sigma", f.sigma);
                item.set("cutoff", f.cutoff);
                item.set("potential_file", f.potential_file.as_str());
                item.set("eam_format", f.eam_format);
                item.set("potential_options", f.potential_options.as_str());
                item
            })
            .collect();
        state.set("fields", Json::Array(fields));
        // Results live in temporary files; they are embedded so the project is self-contained.
        if !self.result.output.as_os_str().is_empty() {
            let embedded = (|| -> Result<(Json, Json), String> {
                let document = Json::parse(&read_text(&self.result.output)?)?;
                let frames = if self.result.structures.as_os_str().is_empty() { Json::Null } else { Json::from(read_text(&self.result.structures)?) };
                Ok((document, frames))
            })();
            if let Ok((document, frames)) = embedded {
                state.set("result", document);
                state.set("frames", frames);
                state.set("property", self.property);
            }
        }
        state.dump(-1)
    }

    /// The trajectory loaded in Trajectory Playback (setTrajectorySource), for
    /// "Use Loaded Trajectory".
    pub fn set_trajectory_source(&mut self, path: &str) {
        if self.trajectory != path {
            self.trajectory = path.to_string();
        }
    }

    /// Restores a snapshot(); empty text keeps the current state.
    pub fn restore(&mut self, text: &str) {
        if text.is_empty() || self.running() {
            return;
        }
        self.step = None;
        let restored = (|| -> Result<(), String> {
            let state = Json::parse(text)?;
            let Some(tool) = state.find("tool").and_then(|t| t.string().ok()).and_then(find_tool) else {
                self.open = false;
                self.tool = "";
                self.result = ToolResult::default();
                return Ok(());
            };
            self.select_tool(tool.id);
            let flag = |item: &Json, key: &str, fallback: bool| item.find(key).and_then(|v| v.boolean().ok()).unwrap_or(fallback);
            let number = |item: &Json, key: &str, fallback: f64| item.find(key).and_then(|v| v.number().ok()).unwrap_or(fallback);
            let string = |item: &Json, key: &str| item.find(key).and_then(|v| v.string().ok()).map(str::to_string);
            // Inputs are restored only when the tool still has the same parameters.
            if let Some(Json::Array(items)) = state.find("fields") {
                if items.len() == self.fields.len() {
                    for (f, item) in self.fields.iter_mut().zip(items) {
                        if !item.is_object() {
                            continue;
                        }
                        if let Some(v) = string(item, "value") {
                            f.value = v;
                        }
                        f.enabled = flag(item, "enabled", f.enabled);
                        f.path = string(item, "path").unwrap_or_default();
                        f.use_file = flag(item, "use_file", f.use_file);
                        f.select_column = flag(item, "select_column", false);
                        f.column = number(item, "column", 0.0) as i32;
                        if let Some(v) = string(item, "field") {
                            f.field = v;
                        }
                        f.potential = (number(item, "potential", 0.0) as usize).min(5);
                        f.epsilon = number(item, "epsilon", f.epsilon);
                        f.sigma = number(item, "sigma", f.sigma);
                        f.cutoff = number(item, "cutoff", f.cutoff);
                        f.potential_file = string(item, "potential_file").unwrap_or_default();
                        f.eam_format = (number(item, "eam_format", 0.0) as usize).min(3);
                        if let Some(v) = string(item, "potential_options") {
                            f.potential_options = v;
                        }
                    }
                }
            }
            if let Some(document) = state.find("result").filter(|d| d.is_object() && d.contains("result")) {
                let directory = new_run_directory()?;
                let result_path = directory.join("result.json");
                write_text(&result_path, &(document.dump(2) + "\n"))?;
                let mut structures = PathBuf::new();
                if let Some(frames) = state.find("frames").and_then(|f| f.string().ok()).filter(|f| !f.is_empty()) {
                    structures = directory.join("frames.extxyz");
                    write_text(&structures, frames)?;
                }
                self.result = ToolResult::from_result(self.tool, result_path, structures, document.at("result")?)?;
                self.property = number(&state, "property", 0.0) as usize;
            }
            self.open = flag(&state, "open", false);
            Ok(())
        })();
        if let Err(e) = restored {
            self.error = format!("Saved scientific tool state could not be restored: {e}");
        }
    }

    // -----------------------------------------------------------------------
    // Drawing
    // -----------------------------------------------------------------------

    fn open_picker(&mut self, target: i32, title: &str, saving: bool, name: &str) {
        self.picker_target = target;
        let path = Path::new(name);
        if let Some(parent) = path.parent().filter(|p| p.is_dir()) {
            self.picker.folder = parent.to_path_buf();
        }
        let file = path.file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default();
        self.picker.show_for(if saving { Purpose::SaveAs } else { Purpose::Open }, title, &file, &[]);
    }

    fn step_input_line(&self, ui: &mut egui::Ui) {
        let input = self.step.as_ref().map(|s| &s.input);
        match input {
            Some(s) if !s.atoms.is_empty() => {
                ui.weak(format!("The structure entering the step: {} atoms, {}", s.atoms.len(), if s.has_unit_cell { "periodic cell" } else { "no cell (open boundaries)" }));
            }
            _ => {
                ui.weak("The structure entering the step (none yet: earlier steps produce no atoms).");
            }
        }
    }

    fn draw_field(&mut self, ui: &mut egui::Ui, i: usize, structure: &Structure) {
        let Some(def) = self.def() else { return };
        let p = &def.parameters[i];
        let kind = p.kind;
        let label = unit_display(p.label);
        ui.push_id(i, |ui| {
            if kind == "structure" && self.step.is_some() {
                // The pipeline supplies the structure of a step.
                ui.label(&label);
                self.step_input_line(ui);
                ui.add_space(4.0);
                return;
            }
            if optional(p) {
                crate::dialog_layout::labelled_check(ui, &mut self.fields[i].enabled, &&label);
            } else if kind != "bool" {
                ui.label(&label);
            }
            if optional(p) && !self.fields[i].enabled {
                ui.add_space(4.0);
                return;
            }
            let enabled = self.fields[i].enabled;
            ui.add_enabled_ui(enabled, |ui| self.field_widget(ui, i, p, &label, structure));
            ui.add_space(4.0);
        });
    }

    fn field_widget(&mut self, ui: &mut egui::Ui, i: usize, p: &'static ParameterDef, label: &str, structure: &Structure) {
        let kind = p.kind;
        let name = p.name;
        match kind {
            "bool" => {
                let mut value = self.fields[i].value.trim() == "true";
                if crate::dialog_layout::labelled_check(ui, &mut value, &label).changed() {
                    self.fields[i].value = (if value { "true" } else { "false" }).into();
                }
            }
            "axes" => {
                let text = self.fields[i].value.clone();
                let inner = text.trim().trim_start_matches('[').trim_end_matches(']');
                let mut values: Vec<bool> = inner.split(',').map(|t| t.contains("true")).collect();
                values.resize(3, false);
                let mut changed = false;
                ui.horizontal(|ui| {
                    for (axis, l) in ["a", "b", "c"].iter().enumerate() {
                        changed |= crate::dialog_layout::labelled_check(ui, &mut values[axis], &*l).changed();
                    }
                });
                if changed {
                    self.fields[i].value = format!("[{},{},{}]", values[0], values[1], values[2]);
                }
            }
            "range" | "vector" => {
                let count = if kind == "range" { 2 } else { 3 };
                let text: String = self.fields[i].value.chars().map(|c| if matches!(c, '[' | ']' | ',') { ' ' } else { c }).collect();
                let mut values = [0.0f64; 3];
                for (k, w) in text.split_whitespace().take(count).enumerate() {
                    match w.parse() {
                        Ok(v) => values[k] = v,
                        Err(_) => break,
                    }
                }
                let mut changed = false;
                ui.horizontal(|ui| {
                    let gap = ui.spacing().item_spacing.x;
                    let each = ((ui.available_width() - gap * (count - 1) as f32) / count as f32).max(40.0);
                    for (k, v) in values.iter_mut().take(count).enumerate() {
                        changed |= component_field(ui, (name, k), v, each);
                    }
                });
                if changed {
                    self.fields[i].value = format!("[{}]", values[..count].iter().map(|&v| format_g(v, 17)).collect::<Vec<_>>().join(","));
                }
                ui.weak(if kind == "range" { "Start / End" } else { "kx / ky / kz" });
            }
            "window" => {
                let mut selected = if self.fields[i].value.trim() == "\"hann\"" { 0 } else { 1 };
                let before = selected;
                egui::ComboBox::from_id_salt("window").selected_text(["Hann", "Rectangular"][selected]).width(ui.available_width().min(340.0) - 10.0).show_ui(ui, |ui| {
                    ui.selectable_value(&mut selected, 0, "Hann");
                    ui.selectable_value(&mut selected, 1, "Rectangular");
                });
                if selected != before {
                    self.fields[i].value = (if selected == 0 { "\"hann\"" } else { "\"none\"" }).into();
                }
            }
            "float" | "int" => {
                let width = ui.available_width().min(340.0);
                ui.add_sized(egui::vec2(width, crate::theme::FRAME_HEIGHT), egui::TextEdit::singleline(&mut self.fields[i].value).margin(egui::vec2(8.0, 4.0)).align(egui::Align2::LEFT_CENTER));
            }
            _ => {
                if matches!(kind, "data" | "structure" | "file") {
                    self.source_widgets(ui, i, p, structure);
                }
                if kind == "calculator" {
                    self.potential_widgets(ui, i);
                } else if !self.fields[i].use_file || kind == "string" {
                    ui.add(egui::TextEdit::singleline(&mut self.fields[i].value).desired_width(f32::INFINITY));
                    if kind == "data" {
                        ui.weak("Inline JSON array, or browse numeric CSV/NPY/JSON/structure data.");
                    }
                }
                let _ = name;
            }
        }
    }

    fn source_widgets(&mut self, ui: &mut egui::Ui, i: usize, p: &'static ParameterDef, structure: &Structure) {
        let kind = p.kind;
        let name = p.name;
        if kind == "structure" && !structure.atoms.is_empty() {
            // Where the structure comes from: the view's, a file or typed values.
            let f = &self.fields[i];
            let mut source = if !f.use_file {
                2
            } else if f.active_label.is_empty() {
                1
            } else {
                0
            };
            let before = source;
            ui.horizontal(|ui| {
                egui::ComboBox::from_id_salt("source").selected_text(["Active Structure", "File", "Enter Values"][source]).width(ui.available_width().min(190.0) - 10.0).show_ui(ui, |ui| {
                    for (k, l) in ["Active Structure", "File", "Enter Values"].iter().enumerate() {
                        ui.selectable_value(&mut source, k, *l);
                    }
                });
                let f = &self.fields[i];
                if source == before && f.use_file && !f.active_label.is_empty() {
                    // "Active structure (4 atoms)" beside the source list: the count.
                    let label = &f.active_label;
                    let count = match (label.find('('), label.rfind(')')) {
                        (Some(a), Some(b)) if b > a => label[a + 1..b].to_string(),
                        _ => label.clone(),
                    };
                    ui.weak(count);
                }
            });
            if source != before {
                match source {
                    0 => {
                        if let Err(e) = self.use_active_structure(i, structure) {
                            self.error = e;
                        }
                    }
                    1 => {
                        let f = &mut self.fields[i];
                        if !f.active_label.is_empty() {
                            f.path.clear(); // the saved copy of the view
                        }
                        f.active_label.clear();
                        f.use_file = true;
                    }
                    _ => self.fields[i].use_file = false,
                }
            }
        } else if kind != "file" {
            let mut source = if self.fields[i].use_file { 0 } else { 1 };
            let before = source;
            egui::ComboBox::from_id_salt("source").selected_text(["File", "Enter Values"][source]).width(ui.available_width().min(170.0) - 10.0).show_ui(ui, |ui| {
                ui.selectable_value(&mut source, 0, "File");
                ui.selectable_value(&mut source, 1, "Enter Values");
            });
            if source != before {
                self.fields[i].use_file = source == 0;
            }
        }
        let f = &self.fields[i];
        if f.use_file && !(kind == "structure" && !f.active_label.is_empty()) {
            let path = f.path.clone();
            if ui.button("Browse").clicked() {
                self.open_picker(i as i32, p.label, false, &path);
            }
            let f = &mut self.fields[i];
            if !f.active_label.is_empty() {
                ui.colored_label(egui::Color32::from_rgb(90, 150, 230), &f.active_label);
            } else {
                ui.label(if f.path.is_empty() { "Choose a file" } else { f.path.as_str() });
            }
            if kind == "data" {
                crate::dialog_layout::collapsing(ui, "file options", &"File options", false, |ui| {
                    ui.add(egui::TextEdit::singleline(&mut f.field).hint_text("JSON field, e.g. result.msd_A2 (optional)").desired_width(f32::INFINITY));
                    crate::dialog_layout::labelled_check(ui, &mut f.select_column, &"Select numeric column");
                    if f.select_column {
                        ui.add(egui::DragValue::new(&mut f.column).range(0..=1_000_000));
                        ui.weak("Column index starts at zero");
                    }
                });
            }
        }
        let cell_input = matches!(name, "cell" | "reference_cell" | "current_cell");
        if kind == "data" && cell_input && structure.has_unit_cell && ui.button("Use Active Cell").clicked() {
            self.fields[i].value = cell_json(structure);
            self.fields[i].use_file = false;
        }
        if !structure.atoms.is_empty() && kind != "structure" && companion_cell(name).is_some() && ui.button("Use Active Structure").clicked() {
            if let Err(e) = self.use_active_structure(i, structure) {
                self.error = e;
            }
        }
        let values_input = matches!(name, "values" | "values_a" | "values_b" | "values_x" | "values_y");
        if kind == "data" && values_input && !structure.atoms.is_empty() && structure.atom_property.len() == structure.atoms.len() && ui.button("Use Active Property").clicked() {
            if let Err(e) = self.use_active_property(i, structure) {
                self.error = e;
            }
        }
        if matches!(kind, "data" | "file") && !self.trajectory.is_empty() && accepts_trajectory(self.tool, name) {
            if ui.button("Use Loaded Trajectory").clicked() {
                let f = &mut self.fields[i];
                f.path = self.trajectory.clone();
                f.use_file = true;
                f.field.clear();
                f.select_column = false;
            }
            ui.weak(format!("Loaded: {}", Path::new(&self.trajectory).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default()));
        }
    }

    fn potential_widgets(&mut self, ui: &mut egui::Ui, i: usize) {
        let mut potential = self.fields[i].potential.min(5);
        egui::ComboBox::from_id_salt("potential").selected_text(POTENTIALS[potential]).width(ui.available_width().min(340.0) - 10.0).show_ui(ui, |ui| {
            for (k, l) in POTENTIALS.iter().enumerate() {
                ui.selectable_value(&mut potential, k, *l);
            }
        });
        self.fields[i].potential = potential;
        match potential {
            TERSOFF | STILLINGER_WEBER => {
                let tersoff = potential == TERSOFF;
                let file = self.fields[i].potential_file.clone();
                ui.horizontal(|ui| {
                    if ui.button("Browse Parameter File").clicked() {
                        self.open_picker(-100 - i as i32, if tersoff { "Tersoff parameter file" } else { "Stillinger-Weber parameter file" }, false, &file);
                    }
                    if !file.is_empty() && ui.button("Use Built-in Si").clicked() {
                        self.fields[i].potential_file.clear();
                    }
                });
                let f = &self.fields[i];
                ui.label(if f.potential_file.is_empty() {
                    if tersoff {
                        "Built-in Si (Tersoff 1988). Choose a LAMMPS .tersoff file for other elements, e.g. SiC.tersoff."
                    } else {
                        "Built-in Si (Stillinger-Weber 1985). Choose a LAMMPS .sw file for other elements."
                    }
                } else {
                    f.potential_file.as_str()
                });
            }
            BUCKINGHAM => {
                ui.label("Pairs A exp(-r/rho) - C/r^6 and fixed charges (e); Coulomb terms are Ewald-summed in periodic cells.");
                ui.add(egui::TextEdit::multiline(&mut self.fields[i].potential_options).desired_rows(6).desired_width(f32::INFINITY).code_editor());
            }
            EAM => {
                let file = self.fields[i].potential_file.clone();
                if ui.button("Browse Potential").clicked() {
                    self.open_picker(-100 - i as i32, "EAM potential file", false, &file);
                }
                let f = &mut self.fields[i];
                ui.label(if f.potential_file.is_empty() { "Choose a LAMMPS-format .eam.alloy, .eam.fs or .eam file" } else { f.potential_file.as_str() });
                egui::ComboBox::from_label("Table format").selected_text(EAM_FORMAT_LABELS[f.eam_format.min(3)]).width(240.0).show_ui(ui, |ui| {
                    for (k, l) in EAM_FORMAT_LABELS.iter().enumerate() {
                        ui.selectable_value(&mut f.eam_format, k, *l);
                    }
                });
                ui.weak("Elements are matched by symbol; the file's cutoff is used.");
            }
            EMT => {
                ui.label("Built-in EMT for Al, Cu, Ag, Au, Ni, Pd and Pt (H, C, N and O are illustrative only).");
            }
            _ => {
                let f = &mut self.fields[i];
                for (value, l) in [(&mut f.epsilon, "Well depth epsilon (eV)"), (&mut f.sigma, "Zero-crossing sigma (Å)"), (&mut f.cutoff, "Cutoff (Å)")] {
                    ui.horizontal(|ui| {
                        ui.add(egui::DragValue::new(value).speed(0.001).max_decimals(6));
                        ui.label(l);
                    });
                }
                ui.weak("Energy is shifted to zero at the cutoff.");
            }
        }
    }

    fn draw_fields(&mut self, ui: &mut egui::Ui, structure: &Structure, group: i32, skip: Option<usize>) {
        for i in self.members(group, skip, true) {
            self.draw_field(ui, i, structure);
        }
    }

    /// A group of fields across columns, short settings on one row.
    fn draw_field_grid(&mut self, ui: &mut egui::Ui, structure: &Structure, group: i32, columns: usize) {
        let members = self.members(group, None, false);
        if members.is_empty() {
            return;
        }
        let mut columns = columns.clamp(1, members.len());
        if ui.available_width() < 300.0 * columns as f32 {
            columns = ((ui.available_width() / 300.0) as usize).max(1);
        }
        let width = ui.available_width() / columns as f32 - 8.0;
        egui::Grid::new(("field grid", group)).num_columns(columns).min_col_width(width).max_col_width(width).show(ui, |ui| {
            for (n, &i) in members.iter().enumerate() {
                ui.vertical(|ui| self.draw_field(ui, i, structure));
                if (n + 1) % columns == 0 {
                    ui.end_row();
                }
            }
        });
    }

    fn draw_status(&mut self, ui: &mut egui::Ui) {
        if self.running() {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label("Calculating...");
                if ui.button("Cancel").on_hover_text("Stops waiting for this run; its result is discarded").clicked() {
                    self.task = None;
                }
            });
        }
        if !self.error.is_empty() {
            ui.colored_label(error_colour(), &self.error);
        }
    }

    fn draw_run_button(&mut self, ui: &mut egui::Ui, label: &str, structure: &Structure, width: f32) {
        let text = if self.step.is_some() { "Update Step" } else { label };
        let width = if width <= 0.0 { ui.available_width() } else { width };
        let running = self.running();
        let pressed = if width >= ui.available_width() - 1.0 {
            crate::dialog_layout::wide_primary(ui, text, 34.0, !running).clicked()
        } else {
            ui.add_enabled_ui(!running, |ui| crate::theme::primary_button(ui, text, Some(egui::vec2(width, 34.0)))).inner.clicked()
        };
        if !pressed {
            return;
        }
        if self.step.is_none() {
            self.start_calculation(structure);
        } else if let Err(e) = self.update_step() {
            // Editing a step: write the inputs into the pipeline instead of running.
            self.error = e;
        }
    }

    fn draw_summary(&self, ui: &mut egui::Ui, id: &str, pairs_per_row: usize) {
        if self.result.summary.is_empty() {
            return;
        }
        let pairs = pairs_per_row.clamp(1, self.result.summary.len());
        let rows = self.result.summary.len().div_ceil(pairs);
        let p = crate::theme::palette(ui.visuals().dark_mode);
        let total = ui.available_width();
        let pair_width = total / pairs as f32;
        let row_h = crate::theme::FRAME_HEIGHT - 8.0;
        let (rect, _) = ui.allocate_exact_size(egui::vec2(total, row_h * rows as f32), egui::Sense::hover());
        let font = egui::TextStyle::Body.resolve(ui.style());
        let alt = if ui.visuals().dark_mode { egui::Color32::from_rgba_unmultiplied(255, 255, 255, 8) } else { egui::Color32::from_rgba_unmultiplied(0, 0, 0, 6) };
        let _ = id;
        for (n, (label, value)) in self.result.summary.iter().enumerate() {
            let (r, c) = (n / pairs, n % pairs);
            let top = rect.top() + r as f32 * row_h;
            let left = rect.left() + c as f32 * pair_width;
            if c == 0 && r % 2 == 1 {
                ui.painter().rect_filled(egui::Rect::from_min_size(egui::pos2(rect.left(), top), egui::vec2(total, row_h)), 0.0, alt);
            }
            let split = left + pair_width * 0.6;
            let label_text = ui.painter().layout(unit_display(label), font.clone(), p.text, pair_width * 0.6 - 8.0);
            ui.painter().galley(egui::pos2(left + 4.0, top + (row_h - label_text.size().y) / 2.0), label_text, p.text);
            ui.painter().text(egui::pos2(split + 4.0, top + row_h / 2.0), egui::Align2::LEFT_CENTER, value, font.clone(), p.text);
            ui.painter().vline(split, top..=top + row_h, egui::Stroke::new(1.0, p.separator));
            if c > 0 {
                ui.painter().vline(left, top..=top + row_h, egui::Stroke::new(1.0, p.separator));
            }
        }
    }

    fn draw_property_controls(&mut self, ui: &mut egui::Ui, structure: &Structure, events: &mut ToolEvents) {
        if self.result.properties.is_empty() {
            return;
        }
        self.property = self.property.min(self.result.properties.len() - 1);
        for (k, p) in self.result.properties.iter().enumerate() {
            let matches = p.values.len() == structure.atoms.len();
            if ui.add_enabled_ui(matches, |ui| crate::dialog_layout::radio(ui, k == self.property, &p.name)).inner.clicked() {
                self.property = k;
                events.colour = Some(p.clone());
            }
        }
        let selected = &self.result.properties[self.property];
        if selected.values.len() != structure.atoms.len() {
            ui.weak(format!("The active structure has {} atoms; this result has {}.", structure.atoms.len(), selected.values.len()));
        } else if ui.button("Colour Atoms in the View").clicked() {
            events.colour = Some(selected.clone());
        }
    }

    /// Plot `index` as shown, with kept runs overlaid.
    fn displayed_plot(&self, index: usize) -> PlotSpec {
        let current = &self.result.plots[index];
        if !self.overlay || self.kept.is_empty() {
            return current.clone();
        }
        let others: Vec<(String, Vec<PlotSpec>)> = self.kept.iter().map(|r| (r.label.clone(), r.plots.clone())).collect();
        overlay_plots(current, "Current", &others)
    }

    fn draw_plots(&mut self, ui: &mut egui::Ui, height: f32) {
        if self.result.plots.is_empty() {
            return;
        }
        self.plot = self.plot.min(self.result.plots.len() - 1);
        if self.result.plots.len() > 1 {
            ui.horizontal_wrapped(|ui| {
                for (k, p) in self.result.plots.iter().enumerate() {
                    let name = if p.title.is_empty() { format!("Plot {}", k + 1) } else { p.title.clone() };
                    ui.selectable_value(&mut self.plot, k, name);
                }
            });
        }
        let spec = self.displayed_plot(self.plot);
        crate::plot_widget::science_plot(ui, &spec, (height - 40.0).max(200.0));
        ui.horizontal(|ui| {
            if ui.button("Export Plot Data").clicked() {
                self.open_picker(-10 - self.plot as i32, "Save plot data", true, "plot.csv");
            }
            if !self.kept.is_empty() {
                let n = self.kept.len();
                crate::dialog_layout::labelled_check(ui, &mut self.overlay, &format!("Overlay {n} kept run{}", if n > 1 { "s" } else { "" }));
            }
        });
    }

    /// Keeps the current result: later results overlay its plots and tabulate its summary.
    fn keep_for_comparison(&mut self) {
        if self.result.output.as_os_str().is_empty() {
            return;
        }
        self.kept.push(KeptRun { label: format!("Run {}", self.kept.len() + 1), plots: self.result.plots.clone(), summary: self.result.summary.clone() });
    }

    fn draw_comparison(&mut self, ui: &mut egui::Ui) {
        if self.kept.is_empty() || self.result.summary.is_empty() {
            return;
        }
        let mut clear = false;
        crate::dialog_layout::collapsing(ui, "comparison", &format!("Compare with kept runs ({})", self.kept.len()), true, |ui| {
            egui::ScrollArea::horizontal().id_salt("comparison scroll").show(ui, |ui| {
                egui::Grid::new("run comparison").striped(true).num_columns(self.kept.len() + 2).show(ui, |ui| {
                    ui.strong("Quantity");
                    for run in &self.kept {
                        ui.strong(&run.label);
                    }
                    ui.strong("Current");
                    ui.end_row();
                    for (label, value) in &self.result.summary {
                        ui.label(unit_display(label));
                        for run in &self.kept {
                            ui.label(run.summary.iter().find(|r| &r.0 == label).map(|r| r.1.as_str()).unwrap_or("-"));
                        }
                        ui.label(value);
                        ui.end_row();
                    }
                });
            });
            clear = ui.button("Clear Kept Runs").clicked();
        });
        if clear {
            self.kept.clear();
        }
    }

    /// The last structure the run produced.
    fn final_structure(&self) -> Result<Structure, String> {
        let reader = |path: &Path| af_io::loader::load_structure_from_file(&path.to_string_lossy());
        let mut frames = read_frames(&self.result.structures, Some(&reader))?;
        frames.pop().map(|f| f.structure).ok_or_else(|| "The run produced no structures".to_string())
    }

    fn draw_save_buttons(&mut self, ui: &mut egui::Ui, events: &mut ToolEvents) {
        self.draw_comparison(ui);
        ui.horizontal_wrapped(|ui| {
            if ui.button("Keep for Comparison").on_hover_text("Keep this result; change inputs and run again to overlay plots and compare values.").clicked() {
                self.keep_for_comparison();
            }
            if ui.button("Save Results").clicked() {
                self.open_picker(-3, "Save scientific results", true, "analysis.json");
            }
            if !self.result.structures.as_os_str().is_empty() {
                if ui.button("Save Structures").clicked() {
                    self.open_picker(-4, "Save result structures or trajectory", true, "frames.extxyz");
                }
                if ui.button("Open Final Structure in a New Tab").clicked() {
                    match self.final_structure() {
                        Ok(s) => events.opened.push((s, format!("{} result", title_case(self.def().map(|d| d.title).unwrap_or("Tool"))))),
                        Err(e) => self.error = e,
                    }
                }
            }
        });
    }

    fn draw_details(&self, ui: &mut egui::Ui) {
        if self.result.report.is_empty() {
            return;
        }
        crate::dialog_layout::collapsing(ui, "full report", &"Full report", false, |ui| {
            ui.label(egui::RichText::new(&self.result.report).monospace());
        });
    }

    /// The first OBJ mesh of the result (surface.obj, regions.obj).
    fn result_mesh(&self) -> Option<&(String, String)> {
        self.result.files.iter().find(|(name, _)| name.len() > 4 && name.ends_with(".obj"))
    }

    fn set_mesh_shown(&mut self, show: bool) -> Result<(), String> {
        let Some((_, text)) = self.result_mesh() else { return Err("The result has no mesh".into()) };
        if show && self.mesh_from != self.result.output {
            self.mesh.vertices = triangles_from_obj(text)?;
            self.mesh_from = self.result.output.clone();
        }
        self.mesh.visible = show;
        Ok(())
    }

    /// A mesh of the result (surface, regions), shown translucent over the structure.
    fn draw_mesh_controls(&mut self, ui: &mut egui::Ui) {
        let Some((name, _)) = self.result_mesh() else { return };
        let label = if name == "surface.obj" { "Show the surface in the view" } else { "Show the mesh in the view" };
        section(ui, "Surface");
        let mut show = self.mesh.visible && self.mesh_from == self.result.output;
        if crate::dialog_layout::labelled_check(ui, &mut show, &label).changed() {
            if let Err(e) = self.set_mesh_shown(show) {
                self.error = e;
                self.mesh.visible = false;
            }
        }
        if show {
            ui.horizontal(|ui| {
                crate::dialog_layout::slider_f32(ui, &mut self.mesh.opacity, 0.05..=1.0, ui.available_width().min(240.0), |v| format!("{v:.2}"));
                ui.label("Opacity");
            });
            ui.horizontal(|ui| {
                ui.color_edit_button_rgb(&mut self.mesh.colour);
                ui.label("Colour");
            });
            ui.weak(format!("{} triangles", self.mesh.vertices.len() / 3));
        }
    }

    fn empty_state(&self, ui: &mut egui::Ui, waiting: &str) {
        if !self.error.is_empty() && !self.running() {
            crate::dialog_layout::note(ui, crate::dialog_layout::Note::Error, &self.error);
        } else {
            crate::dialog_layout::empty_state(ui, if self.running() { "Calculating..." } else { waiting });
        }
    }

    fn has_result(&self) -> bool {
        !self.result.output.as_os_str().is_empty()
    }

    /// Per-atom analyses: one panel beside the view, fed by the active structure.
    fn atom_layout(&mut self, ui: &mut egui::Ui, structure: &Structure, events: &mut ToolEvents) {
        let input = self.active_input();
        let busy = self.running();
        ui.add_enabled_ui(!busy, |ui| {
            section(ui, "Atoms");
            if let Some(input) = input {
                if self.step.is_some() {
                    self.step_input_line(ui);
                } else {
                    crate::dialog_layout::labelled_check(ui, &mut self.use_active, &"Analyse the active structure");
                    if self.use_active {
                        if structure.atoms.is_empty() {
                            ui.weak("No structure is loaded.");
                        } else {
                            ui.weak(format!("{} atoms, {}", structure.atoms.len(), if structure.has_unit_cell { "periodic cell" } else { "no cell (open boundaries)" }));
                        }
                    } else {
                        self.draw_field(ui, input, structure);
                    }
                }
            }
            let skip = if self.use_active || self.step.is_some() { input } else { None };
            if self.has_fields(0, skip) {
                self.draw_fields(ui, structure, 0, skip);
            }
            if self.has_fields(1, skip) {
                section(ui, "Settings");
                self.draw_fields(ui, structure, 1, skip);
            }
            if self.has_fields(2, skip) {
                crate::dialog_layout::collapsing(ui, "more options", &"More options", false, |ui| self.draw_fields(ui, structure, 2, skip));
            }
        });
        ui.add_space(4.0);
        self.draw_run_button(ui, "Analyse", structure, -1.0);
        self.draw_status(ui);
        if self.has_result() {
            section(ui, "Result");
            self.draw_summary(ui, "atom summary", 1);
            if !self.result.properties.is_empty() {
                section(ui, "Colour atoms by");
                self.draw_property_controls(ui, structure, events);
                crate::dialog_layout::labelled_check(ui, &mut self.auto_colour, &"Colour automatically after each analysis");
            }
            self.draw_mesh_controls(ui);
            if !self.result.plots.is_empty() {
                crate::dialog_layout::collapsing(ui, "plots", &"Plots", true, |ui| self.draw_plots(ui, 260.0));
            }
            self.draw_save_buttons(ui, events);
            self.draw_details(ui);
        }
    }

    fn settings_column(&mut self, ui: &mut egui::Ui, structure: &Structure, groups: &[(i32, &str, bool)], run: &str) {
        ui.add_enabled_ui(!self.running(), |ui| {
            for &(group, title, collapsed) in groups {
                if !self.has_fields(group, None) {
                    continue;
                }
                if collapsed {
                    crate::dialog_layout::collapsing(ui, title, &title, false, |ui| self.draw_fields(ui, structure, group, None));
                } else {
                    section(ui, title);
                    self.draw_fields(ui, structure, group, None);
                }
            }
        });
        ui.add_space(4.0);
        self.draw_run_button(ui, run, structure, -1.0);
        self.draw_status(ui);
    }

    /// Two panes side by side, or stacked when the window is narrow.
    fn panes(&mut self, ui: &mut egui::Ui, stack_below: f32, left_width: f32, left: impl FnOnce(&mut Self, &mut egui::Ui), right: impl FnOnce(&mut Self, &mut egui::Ui)) {
        let height = ui.available_height().max(300.0);
        if ui.available_width() < stack_below {
            egui::ScrollArea::vertical().id_salt("stacked").show(ui, |ui| {
                crate::dialog_layout::child_frame(ui, |ui| left(self, ui));
                right(self, ui);
            });
            return;
        }
        let pal = crate::theme::palette(ui.visuals().dark_mode);
        ui.horizontal_top(|ui| {
            ui.allocate_ui_with_layout(egui::vec2(left_width, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                egui::Frame::new().fill(pal.surface).stroke(egui::Stroke::new(1.0, pal.border)).corner_radius(8).inner_margin(egui::Margin::symmetric(10, 8)).show(ui, |ui| {
                    ui.set_width(left_width - 22.0);
                    ui.set_height(height - 18.0);
                    egui::ScrollArea::vertical().id_salt("left pane").auto_shrink([false, false]).show(ui, |ui| left(self, ui));
                });
            });
            let rest = ui.available_width();
            ui.allocate_ui_with_layout(egui::vec2(rest, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                egui::Frame::new().fill(pal.surface).corner_radius(8).inner_margin(egui::Margin::symmetric(10, 8)).show(ui, |ui| {
                    ui.set_width(rest - 20.0);
                    ui.set_height(height - 16.0);
                    egui::ScrollArea::vertical().id_salt("right pane").auto_shrink([false, false]).show(ui, |ui| right(self, ui));
                });
            });
        });
    }

    /// Curves and spectra: settings sidebar and a large plot area.
    fn plot_layout(&mut self, ui: &mut egui::Ui, structure: &Structure, events: &mut ToolEvents) {
        let sidebar = 380.0f32.min(ui.available_width() * 0.42);
        self.panes(
            ui,
            900.0,
            sidebar,
            |s, ui| s.settings_column(ui, structure, &[(0, "Data", false), (1, "Settings", false), (2, "Optional inputs", true)], "Compute"),
            |s, ui| {
                if !s.has_result() {
                    s.empty_state(ui, "Results appear here after you run the tool.");
                    return;
                }
                let height = if s.result.plots.is_empty() { 0.0 } else { 240.0f32.max(ui.available_height() * 0.62) };
                s.draw_plots(ui, height);
                ui.add_space(4.0);
                s.draw_summary(ui, "plot summary", 1);
                s.draw_save_buttons(ui, events);
                s.draw_details(ui);
            },
        );
    }

    /// Trajectory analyses: the source across the top, settings in a row, a wide plot.
    fn trajectory_layout(&mut self, ui: &mut egui::Ui, structure: &Structure, events: &mut ToolEvents) {
        egui::ScrollArea::vertical().id_salt("trajectory").show(ui, |ui| {
            ui.add_enabled_ui(!self.running(), |ui| {
                crate::dialog_layout::child_frame(ui, |ui| {
                    section(ui, "Trajectory");
                    if self.trajectory.is_empty() {
                        ui.weak("Load frames in Trajectory playback to use them here, or choose files below.");
                    } else {
                        ui.label(format!("Playback trajectory: {}", Path::new(&self.trajectory).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default()));
                    }
                    self.draw_field_grid(ui, structure, 0, 2);
                });
                if self.has_fields(1, None) {
                    section(ui, "Settings");
                    self.draw_field_grid(ui, structure, 1, 3);
                }
                if self.has_fields(2, None) {
                    crate::dialog_layout::collapsing(ui, "optional", &"Optional inputs", false, |ui| self.draw_field_grid(ui, structure, 2, 3));
                }
            });
            self.draw_run_button(ui, "Analyse Trajectory", structure, 220.0);
            self.draw_status(ui);
            if !self.has_result() {
                return;
            }
            ui.separator();
            self.draw_plots(ui, 320.0);
            self.draw_summary(ui, "trajectory summary", 1);
            self.draw_save_buttons(ui, events);
            self.draw_details(ui);
        });
    }

    /// Material properties: inputs on one side, the computed quantities as large values and tensors.
    fn properties_layout(&mut self, ui: &mut egui::Ui, structure: &Structure, events: &mut ToolEvents) {
        let form = ui.available_width() * 0.42;
        self.panes(
            ui,
            860.0,
            form,
            |s, ui| s.settings_column(ui, structure, &[(0, "Data", false), (3, "Interatomic potential", false), (1, "Settings", false), (2, "Optional inputs", true)], "Calculate"),
            |s, ui| {
                if !s.has_result() {
                    s.empty_state(ui, "Computed properties appear here after you run the tool.");
                    return;
                }
                // Headline quantities as cards: label above a large value.
                let columns = ((ui.available_width() / 190.0) as usize).max(1);
                if !s.result.summary.is_empty() {
                    egui::Grid::new("property cards").num_columns(columns).min_col_width(170.0).show(ui, |ui| {
                        for (n, (label, value)) in s.result.summary.iter().enumerate() {
                            ui.vertical(|ui| {
                                ui.weak(unit_display(label));
                                ui.label(egui::RichText::new(value).size(20.0));
                            });
                            if (n + 1) % columns == 0 {
                                ui.end_row();
                            }
                        }
                    });
                }
                for (m, (label, rows)) in s.result.matrices.iter().enumerate() {
                    section(ui, &unit_display(label));
                    egui::Grid::new(("matrix", m)).striped(true).num_columns(rows[0].len()).show(ui, |ui| {
                        for row in rows {
                            for v in row {
                                ui.monospace(format_g(*v, 5));
                            }
                            ui.end_row();
                        }
                    });
                }
                if !s.result.plots.is_empty() {
                    crate::dialog_layout::collapsing(ui, "plots", &"Plots", true, |ui| s.draw_plots(ui, 280.0));
                }
                s.draw_save_buttons(ui, events);
                s.draw_details(ui);
            },
        );
    }

    /// Simulations: numbered setup steps, a run bar and a run summary.
    fn simulation_layout(&mut self, ui: &mut egui::Ui, structure: &Structure, events: &mut ToolEvents) {
        let setup_height = if self.has_result() { 220.0 } else { (ui.available_height() - 96.0).max(200.0) };
        crate::dialog_layout::child_frame(ui, |ui| {
            egui::ScrollArea::vertical().id_salt("simulation steps").max_height(setup_height).show(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.add_enabled_ui(!self.running(), |ui| {
                    let mut number = 1;
                    for (group, title) in [(0, "System"), (3, "Interatomic potential"), (1, "Protocol")] {
                        if self.has_fields(group, None) {
                            section(ui, &format!("{number}  {title}"));
                            number += 1;
                            self.draw_fields(ui, structure, group, None);
                        }
                    }
                    if self.has_fields(2, None) {
                        crate::dialog_layout::collapsing(ui, "advanced", &"Advanced options", false, |ui| self.draw_fields(ui, structure, 2, None));
                    }
                });
            });
        });
        let label = if self.running() { "Running..." } else { "Run Simulation" };
        self.draw_run_button(ui, label, structure, 200.0);
        self.draw_status(ui);
        if !self.has_result() {
            return;
        }
        egui::ScrollArea::vertical().id_salt("simulation result").show(ui, |ui| {
            section(ui, "Run summary");
            let wide = ui.available_width() > 760.0;
            if wide && !self.result.plots.is_empty() {
                ui.horizontal_top(|ui| {
                    let w = ui.available_width() * 0.4;
                    ui.allocate_ui_with_layout(egui::vec2(w, 320.0), egui::Layout::top_down(egui::Align::Min), |ui| {
                        egui::ScrollArea::vertical().id_salt("run numbers").max_height(320.0).show(ui, |ui| self.draw_summary(ui, "simulation summary", 1));
                    });
                    ui.vertical(|ui| self.draw_plots(ui, 280.0));
                });
            } else {
                self.draw_summary(ui, "simulation summary", if wide { 2 } else { 1 });
                self.draw_plots(ui, 260.0);
            }
            self.draw_save_buttons(ui, events);
            self.draw_details(ui);
        });
    }

    /// Input generators: a short form and a browser for the generated files.
    fn generator_layout(&mut self, ui: &mut egui::Ui, structure: &Structure) {
        let form_height = 220.0f32.max(ui.available_height() * 0.42);
        crate::dialog_layout::child_frame(ui, |ui| {
            egui::ScrollArea::vertical().id_salt("generator form").max_height(form_height).show(ui, |ui| {
                ui.add_enabled_ui(!self.running(), |ui| {
                    ui.columns(2, |columns| {
                        section(&mut columns[0], "Structure and data");
                        self.draw_fields(&mut columns[0], structure, 0, None);
                        if self.has_fields(3, None) {
                            section(&mut columns[0], "Potential");
                            self.draw_fields(&mut columns[0], structure, 3, None);
                        }
                        section(&mut columns[1], "Options");
                        self.draw_field_grid(&mut columns[1], structure, 1, 2);
                        if self.has_fields(2, None) {
                            crate::dialog_layout::collapsing(&mut columns[1], "optional", &"Optional", false, |ui| self.draw_field_grid(ui, structure, 2, 2));
                        }
                    });
                });
            });
        });
        self.draw_run_button(ui, "Generate Files", structure, 200.0);
        self.draw_status(ui);
        section(ui, "Generated files");
        if self.result.files.is_empty() {
            if !self.error.is_empty() && !self.running() {
                ui.colored_label(error_colour(), &self.error);
            } else {
                ui.weak(if !self.has_result() { "Generate to preview the files here." } else { "This run produced no files." });
            }
            return;
        }
        self.file = self.file.min(self.result.files.len() - 1);
        let height = (ui.available_height() - 40.0).max(200.0);
        let list_width = 220.0f32.min(ui.available_width() * 0.3);
        ui.horizontal_top(|ui| {
            ui.allocate_ui_with_layout(egui::vec2(list_width, height), egui::Layout::top_down(egui::Align::Min), |ui| {
                crate::dialog_layout::child_frame(ui, |ui| {
                    egui::ScrollArea::vertical().id_salt("file list").max_height(height - 14.0).show(ui, |ui| {
                        for (k, (name, text)) in self.result.files.iter().enumerate() {
                            if ui.selectable_label(k == self.file, name).clicked() {
                                self.file = k;
                            }
                            ui.weak(format!("{} lines", text.matches('\n').count()));
                        }
                    });
                });
            });
            let text = &self.result.files[self.file].1;
            egui::ScrollArea::both().id_salt("file preview").max_height(height).show(ui, |ui| {
                ui.add(egui::Label::new(egui::RichText::new(text.as_str()).monospace()).wrap_mode(egui::TextWrapMode::Extend));
            });
        });
        ui.horizontal(|ui| {
            let (name, text) = self.result.files[self.file].clone();
            if ui.button("Save This File").clicked() {
                self.open_picker(-6, "Save generated file", true, &name);
            }
            if ui.button("Save All Files").clicked() {
                let first = self.result.files[0].0.clone();
                self.open_picker(-5, "Choose the folder (save as the first file)", true, &first);
            }
            if ui.button("Copy to Clipboard").clicked() {
                ui.ctx().copy_text(text);
            }
            if ui.button("Save Results").clicked() {
                self.open_picker(-3, "Save scientific results", true, "analysis.json");
            }
        });
    }

    /// Writes what the file picker was opened for.
    fn picked(&mut self, path: &Path) -> Result<(), String> {
        let copy = |from: &Path, to: &Path| std::fs::copy(from, to).map(|_| ()).map_err(|e| format!("Cannot save {}: {e}", to.display()));
        let target = self.picker_target;
        if target == -5 {
            let folder = path.parent().unwrap_or(Path::new(""));
            for (name, text) in &self.result.files {
                write_text(&folder.join(name), text)?;
            }
        } else if target == -6 {
            if let Some((_, text)) = self.result.files.get(self.file) {
                write_text(path, text)?;
            }
        } else if target == -3 {
            copy(&self.result.output, path)?;
        } else if target <= -10 && target > -100 {
            let k = (-10 - target) as usize;
            if k < self.result.plots.len() {
                write_text(path, &plot_csv(&self.displayed_plot(k)))?;
            }
        } else if target == -4 {
            copy(&self.result.structures, path)?;
        } else if target <= -100 {
            if let Some(f) = self.fields.get_mut((-100 - target) as usize) {
                f.potential_file = path.to_string_lossy().into_owned();
            }
        } else if target >= 0 {
            if let Some(f) = self.fields.get_mut(target as usize) {
                f.path = path.to_string_lossy().into_owned();
                f.active_label.clear();
                f.use_file = true;
            }
        }
        Ok(())
    }

    /// Draws the window; returns structures to open and atoms to colour.
    pub fn show(&mut self, ctx: &egui::Context, structure: &Structure) -> ToolEvents {
        let mut events = std::mem::take(&mut self.pending);
        // A tool just opened takes the active structure for its structure
        // inputs that have no file yet (once; a file can still be chosen).
        if self.active_defaults_pending && self.open && self.def().is_some() && !structure.atoms.is_empty() {
            self.active_defaults_pending = false;
            if self.step.is_none() {
                let def = self.def().expect("checked");
                for (i, p) in def.parameters.iter().enumerate() {
                    if p.kind == "structure" && self.fields[i].use_file && self.fields[i].path.is_empty() {
                        let _ = self.use_active_structure(i, structure);
                    }
                }
            }
        }
        if let Some(slot) = &self.task {
            let done = slot.lock().ok().and_then(|mut s| s.take());
            match done {
                Some(done) => {
                    self.task = None;
                    if let Some(p) = self.finish(done, structure) {
                        events.colour = Some(p);
                    }
                }
                None => ctx.request_repaint_after(std::time::Duration::from_millis(100)),
            }
        }
        if !self.open {
            return events;
        }
        let Some(def) = self.def() else {
            self.open = false;
            return events;
        };
        let (window_id, size): (&str, [f32; 2]) = match def.view {
            ToolView::Atoms => ("atom analysis", [430.0, 700.0]),
            ToolView::Plot => ("plot analysis", [1120.0, 700.0]),
            ToolView::Trajectory => ("trajectory analysis", [1000.0, 820.0]),
            ToolView::Properties => ("property calculator", [1040.0, 680.0]),
            ToolView::Simulation => ("simulation run", [900.0, 800.0]),
            ToolView::Generator => ("input generator", [980.0, 760.0]),
        };
        let screen = ctx.content_rect().size();
        let size = egui::vec2(size[0].min(screen.x - 40.0), size[1].min(screen.y - 80.0));
        let mut open = self.open;
        crate::theme::window(title_case(def.title))
            .id(egui::Id::unique(window_id))
            .open(&mut open)
            .default_size(size)
            .min_size([380.0, 320.0])
            .collapsible(false)
            .show(ctx, |ui| {
                if let Some(step) = &self.step {
                    crate::structure_dialog::step_banner(ui, &step.step.clone());
                    if step.step == "standardize-cell" {
                        egui::ComboBox::from_label("Cell to keep").selected_text(["Primitive", "Conventional", "Symmetrized"][self.step_cell.min(2)]).width(220.0).show_ui(ui, |ui| {
                            for (k, l) in ["Primitive", "Conventional", "Symmetrized"].iter().enumerate() {
                                ui.selectable_value(&mut self.step_cell, k, *l);
                            }
                        });
                    }
                }
                ui.weak(title_case(def.category));
                ui.label(def.help.split('\n').next().unwrap_or(""));
                ui.add_space(4.0);
                match def.view {
                    ToolView::Atoms => {
                        egui::ScrollArea::vertical().id_salt("atom layout").auto_shrink([false, false]).show(ui, |ui| self.atom_layout(ui, structure, &mut events));
                    }
                    ToolView::Plot => self.plot_layout(ui, structure, &mut events),
                    ToolView::Trajectory => self.trajectory_layout(ui, structure, &mut events),
                    ToolView::Properties => self.properties_layout(ui, structure, &mut events),
                    ToolView::Simulation => self.simulation_layout(ui, structure, &mut events),
                    ToolView::Generator => self.generator_layout(ui, structure),
                }
            });
        self.open = open && self.open;
        if let Some(path) = self.picker.show(ctx) {
            self.picker.open = false;
            if let Err(e) = self.picked(&path) {
                self.error = e;
            }
        }
        if !self.open {
            self.step = None; // closed with the window's X or Update Step
        }
        events
    }

    /// The result mesh over the view (drawSurfaceMeshOverlay): triangles
    /// shaded by their orientation to the viewer, drawn back to front.
    pub fn draw_overlay(&self, painter: &egui::Painter, rect: egui::Rect, scale: f32, frame: &FrameView) {
        let mesh = &self.mesh;
        if !mesh.visible || mesh.vertices.len() < 3 {
            return;
        }
        let near = 0.03;
        let eye = |p: Vec3| (frame.view * Vec4::new(p.x, p.y, p.z, 1.0)).truncate();
        let mut triangles: Vec<(f32, Vec<egui::Pos2>, f32)> = Vec::with_capacity(mesh.vertices.len() / 3);
        for t in mesh.vertices.chunks_exact(3) {
            let world = [Vec3::from_array(t[0]), Vec3::from_array(t[1]), Vec3::from_array(t[2])];
            let e = world.map(eye);
            let normal = (e[1] - e[0]).cross(e[2] - e[0]);
            let length = normal.length();
            if length <= 0.0 {
                continue;
            }
            // Lambert-like shading from the viewer; both sides lit alike.
            let shade = 0.35 + 0.65 * (normal.z / length).abs();
            let depth = (e[0].z + e[1].z + e[2].z) / 3.0;
            // Clipped at the near plane (parts behind the camera are not drawn).
            let mut clipped = Vec::with_capacity(4);
            for k in 0..3 {
                let (a, b) = (k, (k + 1) % 3);
                let (inside_a, inside_b) = (e[a].z <= -near, e[b].z <= -near);
                if inside_a {
                    clipped.push(world[a]);
                }
                if inside_a != inside_b && (e[b].z - e[a].z).abs() > 1e-6 {
                    let s = ((-near - e[a].z) / (e[b].z - e[a].z)).clamp(0.0, 1.0);
                    clipped.push(world[a] + (world[b] - world[a]) * s);
                }
            }
            if clipped.len() < 3 {
                continue;
            }
            let points: Option<Vec<egui::Pos2>> = clipped.iter().map(|&p| crate::view::project(p, frame).map(|(x, y)| rect.min + egui::vec2(x / scale, y / scale))).collect();
            if let Some(points) = points {
                triangles.push((depth, points, shade));
            }
        }
        // Far (most negative eye z) first.
        triangles.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let alpha = (mesh.opacity.clamp(0.0, 1.0) * 255.0) as u8;
        let channel = |c: f32, shade: f32| (c * shade * 255.0).clamp(0.0, 255.0) as u8;
        // One mesh without anti-aliasing feathering, which would extrude the
        // many thin triangles of a surface into spikes (as ImGui's
        // AddTriangleFilled draws them: plain triangles).
        let mut shape = egui::Mesh::default();
        for (_, points, shade) in triangles {
            let colour = egui::Color32::from_rgba_unmultiplied(channel(mesh.colour[0], shade), channel(mesh.colour[1], shade), channel(mesh.colour[2], shade), alpha);
            let base = shape.vertices.len() as u32;
            for p in &points {
                shape.colored_vertex(*p, colour);
            }
            for k in 1..points.len() as u32 - 1 {
                shape.add_triangle(base, base + k, base + k + 1);
            }
        }
        painter.add(egui::Shape::mesh(shape));
    }

    /// UI-script commands (`toolwin COMMAND ARGS`): surface on|off, expect-mesh
    /// N, keep, expect-kept N, save-results PATH, save-structures PATH,
    /// save-files FOLDER, open-final, expect-summary LABEL = VALUE,
    /// expect-field NAME TEXT, use-active on|off, use-active-structure NAME,
    /// step-cell N, update-step,
    /// expect-step yes|no, expect-error TEXT, restore-snapshot.
    pub fn script(&mut self, command: &str, args: &str, structure: &Structure) -> Result<(), String> {
        let args = args.trim();
        let number = || args.split_whitespace().next().and_then(|w| w.parse::<f64>().ok()).ok_or(format!("{command} N"));
        match command {
            "surface" => self.set_mesh_shown(args == "on")?,
            "expect-mesh" => {
                let n = number()? as usize;
                let shown = if self.mesh.visible { self.mesh.vertices.len() / 3 } else { 0 };
                if shown != n {
                    return Err(format!("the view shows {shown} mesh triangles"));
                }
            }
            "keep" => self.keep_for_comparison(),
            "expect-kept" => {
                if self.kept.len() != number()? as usize {
                    return Err(format!("{} kept runs", self.kept.len()));
                }
            }
            "save-results" => {
                self.picker_target = -3;
                self.picked(Path::new(args))?;
            }
            "save-structures" => {
                if self.result.structures.as_os_str().is_empty() {
                    return Err("the run produced no structures".into());
                }
                self.picker_target = -4;
                self.picked(Path::new(args))?;
            }
            "save-files" => {
                if self.result.files.is_empty() {
                    return Err("the run produced no files".into());
                }
                self.picker_target = -5;
                self.picked(&Path::new(args).join(&self.result.files[0].0))?;
            }
            "open-final" => {
                let s = self.final_structure()?;
                let title = format!("{} result", title_case(self.def().map(|d| d.title).unwrap_or("Tool")));
                self.pending.opened.push((s, title));
            }
            "expect-summary" => {
                let (label, value) = args.split_once('=').ok_or("expect-summary LABEL = VALUE")?;
                let (label, value) = (label.trim(), value.trim());
                match self.result.summary.iter().find(|r| r.0 == label) {
                    Some(r) if r.1 == value => {}
                    Some(r) => return Err(format!("{label} = {}", r.1)),
                    None => return Err(format!("no {label} in {:?}", self.result.summary.iter().map(|r| &r.0).collect::<Vec<_>>())),
                }
            }
            "expect-field" => {
                let (name, text) = args.split_once(' ').unwrap_or((args, ""));
                let def = self.def().ok_or("No tool is open")?;
                let i = def.parameters.iter().position(|p| p.name == name).ok_or(format!("no parameter {name}"))?;
                let f = &self.fields[i];
                let shown = if f.use_file { f.path.clone() } else { f.value.trim().to_string() };
                if shown != text.trim() {
                    return Err(format!("{name} is {shown:?}"));
                }
            }
            "use-active" => self.use_active = args == "on",
            "use-active-structure" => {
                // The input's "Use Active Structure" button (positions, with their cell and axes).
                let def = self.def().ok_or("No tool is open")?;
                let i = def.parameters.iter().position(|p| p.name == args).ok_or(format!("no parameter {args}"))?;
                if structure.atoms.is_empty() {
                    return Err("no structure is open".into());
                }
                self.use_active_structure(i, structure)?;
            }
            "step-cell" => self.step_cell = (number()? as usize).min(2),
            "update-step" => self.update_step()?,
            "expect-step" => {
                if self.step.is_some() != (args == "yes") {
                    return Err(format!("editing a step: {}", self.step.is_some()));
                }
            }
            "expect-error" => {
                if !self.error.contains(args) {
                    return Err(format!("the error is {:?}", self.error));
                }
            }
            "restore-snapshot" => {
                // A project round trip: the state saved and restored.
                let (summary, plots) = (self.result.summary.clone(), self.result.plots.len());
                let text = self.snapshot();
                self.restore(&text);
                if !self.error.is_empty() {
                    return Err(self.error.clone());
                }
                if self.result.summary != summary || self.result.plots.len() != plots {
                    return Err("the restored result differs".into());
                }
            }
            "run-background" => {
                self.start_calculation(structure);
                if !self.error.is_empty() {
                    return Err(self.error.clone());
                }
            }
            other => return Err(format!("unknown toolwin command {other}")),
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn obj_triangles() {
        let quad = "# a quad\nv 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nvn 0 0 1\nf 1//1 2//1 3//1 4//1\n";
        let t = triangles_from_obj(quad).unwrap();
        assert_eq!(t.len(), 6);
        assert_eq!(t[3..6], [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]);
        // Negative indices count back from the vertices read so far.
        let t = triangles_from_obj("v 0 0 0\nv 1 0 0\nv 0 1 0\nf -3 -2 -1\nv 5 5 5\nf 1/1/1 2/2/2 -1\n").unwrap();
        assert_eq!(t.len(), 6);
        assert_eq!(t[5], [5.0, 5.0, 5.0]);
        assert_eq!(triangles_from_obj("v 0 0 0\nf 1 2 3\n").unwrap_err(), "OBJ face refers to a missing vertex");
        assert!(triangles_from_obj("v 0 x 0\n").unwrap_err().starts_with("Malformed OBJ vertex"));
        assert_eq!(triangles_from_obj("v 0 0 0\nv 1 0 0\nf 1 2\n").unwrap().len(), 0, "two-vertex faces make no triangle");
    }

    #[test]
    fn labels() {
        assert_eq!(pretty_label("max_force_eV_per_A"), "Max force (eV/A)");
        assert_eq!(pretty_label("msd_A2"), "Msd (A^2)");
        assert_eq!(pretty_label("energy_eV_per_atom"), "Energy (eV/atom)");
        assert_eq!(pretty_label("atoms_per_cell"), "Atoms per cell");
        assert_eq!(pretty_label("spacing_inv_A"), "Spacing (1/A)");
        assert_eq!(pretty_label("eV"), "EV");
        assert_eq!(title_case("Powder X-ray diffraction"), "Powder X-ray Diffraction");
        assert_eq!(title_case("surface mesh, area and solid volume"), "Surface mesh, Area and Solid Volume");
        assert_eq!(title_case("radial distribution g(r)"), "Radial Distribution g(r)");
    }

    #[test]
    fn summary_and_matrices() {
        let r = Json::parse(r#"{"energy_eV": -3.5123456789, "converged": true, "name": "fcc", "cell": [[1, 1e-12], [0, 2]], "ragged": [[1, 2], [3]], "list": [1, 2]}"#).unwrap();
        assert_eq!(summary_of(&r), vec![("Energy (eV)".to_string(), "-3.51235".to_string()), ("Converged".into(), "yes".into()), ("Name".into(), "fcc".into())]);
        assert_eq!(matrices_of(&r), vec![("Cell".to_string(), vec![vec![1.0, 0.0], vec![0.0, 2.0]])]);
    }

    #[test]
    fn relax_step_round_trip() {
        let mut w = ToolWindow::default();
        let mut parameters = Json::object();
        parameters.set("fmax", 0.05);
        parameters.set("steps", 200.0);
        parameters.set("relax_cell", true);
        parameters.set("calculator", "{\"potential\": \"LennardJones\", \"epsilon\": 0.2, \"sigma\": 2.5, \"cutoff\": 7}");
        assert!(w.edit_step(StepEdit { step: "relax".into(), parameters, ..Default::default() }));
        assert_eq!(w.fields[2].value, "0.05");
        assert_eq!(w.fields[1].potential, LENNARD_JONES);
        w.fields[3].value = "300".into();
        w.update_step().unwrap();
        let out = w.updated.take().unwrap().parameters;
        assert_eq!(out.find("steps").unwrap().number().unwrap(), 300.0);
        assert_eq!(out.find("fmax").unwrap().number().unwrap(), 0.05);
        assert_eq!(out.find("relax_cell").unwrap().boolean().unwrap(), true);
        assert_eq!(out.find("calculator").unwrap().string().unwrap(), "{\"potential\":\"LennardJones\",\"epsilon\":0.20000000000000001,\"sigma\":2.5,\"cutoff\":7}");
        assert!(!w.open && w.step.is_none());
        // standardize-cell keeps its cell choice.
        let mut parameters = Json::object();
        parameters.set("cell", "conventional");
        assert!(w.edit_step(StepEdit { step: "standardize-cell".into(), parameters, ..Default::default() }));
        assert_eq!(w.step_cell, 1);
        w.step_cell = 2;
        w.update_step().unwrap();
        assert_eq!(w.updated.take().unwrap().parameters.find("cell").unwrap().string().unwrap(), "symmetrized");
        assert!(!w.edit_step(StepEdit { step: "replicate".into(), ..Default::default() }));
    }
}
