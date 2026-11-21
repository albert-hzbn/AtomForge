//! AtomForge projects (.afproject, the C++ model/Workspace): every tab's
//! structure, volumetric data, tables, camera, settings, tool and pipeline
//! state and look, in a little-endian binary archive (formats 1-7; numbers
//! as doubles, per-atom values bit-exact). Written to a temporary file and
//! renamed, so a failed save leaves the old project.

use std::collections::BTreeMap;

use af_core::structure::{AtomSite, Structure};
use glam::{DMat3, DVec3};

use crate::volume::{Grid, Mesh, Site, Volume};

/// One tab of a project (portable data only).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Workspace {
    pub structure: Structure,
    pub volume: Volume,
    pub reference: Volume,
    pub surface: Mesh,
    pub table: Vec<f64>,
    pub camera: Vec<f64>,
    pub settings: BTreeMap<String, f64>,
    pub slice_settings: BTreeMap<String, f64>,
    pub columns: i32,
    pub heading: String,
    pub source_path: String,
    pub reference_path: String,
    pub title: String,
    pub history: Vec<String>,
    /// Scientific tools dialog state (JSON), format 2.
    pub science: String,
    /// The structure pipeline (JSON) and, while active, its unmodified input (format 3).
    pub pipeline: String,
    pub has_pipeline_input: bool,
    pub pipeline_input: Structure,
    /// The look (display settings, element colours and radii; format 6).
    pub display: String,
}

pub const FORMAT: i32 = 7;

type R<T> = Result<T, String>;

/// Bytes per element the C++ allows for (sizeof(T)); bounds the array length.
mod size {
    pub const ATOM: u64 = 112;
    pub const SITE: u64 = 32;
    pub const GRID: u64 = 152;
    pub const VEC3D: u64 = 24;
    pub const VEC3F: u64 = 12;
    pub const DOUBLE: u64 = 8;
    pub const INT: u64 = 4;
    pub const STRING: u64 = 32;
    pub const WORKSPACE: u64 = 1024;
}

struct Writer {
    out: Vec<u8>,
    version: i32,
}

impl Writer {
    fn bits(&mut self, v: u64) {
        self.out.extend_from_slice(&v.to_le_bytes());
    }
    fn f64(&mut self, v: f64) -> R<()> {
        if !v.is_finite() {
            return Err("Non-finite project value".into());
        }
        self.bits(v.to_bits());
        Ok(())
    }
    fn f32(&mut self, v: f32) -> R<()> {
        if !v.is_finite() {
            return Err("Float out of range".into());
        }
        self.f64(v as f64)
    }
    fn int(&mut self, v: i32) -> R<()> {
        self.f64(v as f64)
    }
    fn boolean(&mut self, v: bool) -> R<()> {
        self.int(v as i32)
    }
    fn count(&mut self, n: usize, element: u64) -> R<()> {
        if n as u64 > 100_000_000 || n as u64 > 2_147_483_648 / element {
            return Err("Project array too large".into());
        }
        self.bits(n as u64);
        Ok(())
    }
    fn string(&mut self, s: &str, limit: u64) -> R<()> {
        if s.len() as u64 > limit {
            return Err("Project string too large".into());
        }
        self.bits(s.len() as u64);
        self.out.extend_from_slice(s.as_bytes());
        Ok(())
    }
    fn dvec3(&mut self, v: DVec3) -> R<()> {
        self.f64(v.x)?;
        self.f64(v.y)?;
        self.f64(v.z)
    }
    fn map(&mut self, m: &BTreeMap<String, f64>) -> R<()> {
        if m.len() > 10000 {
            return Err("Too many project settings".into());
        }
        self.bits(m.len() as u64);
        for (k, &v) in m {
            self.string(k, 1 << 20)?;
            self.f64(v)?;
        }
        Ok(())
    }
    fn structure(&mut self, s: &Structure) -> R<()> {
        if (!s.grain_colors.is_empty() && s.grain_colors.len() != s.atoms.len()) || (!s.grain_region_ids.is_empty() && s.grain_region_ids.len() != s.atoms.len()) {
            return Err("Misaligned project atom metadata".into());
        }
        self.count(s.atoms.len(), size::ATOM)?;
        for a in &s.atoms {
            self.string(&a.symbol, 1 << 20)?;
            self.int(a.atomic_number)?;
            self.f64(a.x)?;
            self.f64(a.y)?;
            self.f64(a.z)?;
            self.f32(a.r)?;
            self.f32(a.g)?;
            self.f32(a.b)?;
        }
        self.boolean(s.has_unit_cell)?;
        for row in &s.cell_vectors {
            for &v in row {
                self.f64(v)?;
            }
        }
        for &v in &s.cell_offset {
            self.f64(v)?;
        }
        self.f32(s.pbc_boundary_tol)?;
        self.count(s.grain_colors.len(), size::VEC3F)?;
        for c in &s.grain_colors {
            for &v in c {
                self.f32(v)?;
            }
        }
        self.count(s.grain_region_ids.len(), size::INT)?;
        for &v in &s.grain_region_ids {
            self.int(v)?;
        }
        self.string(&s.ipf_load_status, 1 << 20)?;
        self.count(s.dislocation_loop_points.len(), size::VEC3F)?;
        for p in &s.dislocation_loop_points {
            for &v in p {
                self.f32(v)?;
            }
        }
        self.boolean(s.dislocation_detection_done)
    }
    fn fixed(&mut self, s: &Structure) -> R<()> {
        self.bits(s.atoms.len() as u64);
        for a in &s.atoms {
            self.int(a.fixed as i32)?;
        }
        Ok(())
    }
    fn periodicity(&mut self, s: &Structure) -> R<()> {
        for &p in &s.pbc {
            self.boolean(p)?;
        }
        Ok(())
    }
    fn exact(&mut self, values: &[f64]) {
        for &v in values {
            self.bits(v.to_bits());
        }
    }
    fn frozen(&mut self, s: &Structure) -> R<()> {
        if !s.frozen_property.is_empty() && s.frozen_property.len() != s.atoms.len() {
            return Err("Misaligned project frozen property".into());
        }
        self.bits(s.frozen_property.len() as u64);
        self.exact(&s.frozen_property);
        self.string(&s.frozen_property_name, 1 << 20)
    }
    fn atom_data(&mut self, s: &Structure) {
        self.bits(s.atoms.len() as u64);
        for a in &s.atoms {
            self.exact(&[a.magmom, a.charge, a.velocity[0], a.velocity[1], a.velocity[2]]);
        }
    }
    fn grid(&mut self, g: &Grid) -> R<()> {
        for &n in &g.shape {
            self.int(n)?;
        }
        for c in [g.cell.x_axis, g.cell.y_axis, g.cell.z_axis] {
            self.dvec3(c)?;
        }
        self.dvec3(g.origin)?;
        self.boolean(g.periodic)?;
        self.string(&g.name, 1 << 20)?;
        self.string(&g.unit, 1 << 20)?;
        self.count(g.values.len(), size::DOUBLE)?;
        for &v in &g.values {
            self.f64(v)?;
        }
        g.validate()
    }
    fn volume(&mut self, v: &Volume) -> R<()> {
        self.count(v.sites.len(), size::SITE)?;
        for s in &v.sites {
            self.int(s.number)?;
            self.dvec3(s.position)?;
        }
        self.count(v.fields.len(), size::GRID)?;
        for g in &v.fields {
            self.grid(g)?;
        }
        Ok(())
    }
    fn mesh(&mut self, m: &Mesh) -> R<()> {
        self.count(m.vertices.len(), size::VEC3D)?;
        for &v in &m.vertices {
            self.dvec3(v)?;
        }
        self.count(m.colors.len(), size::DOUBLE)?;
        for &v in &m.colors {
            self.f64(v)?;
        }
        if m.vertices.len() % 3 != 0 || (!m.colors.is_empty() && m.colors.len() != m.vertices.len()) {
            return Err("Invalid project mesh".into());
        }
        Ok(())
    }
    fn doubles(&mut self, v: &[f64]) -> R<()> {
        self.count(v.len(), size::DOUBLE)?;
        for &x in v {
            self.f64(x)?;
        }
        Ok(())
    }
    fn workspace(&mut self, w: &Workspace) -> R<()> {
        self.structure(&w.structure)?;
        self.volume(&w.volume)?;
        self.volume(&w.reference)?;
        self.mesh(&w.surface)?;
        self.doubles(&w.table)?;
        self.int(w.columns)?;
        self.string(&w.heading, 1 << 20)?;
        self.map(&w.settings)?;
        self.map(&w.slice_settings)?;
        self.doubles(&w.camera)?;
        self.string(&w.source_path, 1 << 20)?;
        self.string(&w.reference_path, 1 << 20)?;
        self.string(&w.title, 1 << 20)?;
        self.count(w.history.len(), size::STRING)?;
        for h in &w.history {
            self.string(h, 1 << 20)?;
        }
        if self.version >= 2 {
            self.string(&w.science, 1 << 30)?;
            let s = &w.structure;
            if !s.atom_property.is_empty() && s.atom_property.len() != s.atoms.len() {
                return Err("Misaligned project atom property".into());
            }
            self.bits(s.atom_property.len() as u64);
            self.exact(&s.atom_property);
            self.string(&s.atom_property_name, 1 << 20)?;
        }
        if self.version >= 3 {
            self.string(&w.pipeline, 1 << 30)?;
            self.boolean(w.has_pipeline_input)?;
            if w.has_pipeline_input {
                self.structure(&w.pipeline_input)?;
            }
        }
        if self.version >= 4 {
            self.fixed(&w.structure)?;
            if w.has_pipeline_input {
                self.fixed(&w.pipeline_input)?;
            }
        }
        if self.version >= 5 {
            self.periodicity(&w.structure)?;
            if w.has_pipeline_input {
                self.periodicity(&w.pipeline_input)?;
            }
            self.frozen(&w.structure)?;
        }
        if self.version >= 6 {
            self.string(&w.display, 1 << 24)?;
        }
        if self.version >= 7 {
            self.atom_data(&w.structure);
            if w.has_pipeline_input {
                self.atom_data(&w.pipeline_input);
            }
        }
        check_workspace(w)
    }
}

fn check_workspace(w: &Workspace) -> R<()> {
    if !w.camera.is_empty() && w.camera.len() != 7 {
        return Err("Invalid project camera".into());
    }
    if w.camera.iter().any(|c| c.abs() > 1e9) {
        return Err("Project camera outside supported range".into());
    }
    if w.columns < 0 || w.columns > 10000 || (!w.table.is_empty() && (w.columns == 0 || w.table.len() % w.columns as usize != 0)) {
        return Err("Invalid project table".into());
    }
    Ok(())
}

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
    version: i32,
}

impl Reader<'_> {
    fn remaining(&self) -> u64 {
        (self.data.len() - self.pos) as u64
    }
    fn bytes(&mut self, n: usize) -> R<&[u8]> {
        if n as u64 > self.remaining() {
            return Err("Truncated project".into());
        }
        let slice = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(slice)
    }
    fn bits(&mut self) -> R<u64> {
        let b = self.bytes(8)?;
        Ok(u64::from_le_bytes(b.try_into().unwrap()))
    }
    fn f64(&mut self) -> R<f64> {
        let v = f64::from_bits(self.bits()?);
        if !v.is_finite() {
            return Err("Non-finite project value".into());
        }
        Ok(v)
    }
    fn f32(&mut self) -> R<f32> {
        let v = self.f64()? as f32;
        if !v.is_finite() {
            return Err("Float out of range".into());
        }
        Ok(v)
    }
    fn int(&mut self) -> R<i32> {
        let v = self.f64()?;
        if v != v.floor() || v < i32::MIN as f64 || v > i32::MAX as f64 {
            return Err("Invalid project integer".into());
        }
        Ok(v as i32)
    }
    fn boolean(&mut self) -> R<bool> {
        match self.int()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err("Invalid project boolean".into()),
        }
    }
    fn count(&mut self, element: u64) -> R<usize> {
        let n = self.bits()?;
        if n > 100_000_000 || n > 2_147_483_648 / element || n > self.remaining() / 8 {
            return Err("Project array too large".into());
        }
        Ok(n as usize)
    }
    fn string(&mut self, limit: u64) -> R<String> {
        let n = self.bits()?;
        if n > limit || n > self.remaining() {
            return Err("Project string too large".into());
        }
        let b = self.bytes(n as usize)?;
        Ok(String::from_utf8_lossy(b).into_owned())
    }
    fn dvec3(&mut self) -> R<DVec3> {
        Ok(DVec3::new(self.f64()?, self.f64()?, self.f64()?))
    }
    fn map(&mut self) -> R<BTreeMap<String, f64>> {
        let n = self.bits()?;
        if n > 10000 {
            return Err("Too many project settings".into());
        }
        let mut m = BTreeMap::new();
        for _ in 0..n {
            let k = self.string(1 << 20)?;
            let v = self.f64()?;
            if m.insert(k, v).is_some() {
                return Err("Duplicate project setting".into());
            }
        }
        Ok(m)
    }
    fn structure(&mut self) -> R<Structure> {
        let mut s = Structure::default();
        let n = self.count(size::ATOM)?;
        for _ in 0..n {
            let symbol = self.string(1 << 20)?;
            let atomic_number = self.int()?;
            let (x, y, z) = (self.f64()?, self.f64()?, self.f64()?);
            let (r, g, b) = (self.f32()?, self.f32()?, self.f32()?);
            s.atoms.push(AtomSite { symbol, atomic_number, x, y, z, r, g, b, ..Default::default() });
        }
        s.has_unit_cell = self.boolean()?;
        for row in s.cell_vectors.iter_mut() {
            for v in row.iter_mut() {
                *v = self.f64()?;
            }
        }
        for v in s.cell_offset.iter_mut() {
            *v = self.f64()?;
        }
        s.pbc_boundary_tol = self.f32()?;
        let n = self.count(size::VEC3F)?;
        for _ in 0..n {
            s.grain_colors.push([self.f32()?, self.f32()?, self.f32()?]);
        }
        let n = self.count(size::INT)?;
        for _ in 0..n {
            s.grain_region_ids.push(self.int()?);
        }
        s.ipf_load_status = self.string(1 << 20)?;
        let n = self.count(size::VEC3F)?;
        for _ in 0..n {
            s.dislocation_loop_points.push([self.f32()?, self.f32()?, self.f32()?]);
        }
        s.dislocation_detection_done = self.boolean()?;
        if (!s.grain_colors.is_empty() && s.grain_colors.len() != s.atoms.len()) || (!s.grain_region_ids.is_empty() && s.grain_region_ids.len() != s.atoms.len()) {
            return Err("Misaligned project atom metadata".into());
        }
        Ok(s)
    }
    fn fixed(&mut self, s: &mut Structure) -> R<()> {
        if self.bits()? != s.atoms.len() as u64 {
            return Err("Misaligned project fixed atoms".into());
        }
        for a in s.atoms.iter_mut() {
            let f = self.int()?;
            if !(0..=7).contains(&f) {
                return Err("Invalid project fixed atom".into());
            }
            a.fixed = f as u8;
        }
        Ok(())
    }
    fn exact(&mut self) -> R<f64> {
        Ok(f64::from_bits(self.bits()?))
    }
    fn atom_data(&mut self, s: &mut Structure) -> R<()> {
        if self.bits()? != s.atoms.len() as u64 {
            return Err("Misaligned project atom data".into());
        }
        for a in s.atoms.iter_mut() {
            a.magmom = self.exact()?;
            a.charge = self.exact()?;
            for v in a.velocity.iter_mut() {
                *v = self.exact()?;
            }
        }
        Ok(())
    }
    fn grid(&mut self) -> R<Grid> {
        let mut g = Grid { shape: [self.int()?, self.int()?, self.int()?], ..Default::default() };
        g.cell = DMat3::from_cols(self.dvec3()?, self.dvec3()?, self.dvec3()?);
        g.origin = self.dvec3()?;
        g.periodic = self.boolean()?;
        g.name = self.string(1 << 20)?;
        g.unit = self.string(1 << 20)?;
        let n = self.count(size::DOUBLE)?;
        g.values = (0..n).map(|_| self.f64()).collect::<R<_>>()?;
        g.validate()?;
        Ok(g)
    }
    fn volume(&mut self) -> R<Volume> {
        let mut v = Volume::default();
        let n = self.count(size::SITE)?;
        for _ in 0..n {
            v.sites.push(Site { number: self.int()?, position: self.dvec3()? });
        }
        let n = self.count(size::GRID)?;
        for _ in 0..n {
            v.fields.push(self.grid()?);
        }
        Ok(v)
    }
    fn mesh(&mut self) -> R<Mesh> {
        let n = self.count(size::VEC3D)?;
        let vertices = (0..n).map(|_| self.dvec3()).collect::<R<Vec<_>>>()?;
        let n = self.count(size::DOUBLE)?;
        let colors = (0..n).map(|_| self.f64()).collect::<R<Vec<_>>>()?;
        if vertices.len() % 3 != 0 || (!colors.is_empty() && colors.len() != vertices.len()) {
            return Err("Invalid project mesh".into());
        }
        Ok(Mesh { vertices, colors })
    }
    fn doubles(&mut self) -> R<Vec<f64>> {
        let n = self.count(size::DOUBLE)?;
        (0..n).map(|_| self.f64()).collect()
    }
    fn workspace(&mut self) -> R<Workspace> {
        let mut w = Workspace { structure: self.structure()?, volume: self.volume()?, reference: self.volume()?, surface: self.mesh()?, ..Default::default() };
        w.table = self.doubles()?;
        w.columns = self.int()?;
        w.heading = self.string(1 << 20)?;
        w.settings = self.map()?;
        w.slice_settings = self.map()?;
        w.camera = self.doubles()?;
        w.source_path = self.string(1 << 20)?;
        w.reference_path = self.string(1 << 20)?;
        w.title = self.string(1 << 20)?;
        let n = self.count(size::STRING)?;
        for _ in 0..n {
            w.history.push(self.string(1 << 20)?);
        }
        if self.version >= 2 {
            w.science = self.string(1 << 30)?;
            let n = self.bits()?;
            if n != 0 && n != w.structure.atoms.len() as u64 {
                return Err("Misaligned project atom property".into());
            }
            w.structure.atom_property = (0..n).map(|_| self.exact()).collect::<R<_>>()?;
            w.structure.atom_property_name = self.string(1 << 20)?;
        }
        if self.version >= 3 {
            w.pipeline = self.string(1 << 30)?;
            w.has_pipeline_input = self.boolean()?;
            if w.has_pipeline_input {
                w.pipeline_input = self.structure()?;
            }
        }
        if self.version >= 4 {
            self.fixed(&mut w.structure)?;
            if w.has_pipeline_input {
                self.fixed(&mut w.pipeline_input)?;
            }
        }
        if self.version >= 5 {
            for p in w.structure.pbc.iter_mut() {
                *p = self.boolean()?;
            }
            if w.has_pipeline_input {
                for p in w.pipeline_input.pbc.iter_mut() {
                    *p = self.boolean()?;
                }
            }
            let n = self.bits()?;
            if n != 0 && n != w.structure.atoms.len() as u64 {
                return Err("Misaligned project frozen property".into());
            }
            w.structure.frozen_property = (0..n).map(|_| self.exact()).collect::<R<_>>()?;
            w.structure.frozen_property_name = self.string(1 << 20)?;
        }
        if self.version >= 6 {
            w.display = self.string(1 << 24)?;
        }
        if self.version >= 7 {
            self.atom_data(&mut w.structure)?;
            if w.has_pipeline_input {
                self.atom_data(&mut w.pipeline_input)?;
            }
        }
        check_workspace(&w)?;
        Ok(w)
    }
}

/// The bytes of a project in a format (1-7).
pub fn project_bytes(tabs: &[Workspace], format: i32) -> R<Vec<u8>> {
    if tabs.is_empty() || tabs.len() > 1000 {
        return Err("A project requires 1-1000 tabs".into());
    }
    if !(1..=7).contains(&format) {
        return Err("Unsupported project format".into());
    }
    let mut w = Writer { out: Vec::new(), version: format };
    w.string(&format!("ATOMFORGE_PROJECT_{format}"), 1 << 20)?;
    w.count(tabs.len(), size::WORKSPACE)?;
    for t in tabs {
        w.workspace(t)?;
    }
    Ok(w.out)
}

/// Saves a project (saveWorkspace): written beside the file, then renamed over it.
pub fn save_project(tabs: &[Workspace], path: &str, format: i32) -> R<()> {
    let bytes = project_bytes(tabs, format)?;
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0);
    let temporary = format!("{path}.tmp-{stamp}");
    let result = std::fs::write(&temporary, &bytes).map_err(|_| "Cannot write project".to_string()).and_then(|_| std::fs::rename(&temporary, path).map_err(|_| "Cannot replace project file".to_string()));
    if result.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    result
}

/// Reads a project's bytes (loadWorkspace).
pub fn parse_project(data: &[u8]) -> R<Vec<Workspace>> {
    let mut r = Reader { data, pos: 0, version: 3 };
    let magic = r.string(1 << 20)?;
    r.version = match magic.strip_prefix("ATOMFORGE_PROJECT_").and_then(|v| v.parse::<i32>().ok()) {
        Some(v @ 1..=7) if magic == format!("ATOMFORGE_PROJECT_{v}") => v,
        _ => return Err("Unsupported project version (saved by a newer AtomForge?)".into()),
    };
    let n = r.bits()?;
    if n == 0 || n > 1000 {
        return Err("Invalid project tab count".into());
    }
    let mut tabs = Vec::with_capacity(n as usize);
    for _ in 0..n {
        tabs.push(r.workspace()?);
    }
    if r.remaining() != 0 {
        return Err("Unexpected trailing project data".into());
    }
    Ok(tabs)
}

pub fn load_project(path: &str) -> R<Vec<Workspace>> {
    let size = std::fs::metadata(path).map_err(|e| e.to_string())?.len();
    if size > 8 * 1024 * 1024 * 1024 {
        return Err("Project exceeds 8 GiB limit".into());
    }
    let data = std::fs::read(path).map_err(|_| "Cannot read project".to_string())?;
    parse_project(&data)
}
