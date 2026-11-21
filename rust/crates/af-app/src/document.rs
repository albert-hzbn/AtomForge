//! One open structure (a tab of the C++ EditorApplication): the structure,
//! its file, history, selection, camera and the scene drawn from it.

use af_core::display::{DisplaySettings, Look};
use af_core::occlusion::darken_occluded;
use af_core::structure::Structure;
use af_render::colouring::{apply_colour_mode, ColourMode};
use af_render::scene::{build_instance_data, BondStyle, OrbitCamera, Scene};
use af_science::atom_properties::colour_by_property;

use crate::undo::{Snapshot, UndoRedo};
use crate::view::{has_property, ViewOptions};

pub struct Document {
    /// Unique for the run (pipeline step edits find their tab by it).
    pub id: u64,
    pub structure: Structure,
    pub path: Option<String>,
    pub title: String,
    pub history: UndoRedo,
    saved_revision: u64,
    /// Selected instances (periodic images included), as the C++ selectedInstanceIndices.
    pub selection: Vec<usize>,
    pub camera: OrbitCamera,
    /// Pending scroll zoom (Camera::scrollVelocity).
    pub scroll_velocity: f32,
    pub pending_default_view: bool,
    pub element_colours: Vec<[f32; 3]>,
    pub element_radii: Vec<f32>,
    pub element_shininess: Vec<f32>,
    pub scene: Scene,
    /// The scene must be rebuilt from the structure.
    pub scene_dirty: bool,
    /// The renderer holds another document's (or an older) scene.
    pub upload_needed: bool,
    /// The colour mode and history step the atom colours were last set for.
    colours_for: Option<(ColourMode, u64)>,
    /// The project tab this document was restored from: the parts this app
    /// does not edit yet are saved again with the project.
    pub project: Option<af_electronic::project::Workspace>,
    /// Edit > Structure Pipeline: its steps and input (while it runs).
    pub pipeline: af_pipeline::editor::PipelineEditor,
    /// A large dataset (view only): drawn by the cloud renderer instead of
    /// the structure, which stays empty.
    pub cloud: Option<Box<af_render::cloud::CloudRenderer>>,
    /// Frame the cloud on the next frame.
    pub pending_cloud_fit: bool,
}

pub const IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

impl Document {
    /// A tab for a structure, its element tables `look` (the saved default look's).
    pub fn new(structure: Structure, path: Option<String>, look: &Look) -> Self {
        let title = match &path {
            Some(p) => std::path::Path::new(p).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| p.clone()),
            None => "Untitled".to_string(),
        };
        static NEXT_ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
        let mut doc = Self {
            id: NEXT_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
            structure: Structure::default(),
            path,
            title,
            history: UndoRedo::default(),
            saved_revision: 0,
            selection: Vec::new(),
            camera: OrbitCamera::default(),
            scroll_velocity: 0.0,
            pending_default_view: true,
            element_colours: look.colours.clone(),
            element_radii: look.radii.clone(),
            element_shininess: vec![32.0; 119],
            scene: Scene::default(),
            scene_dirty: true,
            upload_needed: true,
            colours_for: None,
            project: None,
            pipeline: af_pipeline::editor::PipelineEditor::new(),
            cloud: None,
            pending_cloud_fit: false,
        };
        // The history starts empty, so the first step is "load the structure".
        doc.history.reset(doc.snapshot());
        doc.structure = structure;
        doc.apply_element_colours();
        doc.history.commit(doc.snapshot());
        doc.saved_revision = doc.history.current_id();
        doc
    }

    /// Atoms take their element's colour (as loaded files do in the C++).
    pub fn apply_element_colours(&mut self) {
        for atom in self.structure.atoms.iter_mut() {
            let z = atom.atomic_number;
            if z >= 0 && (z as usize) < self.element_colours.len() {
                [atom.r, atom.g, atom.b] = self.element_colours[z as usize];
            }
        }
    }

    /// Element colours or the property display changed: colour the atoms again.
    pub fn invalidate_colours(&mut self) {
        self.colours_for = None;
        self.scene_dirty = true;
    }

    /// The colour mode in effect: Atom Property falls back to Element
    /// without a value for every atom.
    pub fn effective_colour_mode(&self, mode: ColourMode) -> ColourMode {
        if mode == ColourMode::AtomProperty && !has_property(&self.structure) {
            ColourMode::Element
        } else {
            mode
        }
    }

    pub fn snapshot(&self) -> Snapshot {
        Snapshot {
            structure: self.structure.clone(),
            element_radii: self.element_radii.clone(),
            element_colours: self.element_colours.clone(),
            element_shininess: self.element_shininess.clone(),
        }
    }

    fn restore(&mut self, s: Snapshot) {
        self.structure = s.structure;
        self.element_radii = s.element_radii;
        self.element_colours = s.element_colours;
        self.element_shininess = s.element_shininess;
        self.selection.clear();
        self.scene_dirty = true;
    }

    /// The history starts again at the current state, saved (a restored project).
    pub fn reset_history(&mut self) {
        self.history.reset(self.snapshot());
        self.saved_revision = self.history.current_id();
        self.selection.clear();
        self.colours_for = None;
        self.scene_dirty = true;
    }

    /// Records the current state as an undoable step.
    pub fn commit(&mut self) {
        let snapshot = self.snapshot();
        self.history.commit(snapshot);
        self.scene_dirty = true;
    }

    /// Replaces the structure (a builder's or an edit's result) as one undoable step.
    pub fn replace_structure(&mut self, structure: Structure) {
        self.structure = structure;
        self.selection.clear();
        self.commit();
    }

    pub fn undo(&mut self) {
        if let Some(s) = self.history.undo().cloned() {
            self.restore(s);
        }
    }

    pub fn redo(&mut self) {
        if let Some(s) = self.history.redo().cloned() {
            self.restore(s);
        }
    }

    pub fn jump_to(&mut self, step: usize) {
        if let Some(s) = self.history.jump_to(step).cloned() {
            self.restore(s);
        }
    }

    /// Unsaved edits: the history moved away from the state last loaded or saved.
    pub fn is_modified(&self) -> bool {
        !self.structure.atoms.is_empty() && self.history.current_id() != self.saved_revision
    }

    pub fn mark_saved(&mut self, path: &str) {
        self.saved_revision = self.history.current_id();
        self.path = Some(path.to_string());
        self.title = std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string());
    }

    /// Saved as part of a project: no unsaved changes, the file name kept.
    pub fn mark_saved_in_place(&mut self) {
        self.saved_revision = self.history.current_id();
    }

    /// Atom indices (in Structure::atoms) of the selected instances, sorted and unique.
    pub fn selected_atoms(&self) -> Vec<usize> {
        let mut atoms: Vec<usize> =
            self.selection.iter().filter_map(|&i| self.scene.instances.atom_indices.get(i).copied()).filter(|&a| a >= 0).map(|a| a as usize).collect();
        atoms.sort_unstable();
        atoms.dedup();
        atoms
    }

    /// Deletes the selected atoms (deleteSelectedAtoms).
    pub fn delete_selected(&mut self) -> usize {
        let atoms = self.selected_atoms();
        if atoms.is_empty() {
            return 0;
        }
        for &a in atoms.iter().rev() {
            if a < self.structure.atoms.len() {
                self.structure.atoms.remove(a);
                if self.structure.grain_colors.len() > a {
                    self.structure.grain_colors.remove(a);
                }
                if self.structure.grain_region_ids.len() > a {
                    self.structure.grain_region_ids.remove(a);
                }
                if self.structure.atom_property.len() > a {
                    self.structure.atom_property.remove(a);
                }
            }
        }
        self.selection.clear();
        self.commit();
        atoms.len()
    }

    /// The scene from the structure (buildRenderData): atoms coloured by the
    /// view's mode, hidden outside a property range, sized by the atom style
    /// and the display's atom size, ambient occlusion, selection tinted.
    pub fn rebuild_scene(&mut self, display: &DisplaySettings, view: &ViewOptions) {
        if let Some(cloud) = &self.cloud {
            // A cloud tab: the scene only carries the cloud's centre and outline.
            self.scene = crate::large_data::cloud_scene(cloud);
            self.scene_dirty = false;
            self.upload_needed = true;
            return;
        }
        let mode = self.effective_colour_mode(view.colour_mode);
        let property = if mode == ColourMode::AtomProperty { Some(colour_by_property(&self.structure.atom_property, &view.property)) } else { None };
        let key = (mode, self.history.current_id());
        if self.colours_for != Some(key) {
            apply_colour_mode(&mut self.structure, &self.element_colours, mode, property.as_ref().map(|(c, _)| c.as_slice()));
            self.colours_for = Some(key);
        }
        let mut data = build_instance_data(&self.structure, false, &IDENTITY, &self.element_radii, &self.element_shininess);
        if let Some((_, visible)) = &property {
            filter_instances(&mut data, visible);
        }
        let scale = view.style.radius_scale() * display.atom_scale;
        if (scale - 1.0).abs() > 1e-4 {
            for s in data.scales.iter_mut() {
                *s *= scale;
            }
        }
        darken_occluded(&data.positions, &data.scales, &mut data.colors, display.ambient_occlusion);
        let bonds = BondStyle::from_display(display);
        let mut scene = Scene::new(data, &bonds);
        // Selected atoms keep their colours and get a halo (Renderer::drawSelectionHalo).
        let selected: std::collections::HashSet<usize> = self.selected_atoms().into_iter().collect();
        scene.selected = scene.instances.atom_indices.iter().enumerate().filter(|(_, &a)| a >= 0 && selected.contains(&(a as usize))).map(|(i, _)| i).collect();
        self.scene = scene;
        self.scene_dirty = false;
        self.upload_needed = true;
    }
}

/// Drops the instances of hidden atoms (filterInstances); atom indices keep picking aligned.
fn filter_instances(data: &mut af_render::scene::InstanceData, visible: &[bool]) {
    let mut kept = 0;
    for i in 0..data.positions.len() {
        let atom = data.atom_indices[i];
        if atom >= 0 && (atom as usize) < visible.len() && !visible[atom as usize] {
            continue;
        }
        data.positions[kept] = data.positions[i];
        data.colors[kept] = data.colors[i];
        data.scales[kept] = data.scales[i];
        data.bond_radii[kept] = data.bond_radii[i];
        data.shininess[kept] = data.shininess[i];
        data.atomic_numbers[kept] = data.atomic_numbers[i];
        data.atom_indices[kept] = data.atom_indices[i];
        kept += 1;
    }
    data.positions.truncate(kept);
    data.colors.truncate(kept);
    data.scales.truncate(kept);
    data.bond_radii.truncate(kept);
    data.shininess.truncate(kept);
    data.atomic_numbers.truncate(kept);
    data.atom_indices.truncate(kept);
}
