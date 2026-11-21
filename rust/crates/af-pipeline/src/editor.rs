//! Editing state behind the desktop's Structure Pipeline panel (the C++
//! pipeline/PipelineEditor), without any GUI: the captured input structure,
//! the modifier list and the per-stage cache; each edit invalidates exactly
//! the stages it affects.

use af_core::json::Json;
use af_core::structure::Structure;

use crate::pipeline::{make_modifier, parameter_hint, require_type, Modifier, Pipeline, PipelineCache, PipelineData, StageResult};

#[derive(Clone, Debug, Default)]
pub struct PipelineEditor {
    input: Option<Structure>,
    pipeline: Pipeline,
    cache: PipelineCache,
    output: PipelineData,
    output_print: u64,
    dirty: bool,
}

impl PipelineEditor {
    pub fn new() -> Self {
        Self { dirty: true, ..Default::default() }
    }

    pub fn active(&self) -> bool {
        self.input.is_some()
    }

    /// Starts (or restarts) the pipeline on a copy of `structure`; the modifiers are kept.
    pub fn set_input(&mut self, structure: &Structure) {
        self.input = Some(structure.clone());
        self.changed(0);
    }

    pub fn input(&self) -> Option<&Structure> {
        self.input.as_ref()
    }

    /// Ends the pipeline; returns its last output (the structure to keep).
    pub fn bake(&mut self) -> Structure {
        let out = if self.active() { self.evaluate().map(|d| d.structure.clone()).unwrap_or_default() } else { Structure::default() };
        self.reset();
        out
    }

    pub fn reset(&mut self) {
        *self = Self::new();
    }

    fn changed(&mut self, from: usize) {
        self.cache.invalidate_from(from);
        self.dirty = true;
    }

    pub fn pipeline(&self) -> &Pipeline {
        &self.pipeline
    }

    pub fn size(&self) -> usize {
        self.pipeline.modifiers.len()
    }

    pub fn modifier(&self, index: usize) -> &Modifier {
        &self.pipeline.modifiers[index]
    }

    fn check(&self, index: usize) -> Result<(), String> {
        if index >= self.size() {
            return Err(format!("No modifier {index}"));
        }
        Ok(())
    }

    /// Inserts a modifier of `type_id` with default parameters at `index`
    /// (the end when past it); returns its position.
    pub fn add(&mut self, type_id: &str, index: usize) -> Result<usize, String> {
        let index = index.min(self.size());
        self.pipeline.modifiers.insert(index, make_modifier(type_id)?);
        self.changed(index);
        Ok(index)
    }

    pub fn remove(&mut self, index: usize) -> Result<(), String> {
        self.check(index)?;
        self.pipeline.modifiers.remove(index);
        self.changed(index);
        Ok(())
    }

    /// Moves the modifier at `from` so that it ends up at `to`.
    pub fn move_step(&mut self, from: usize, to: usize) -> Result<(), String> {
        if from >= self.size() || to >= self.size() {
            return Err("Modifier position out of range".into());
        }
        if from == to {
            return Ok(());
        }
        let moving = self.pipeline.modifiers.remove(from);
        self.pipeline.modifiers.insert(to, moving);
        self.changed(from.min(to));
        Ok(())
    }

    pub fn duplicate(&mut self, index: usize) -> Result<usize, String> {
        self.check(index)?;
        let copy = self.pipeline.modifiers[index].clone();
        self.pipeline.modifiers.insert(index + 1, copy);
        self.changed(index + 1);
        Ok(index + 1)
    }

    /// Replaces the modifier at `index` by one of `type_id`, keeping its
    /// enabled state, its name and the parameters both types have (same name and kind).
    pub fn replace(&mut self, index: usize, type_id: &str) -> Result<(), String> {
        self.check(index)?;
        let old = self.pipeline.modifiers[index].clone();
        if old.type_id == type_id {
            return Ok(());
        }
        let mut replacement = make_modifier(type_id)?;
        replacement.enabled = old.enabled;
        replacement.label = old.label.clone();
        let from = require_type(&old.type_id)?;
        let to = require_type(type_id)?;
        for p in &to.parameters {
            for q in &from.parameters {
                if p.name == q.name && p.kind == q.kind {
                    if let Some(v) = old.parameters.find(&p.name) {
                        replacement.parameters.set(&p.name, v.clone());
                    }
                }
            }
        }
        self.pipeline.modifiers[index] = replacement;
        self.changed(index);
        Ok(())
    }

    pub fn set_enabled(&mut self, index: usize, enabled: bool) -> Result<(), String> {
        self.check(index)?;
        if self.pipeline.modifiers[index].enabled != enabled {
            self.pipeline.modifiers[index].enabled = enabled;
            self.changed(index);
        }
        Ok(())
    }

    /// The display name (does not change the output).
    pub fn set_label(&mut self, index: usize, label: &str) -> Result<(), String> {
        self.check(index)?;
        self.pipeline.modifiers[index].label = label.to_string();
        Ok(())
    }

    pub fn set_parameter(&mut self, index: usize, name: &str, value: Json) -> Result<(), String> {
        self.check(index)?;
        let t = require_type(&self.pipeline.modifiers[index].type_id)?;
        if t.find_parameter(name).is_none() {
            return Err(format!("{} has no parameter {name}{}", t.id, parameter_hint(&t, name)));
        }
        self.pipeline.modifiers[index].parameters.set(name, value);
        self.changed(index);
        Ok(())
    }

    /// Replaces the whole modifier list (text syntax or loaded JSON).
    pub fn set_pipeline(&mut self, pipeline: Pipeline) {
        self.pipeline = pipeline;
        self.cache.clear();
        self.changed(0);
    }

    /// Re-runs the stages that changed since the last call; returns the output.
    pub fn evaluate(&mut self) -> Result<&PipelineData, String> {
        let input = self.input.as_ref().ok_or("The pipeline has no input structure")?;
        if self.dirty {
            self.output = self.pipeline.evaluate_cached(input, &mut self.cache)?;
            self.output_print = fingerprint(&self.output.structure);
            self.dirty = false;
        }
        Ok(&self.output)
    }

    pub fn stages(&self) -> &[StageResult] {
        &self.cache.stages[..self.cache.valid.min(self.cache.stages.len())]
    }

    pub fn dirty(&self) -> bool {
        self.dirty
    }

    /// The data after stage `index`, once evaluated.
    pub fn stage_output(&self, index: usize) -> Option<&PipelineData> {
        (index < self.cache.valid).then(|| &self.cache.outputs[index])
    }

    /// Whether `structure` is still the last output (edits made outside the
    /// pipeline, which its next evaluation would overwrite, change it).
    pub fn is_current_output(&self, structure: &Structure) -> bool {
        !self.dirty && fingerprint(structure) == self.output_print
    }

    /// The pipeline and whether it runs, for project files (the input is stored separately).
    pub fn to_json(&self) -> Json {
        let mut json = self.pipeline.to_json();
        json.set("active", Json::from(self.active()));
        json
    }

    pub fn from_json(&mut self, json: &Json) -> Result<(), String> {
        self.set_pipeline(Pipeline::from_json(json)?);
        Ok(())
    }
}

/// FNV-1a over the species, positions and cell (PipelineEditor::fingerprint).
pub fn fingerprint(s: &Structure) -> u64 {
    let mut hash: u64 = 1469598103934665603;
    let mut mix = |bytes: &[u8]| {
        for &b in bytes {
            hash ^= b as u64;
            hash = hash.wrapping_mul(1099511628211);
        }
    };
    mix(&(s.atoms.len() as u64).to_le_bytes());
    for a in &s.atoms {
        mix(a.symbol.as_bytes());
        mix(&a.x.to_le_bytes());
        mix(&a.y.to_le_bytes());
        mix(&a.z.to_le_bytes());
    }
    mix(&[s.has_unit_cell as u8]);
    for row in &s.cell_vectors {
        for v in row {
            mix(&v.to_le_bytes());
        }
    }
    hash
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::structure::AtomSite;

    /// Conventional fcc cells of copper, a = 3.615.
    fn copper(repeat: usize) -> Structure {
        let a = 3.615;
        let r = repeat as f64;
        let mut s = Structure { has_unit_cell: true, cell_vectors: [[a * r, 0.0, 0.0], [0.0, a * r, 0.0], [0.0, 0.0, a * r]], ..Default::default() };
        for i in 0..repeat {
            for j in 0..repeat {
                for k in 0..repeat {
                    for f in [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]] {
                        s.atoms.push(AtomSite {
                            symbol: "Cu".into(),
                            atomic_number: 29,
                            x: (i as f64 + f[0]) * a,
                            y: (j as f64 + f[1]) * a,
                            z: (k as f64 + f[2]) * a,
                            ..Default::default()
                        });
                    }
                }
            }
        }
        s
    }

    #[test]
    fn replace_keeps_state_and_shared_parameters() {
        let mut e = PipelineEditor::new();
        e.set_input(&copper(2));
        e.set_pipeline(Pipeline::parse("select-sphere 0 0 0 2.6 mode=add | delete-selected").unwrap());
        e.set_enabled(0, false).unwrap();
        e.set_label(0, "corner").unwrap();
        e.replace(0, "select-slab").unwrap();
        let m = e.modifier(0);
        assert!(m.type_id == "select-slab" && !m.enabled && m.label == "corner", "type, enabled state and name");
        assert_eq!(m.parameters.at("mode").unwrap().string().unwrap(), "add", "a parameter both types have is kept");
        assert_eq!(m.parameters.at("axis").unwrap().string().unwrap(), "z", "the new type's own parameters take their defaults");
        e.set_enabled(0, true).unwrap();
        e.replace(0, "select-element").unwrap();
        assert!(e.evaluate().unwrap().structure.atoms.is_empty(), "the replaced step runs (all Cu selected and deleted)");
        assert!(e.replace(0, "nope").is_err(), "unknown type");
        assert!(e.replace(5, "wrap").is_err(), "index out of range");
    }

    #[test]
    fn add_move_duplicate_disable_delete_bake() {
        let mut e = PipelineEditor::new();
        assert!(!e.active(), "inactive until it has an input");
        assert!(e.evaluate().is_err(), "evaluation needs an input");
        let input = copper(2);
        e.set_input(&input);
        let same = |e: &mut PipelineEditor, what: &str| {
            let fresh = e.pipeline().evaluate(&input).unwrap().0;
            let cached = e.evaluate().unwrap();
            assert_eq!(cached.structure.atoms.len(), fresh.structure.atoms.len(), "{what}: atom count");
            for (a, b) in cached.structure.atoms.iter().zip(&fresh.structure.atoms) {
                assert!(a.symbol == b.symbol && a.x == b.x, "{what}: species and positions");
            }
            assert_eq!(cached.selected, fresh.selected, "{what}: selection");
        };
        assert_eq!(e.evaluate().unwrap().structure.atoms.len(), 32, "empty pipeline passes the input through");
        let replicate = e.add("replicate", usize::MAX).unwrap();
        e.set_parameter(replicate, "counts", Json::Array(vec![2.0.into(), 1.0.into(), 1.0.into()])).unwrap();
        same(&mut e, "add replicate");
        assert_eq!(e.evaluate().unwrap().structure.atoms.len(), 64, "replicated");
        let sel = e.add("select-expression", usize::MAX).unwrap();
        e.set_parameter(sel, "expression", Json::from("fx < 0.25")).unwrap();
        e.add("delete-selected", usize::MAX).unwrap();
        same(&mut e, "select and delete");
        let after_delete = e.evaluate().unwrap().structure.atoms.len();
        e.add("assign-element", 0).unwrap();
        e.set_parameter(0, "element", Json::from("Ni")).unwrap();
        e.set_parameter(0, "target", Json::from("all")).unwrap();
        same(&mut e, "insert at the front");
        assert_eq!(e.evaluate().unwrap().structure.atoms[0].symbol, "Ni", "assigned before replication");
        e.move_step(0, 3).unwrap();
        same(&mut e, "move to the end");
        assert!(e.modifier(3).type_id == "assign-element" && e.modifier(0).type_id == "replicate", "order after the move");
        e.duplicate(0).unwrap();
        same(&mut e, "duplicate replicate");
        assert_eq!(e.evaluate().unwrap().structure.atoms.len(), 2 * after_delete, "duplicated replicate doubles the cell again");
        e.set_enabled(1, false).unwrap();
        same(&mut e, "disable the copy");
        assert_eq!(e.evaluate().unwrap().structure.atoms.len(), after_delete, "bypassed step has no effect");
        e.remove(1).unwrap();
        same(&mut e, "remove the copy");
        assert_eq!(e.size(), 4);
        assert!(e.set_parameter(0, "bogus", Json::from(1.0)).is_err(), "unknown parameter");
        assert!(e.move_step(0, 9).is_err() && e.remove(9).is_err());
        let output = e.evaluate().unwrap().structure.clone();
        assert!(e.is_current_output(&output), "output recognised");
        let mut edited = output.clone();
        edited.atoms[0].x += 0.5;
        assert!(!e.is_current_output(&edited), "outside edits detected");
        assert_eq!(e.input().unwrap().atoms.len(), 32, "input kept");
        let mut copy = PipelineEditor::new();
        copy.from_json(&e.to_json()).unwrap();
        copy.set_input(&input);
        assert_eq!(copy.evaluate().unwrap().structure.atoms.len(), output.atoms.len(), "restored pipeline gives the same output");
        let baked = e.bake();
        assert!(baked.atoms.len() == output.atoms.len() && !e.active() && e.size() == 0, "bake ends the pipeline with its output");
    }
}
