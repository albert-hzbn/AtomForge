//! The non-destructive modification pipeline (the C++ pipeline/Pipeline): an
//! input structure flows through an ordered list of modifiers, like a pipe
//!     replicate 2 2 1 | select-expression "fz > 0.5" | delete-selected
//! Modifiers act on the whole structure or on the current selection.

use std::sync::{Arc, OnceLock, RwLock};
use std::time::Instant;

use af_core::format::{format_g, stod};
use af_core::json::Json;
use af_core::suggest;
use af_core::task::checkpoint;
use af_core::Structure;

/// What flows between modifiers.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PipelineData {
    pub structure: Structure,
    /// Per atom, 1 = selected (same size as atoms).
    pub selected: Vec<u8>,
    /// Messages from the last modifier.
    pub notes: Vec<String>,
}

impl PipelineData {
    pub fn new(structure: Structure) -> Self {
        let n = structure.atoms.len();
        Self { structure, selected: vec![0; n], notes: Vec::new() }
    }
    pub fn selected_count(&self) -> usize {
        self.selected.iter().filter(|&&s| s == 1).count()
    }
    pub fn clear_selection(&mut self) {
        self.selected = vec![0; self.structure.atoms.len()];
    }
}

/// Parameter kinds: int, float, bool, vector (3 numbers), ivec3 (3 integers),
/// string, element, elements (a list; takes the rest of its step in the text
/// syntax), list, expression, options (likewise), choice (one of `options`,
/// '|'-separated).
#[derive(Clone, Debug, PartialEq)]
pub struct ModifierParameter {
    pub name: String,
    pub label: String,
    pub kind: String,
    /// Default, in the text syntax ("2 2 1", "true", "x > 0").
    pub value: String,
    /// Choice values, else "".
    pub options: String,
}

impl ModifierParameter {
    pub fn new(name: &str, label: &str, kind: &str, value: &str, options: &str) -> Self {
        Self { name: name.into(), label: label.into(), kind: kind.into(), value: value.into(), options: options.into() }
    }
}

/// Shorthand for parameter tables: (name, label, kind, default, options).
pub fn params(list: &[(&str, &str, &str, &str, &str)]) -> Vec<ModifierParameter> {
    list.iter().map(|&(n, l, k, v, o)| ModifierParameter::new(n, l, k, v, o)).collect()
}

pub type ApplyFn = dyn Fn(&Json, &mut PipelineData) -> Result<(), String> + Send + Sync;

#[derive(Clone)]
pub struct ModifierType {
    pub id: String,
    pub title: String,
    /// Selection, Modification, Cell, Analysis or Simulation.
    pub category: String,
    pub help: String,
    pub parameters: Vec<ModifierParameter>,
    pub apply: Arc<ApplyFn>,
}

impl std::fmt::Debug for ModifierType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ModifierType").field("id", &self.id).field("category", &self.category).finish()
    }
}

impl ModifierType {
    pub fn new(
        id: &str,
        title: &str,
        category: &str,
        help: &str,
        parameters: Vec<ModifierParameter>,
        apply: impl Fn(&Json, &mut PipelineData) -> Result<(), String> + Send + Sync + 'static,
    ) -> Self {
        Self { id: id.into(), title: title.into(), category: category.into(), help: help.into(), parameters, apply: Arc::new(apply) }
    }

    pub fn find_parameter(&self, name: &str) -> Option<&ModifierParameter> {
        self.parameters.iter().find(|p| p.name == name)
    }
}

/// The families of built-in modifier types, gathered by the registry. Each
/// family registers here as its crate is ported (analysis, generator, edit,
/// shape, Python and science steps in later phases).
pub type Family = fn() -> Vec<ModifierType>;

static FAMILIES: RwLock<Vec<Family>> = RwLock::new(Vec::new());
static REGISTRY: OnceLock<RwLock<Vec<Arc<ModifierType>>>> = OnceLock::new();

fn category_rank(category: &str) -> usize {
    ["Selection", "Modification", "Cell", "Analysis", "Simulation"].iter().position(|&c| c == category).unwrap_or(5)
}

fn registry() -> &'static RwLock<Vec<Arc<ModifierType>>> {
    REGISTRY.get_or_init(|| {
        let mut types = crate::steps::basic::types();
        types.extend(crate::steps::analysis::types());
        // The C++ order: basic, analysis, generator, edit, shape, Python, science.
        types.extend(crate::steps::generator::types());
        types.extend(crate::steps::edit::types());
        types.extend(crate::steps::shape::types());
        types.extend(crate::steps::python::types());
        types.extend(crate::steps::science::types());
        for family in FAMILIES.read().unwrap().iter() {
            types.extend(family());
        }
        // Menus list steps by category in this order, each category together.
        types.sort_by_key(|t| category_rank(&t.category));
        RwLock::new(types.into_iter().map(Arc::new).collect())
    })
}

/// Adds a family of built-in types; call before the first registry use.
pub fn add_family(family: Family) {
    assert!(REGISTRY.get().is_none(), "modifier families must be added before the registry is used");
    FAMILIES.write().unwrap().push(family);
}

/// All modifier types, in menu order.
pub fn modifier_types() -> Vec<Arc<ModifierType>> {
    registry().read().unwrap().clone()
}

pub fn find_modifier_type(id: &str) -> Option<Arc<ModifierType>> {
    registry().read().unwrap().iter().find(|t| t.id == id).cloned()
}

/// Adds (or replaces, by id) a modifier type, as the desktop and command line
/// register the Build and Edit operations.
pub fn register_modifier_type(t: ModifierType) {
    let mut types = registry().write().unwrap();
    if let Some(existing) = types.iter_mut().find(|e| e.id == t.id) {
        *existing = Arc::new(t);
    } else {
        types.push(Arc::new(t));
    }
}

/// For an error about an unknown parameter: " (did you mean x?); it takes a, b".
pub fn parameter_hint(t: &ModifierType, name: &str) -> String {
    let names: Vec<&str> = t.parameters.iter().map(|p| p.name.as_str()).collect();
    if names.is_empty() {
        return "; it takes none".into();
    }
    format!("{}; it takes {}", suggest::did_you_mean(name, &names), suggest::list(&names))
}

pub fn require_type(id: &str) -> Result<Arc<ModifierType>, String> {
    if let Some(t) = find_modifier_type(id) {
        return Ok(t);
    }
    let ids: Vec<String> = modifier_types().iter().map(|t| t.id.clone()).collect();
    Err(format!("Unknown modifier '{id}'{}", suggest::did_you_mean(id, &ids)))
}

fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn words(text: &str) -> Vec<&str> {
    text.split(|c: char| c.is_ascii() && is_space(c as u8)).filter(|w| !w.is_empty()).collect()
}

/// Parses one value in the text syntax for a parameter kind.
pub fn value_for(p: &ModifierParameter, text: &str) -> Result<Json, String> {
    let kind = p.kind.as_str();
    let number = |token: &str| -> Result<f64, String> {
        match stod(token) {
            Ok((value, used)) if used == token.len() => Ok(value),
            _ => Err(format!("Parameter {}: '{token}' is not a number", p.name)),
        }
    };
    match kind {
        "int" | "float" => {
            let value = number(text)?;
            if kind == "int" && value != value.floor() {
                return Err(format!("Parameter {} must be an integer", p.name));
            }
            Ok(Json::Number(value))
        }
        "bool" => match text {
            "true" | "1" | "yes" | "on" => Ok(Json::Bool(true)),
            "false" | "0" | "no" | "off" => Ok(Json::Bool(false)),
            _ => Err(format!("Parameter {} must be true or false", p.name)),
        },
        "vector" | "ivec3" => {
            let spaced = text.replace(',', " ");
            let parts = words(&spaced);
            if parts.len() != 3 {
                return Err(format!("Parameter {} needs three numbers", p.name));
            }
            let mut list = Vec::new();
            for part in parts {
                let value = number(part)?;
                if kind == "ivec3" && value != value.floor() {
                    return Err(format!("Parameter {} needs integers", p.name));
                }
                list.push(Json::Number(value));
            }
            Ok(Json::Array(list))
        }
        _ => {
            if kind == "choice" && !format!("|{}|", p.options).contains(&format!("|{text}|")) {
                return Err(format!("Parameter {} must be one of {}", p.name, p.options));
            }
            Ok(Json::String(text.to_string()))
        }
    }
}

/// The text form of a value (numbers as a stream with precision 12).
pub fn text_for(p: &ModifierParameter, value: &Json) -> Result<String, String> {
    Ok(match p.kind.as_str() {
        "bool" => (if value.boolean()? { "true" } else { "false" }).to_string(),
        "int" | "float" => format_g(value.number()?, 12),
        "vector" | "ivec3" => {
            let items = value.items()?;
            // A missing component reads as a non-number ("Expected a number").
            let item = |k: usize| -> Result<f64, String> { items.get(k).unwrap_or(&Json::Null).number() };
            format!("{},{},{}", format_g(item(0)?, 12), format_g(item(1)?, 12), format_g(item(2)?, 12))
        }
        _ => value.string()?.to_string(),
    })
}

/// Splits on '|' outside quotes; "||" is the logical operator, not a pipe.
fn steps(text: &str) -> Result<Vec<String>, String> {
    let bytes = text.as_bytes();
    let mut out: Vec<Vec<u8>> = vec![Vec::new()];
    let mut quote = 0u8;
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        let last = out.last_mut().unwrap();
        if quote != 0 {
            if c == quote {
                quote = 0;
            }
            last.push(c);
        } else if c == b'"' || c == b'\'' {
            quote = c;
            last.push(c);
        } else if c == b'|' {
            if bytes.get(i + 1) == Some(&b'|') {
                last.extend_from_slice(b"||");
                i += 1;
            } else {
                out.push(Vec::new());
            }
        } else {
            last.push(c);
        }
        i += 1;
    }
    if quote != 0 {
        return Err("Unterminated quote in the pipeline text".into());
    }
    Ok(out.into_iter().map(|s| String::from_utf8_lossy(&s).into_owned()).collect())
}

/// Tokens of one step: whitespace-separated, quotes group (and are removed).
fn tokens(step: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut current = Vec::new();
    let mut in_token = false;
    let mut quote = 0u8;
    for &c in step.as_bytes() {
        if quote != 0 {
            if c == quote {
                quote = 0;
            } else {
                current.push(c);
            }
            continue;
        }
        if c == b'"' || c == b'\'' {
            quote = c;
            in_token = true;
            continue;
        }
        if is_space(c) {
            if in_token {
                out.push(String::from_utf8_lossy(&current).into_owned());
                current.clear();
                in_token = false;
            }
            continue;
        }
        current.push(c);
        in_token = true;
    }
    if in_token {
        out.push(String::from_utf8_lossy(&current).into_owned());
    }
    out
}

/// One step of a pipeline: its type, parameter values (JSON keyed by name),
/// whether it is applied and an optional display name.
#[derive(Clone, Debug, PartialEq)]
pub struct Modifier {
    pub type_id: String,
    pub parameters: Json,
    pub enabled: bool,
    pub label: String,
}

/// A modifier of `id` with default parameters.
pub fn make_modifier(id: &str) -> Result<Modifier, String> {
    let t = require_type(id)?;
    let mut parameters = Json::object();
    for p in &t.parameters {
        parameters.set(&p.name, value_for(p, &p.value)?);
    }
    Ok(Modifier { type_id: id.into(), parameters, enabled: true, label: String::new() })
}

/// Parameter value or its default; an error for unknown names.
pub fn parameter(parameters: &Json, t: &ModifierType, name: &str) -> Result<Json, String> {
    if let Some(value) = parameters.find(name) {
        return Ok(value.clone());
    }
    let Some(p) = t.find_parameter(name) else {
        return Err(format!("{} has no parameter {name}{}", t.id, parameter_hint(t, name)));
    };
    value_for(p, &p.value)
}

/// Outcome of one modifier during evaluation.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StageResult {
    pub atoms: usize,
    pub selected: usize,
    pub milliseconds: f64,
    pub error: String,
    pub notes: Vec<String>,
    /// Disabled, or after an error.
    pub skipped: bool,
}

/// Outputs of each stage kept between evaluations: editing modifier i only
/// re-runs modifiers i and later.
#[derive(Clone, Debug, Default)]
pub struct PipelineCache {
    pub outputs: Vec<PipelineData>,
    pub stages: Vec<StageResult>,
    /// Leading outputs that are still correct.
    pub valid: usize,
}

impl PipelineCache {
    pub fn invalidate_from(&mut self, index: usize) {
        self.valid = self.valid.min(index);
    }
    pub fn clear(&mut self) {
        *self = Self::default();
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Pipeline {
    pub modifiers: Vec<Modifier>,
}

fn takes_rest(kind: &str) -> bool {
    matches!(kind, "expression" | "elements" | "options" | "list")
}

/// The rest of a step as one value; words that were quoted (they contain
/// spaces) keep their quotes for command-line options.
fn join_rest(kind: &str, mut value: String, words: &[String]) -> String {
    for word in words {
        let spaced = kind == "options" && word.contains([' ', '\t']);
        let mark = if word.contains('"') { '\'' } else { '"' };
        if !value.is_empty() {
            value.push(' ');
        }
        if spaced {
            value.push(mark);
            value.push_str(word);
            value.push(mark);
        } else {
            value.push_str(word);
        }
    }
    value
}

impl Pipeline {
    /// Runs the enabled modifiers on `input`. A failing modifier stops the
    /// pipeline: its stage records the error and the output is the data
    /// before it. Only cancellation fails the whole evaluation.
    pub fn evaluate(&self, input: &Structure) -> Result<(PipelineData, Vec<StageResult>), String> {
        let mut cache = PipelineCache::default();
        let out = self.evaluate_cached(input, &mut cache)?;
        Ok((out, cache.stages))
    }

    pub fn evaluate_cached(&self, input: &Structure, cache: &mut PipelineCache) -> Result<PipelineData, String> {
        let n = self.modifiers.len();
        cache.outputs.resize_with(n, PipelineData::default);
        cache.stages.resize_with(n, StageResult::default);
        cache.valid = cache.valid.min(n);
        let mut data = if cache.valid > 0 { cache.outputs[cache.valid - 1].clone() } else { PipelineData::new(input.clone()) };
        let mut failed = cache.stages[..cache.valid].iter().any(|s| !s.error.is_empty());
        for i in cache.valid..n {
            if let Err(e) = checkpoint() {
                // Stages from here on are not valid.
                cache.valid = i;
                return Err(e);
            }
            let mut stage = StageResult::default();
            let modifier = &self.modifiers[i];
            if failed || !modifier.enabled {
                stage.skipped = true;
            } else {
                let start = Instant::now();
                let result = require_type(&modifier.type_id).and_then(|t| {
                    let mut next = data.clone();
                    next.notes.clear();
                    (t.apply)(&modifier.parameters, &mut next)?;
                    let atoms = next.structure.atoms.len();
                    if next.selected.len() != atoms {
                        next.selected.resize(atoms, 0);
                    }
                    Ok(next)
                });
                match result {
                    Ok(next) => {
                        data = next;
                        stage.notes = data.notes.clone();
                    }
                    Err(error) => {
                        stage.error = error;
                        failed = true;
                    }
                }
                stage.milliseconds = start.elapsed().as_secs_f64() * 1000.0;
            }
            stage.atoms = data.structure.atoms.len();
            stage.selected = data.selected_count();
            cache.outputs[i] = data.clone();
            cache.stages[i] = stage;
        }
        cache.valid = n;
        Ok(if n == 0 { data } else { cache.outputs[n - 1].clone() })
    }

    pub fn to_json(&self) -> Json {
        let mut list = Vec::new();
        for m in &self.modifiers {
            let mut item = Json::object();
            item.set("type", m.type_id.as_str());
            item.set("enabled", m.enabled);
            if !m.label.is_empty() {
                item.set("label", m.label.as_str());
            }
            item.set("parameters", m.parameters.clone());
            list.push(item);
        }
        let mut json = Json::object();
        json.set("format", "atomforge-pipeline");
        json.set("version", 1);
        json.set("modifiers", Json::Array(list));
        json
    }

    pub fn from_json(json: &Json) -> Result<Pipeline, String> {
        if !json.is_object() || !json.contains("modifiers") {
            return Err("A pipeline needs a \"modifiers\" list".into());
        }
        let mut pipeline = Pipeline::default();
        for item in json.at("modifiers")?.items()? {
            let mut modifier = make_modifier(item.at("type")?.string()?)?;
            let t = require_type(&modifier.type_id)?;
            if let Some(Json::Object(members)) = item.find("parameters") {
                for (name, value) in members {
                    let Some(p) = t.find_parameter(name) else {
                        return Err(format!("{} has no parameter {name}{}", modifier.type_id, parameter_hint(&t, name)));
                    };
                    // Values are validated through the text form, which checks kinds and choices.
                    let text = match value {
                        Json::String(s) => s.clone(),
                        _ => text_for(p, value)?,
                    };
                    modifier.parameters.set(name, value_for(p, &text)?);
                }
            }
            if let Some(enabled) = item.find("enabled") {
                modifier.enabled = enabled.boolean()?;
            }
            if let Some(label) = item.find("label") {
                modifier.label = label.string()?.to_string();
            }
            pipeline.modifiers.push(modifier);
        }
        Ok(pipeline)
    }

    /// The shell-pipe syntax: modifiers separated by '|', each
    /// "id positional... key=value..."; vector parameters take three numbers,
    /// an expression takes the rest of its step.
    pub fn parse(text: &str) -> Result<Pipeline, String> {
        let mut pipeline = Pipeline::default();
        let all = steps(text)?;
        for step in &all {
            let parts = tokens(step);
            if parts.is_empty() {
                if step.bytes().all(|c| matches!(c, b' ' | b'\t' | b'\r' | b'\n')) && all.len() == 1 {
                    break;
                }
                return Err("Empty step in the pipeline text".into());
            }
            let mut modifier = make_modifier(&parts[0])?;
            let t = require_type(&modifier.type_id)?;
            let mut positional = 0usize;
            let mut pending: Vec<String> = Vec::new();
            let mut named: Vec<String> = Vec::new();
            // key=value (but not comparisons such as x==1 or x<=1 inside expressions)
            let key_of = |token: &str| -> Option<&ModifierParameter> {
                let b = token.as_bytes();
                let equals = token.find('=')?;
                if equals == 0 || matches!(b[equals - 1], b'<' | b'>' | b'!') || b.get(equals + 1) == Some(&b'=') {
                    return None;
                }
                t.find_parameter(&token[..equals])
            };
            let mut index = 1;
            while index < parts.len() {
                let token = &parts[index];
                if let Some(p) = key_of(token) {
                    named.push(p.name.clone());
                    let mut value = token[token.find('=').unwrap() + 1..].to_string();
                    if takes_rest(&p.kind) {
                        // options=--core-radius 8 --core-element Cu: the following
                        // words too, up to the next name=value of this step.
                        let mut end = index + 1;
                        while end < parts.len() && key_of(&parts[end]).is_none() {
                            end += 1;
                        }
                        value = join_rest(&p.kind, value, &parts[index + 1..end]);
                        index = end - 1;
                    }
                    modifier.parameters.set(&p.name, value_for(p, &value)?);
                } else {
                    pending.push(token.clone());
                }
                index += 1;
            }
            let mut t_index = 0;
            while t_index < pending.len() {
                // Positional values fill the parameters not given by name, in order.
                while positional < t.parameters.len() && named.contains(&t.parameters[positional].name) {
                    positional += 1;
                }
                if positional >= t.parameters.len() {
                    return Err(format!("{}: too many values ({})", t.id, pending[t_index]));
                }
                let p = &t.parameters[positional];
                positional += 1;
                let kind = p.kind.as_str();
                // name=value with a name the step does not have is a misspelt
                // parameter, not a value.
                let word = &pending[t_index];
                if let Some(equals) = word.find('=') {
                    let first_other = word.bytes().position(|c| !(c.is_ascii_lowercase() || c == b'_'));
                    if !takes_rest(kind) && kind != "string" && equals > 0 && first_other == Some(equals) {
                        let name = &word[..equals];
                        return Err(format!("{} has no parameter {name}{}", t.id, parameter_hint(&t, name)));
                    }
                }
                let value;
                if takes_rest(kind) {
                    value = join_rest(kind, String::new(), &pending[t_index..]);
                    t_index = pending.len();
                } else if (kind == "vector" || kind == "ivec3") && !word.contains(',') {
                    if t_index + 3 > pending.len() {
                        return Err(format!("{}: {} needs three numbers", t.id, p.name));
                    }
                    value = format!("{} {} {}", pending[t_index], pending[t_index + 1], pending[t_index + 2]);
                    t_index += 3;
                } else {
                    value = word.clone();
                    t_index += 1;
                }
                modifier.parameters.set(&p.name, value_for(p, &value)?);
            }
            pipeline.modifiers.push(modifier);
        }
        Ok(pipeline)
    }

    pub fn to_text(&self) -> Result<String, String> {
        let mut out = String::new();
        for m in self.modifiers.iter().filter(|m| m.enabled) {
            let t = require_type(&m.type_id)?;
            if !out.is_empty() {
                out.push_str(" | ");
            }
            out.push_str(&m.type_id);
            for p in &t.parameters {
                let value = parameter(&m.parameters, &t, &p.name)?;
                let text = text_for(p, &value)?;
                if text == text_for(p, &value_for(p, &p.value)?)? {
                    continue; // defaults are implied
                }
                let quote = text.is_empty() || text.contains([' ', '\t', '|', '"', '\'']);
                let mark = if text.contains('"') { '\'' } else { '"' };
                if quote {
                    out.push_str(&format!(" {}={mark}{text}{mark}", p.name));
                } else {
                    out.push_str(&format!(" {}={text}", p.name));
                }
            }
        }
        Ok(out)
    }
}
