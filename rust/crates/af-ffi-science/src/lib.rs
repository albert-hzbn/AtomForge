//! C interface of the native scientific tools for the Python package
//! (atomforge.native, via ctypes): the C++ science/NativeAPI under the same
//! names. Requests and results are JSON text, the same documents used by
//! `AtomForge --science`. Returned strings are owned by the caller and
//! released with afs_free; panics never cross the boundary.

use std::ffi::{c_char, c_int, CStr, CString};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;

use af_core::json::Json;
use af_core::structure::Structure;
use af_pipeline::pipeline::{make_modifier, modifier_types, Pipeline};
use af_science::catalog::catalog;
use af_science::data::{structure_from_json, structure_json};
use af_science::tools::run_tool;

fn copy_out(text: &str) -> *mut c_char {
    let bytes: Vec<u8> = text.bytes().map(|b| if b == 0 { b'?' } else { b }).collect();
    CString::new(bytes).map(CString::into_raw).unwrap_or(std::ptr::null_mut())
}

fn input(text: *const c_char, default: &str) -> String {
    if text.is_null() {
        default.to_string()
    } else {
        unsafe { CStr::from_ptr(text) }.to_string_lossy().into_owned()
    }
}

/// Runs `work`, storing *status = 0 and its text, or *status = 1 and the error.
fn with_status(status: *mut c_int, work: impl FnOnce() -> Result<String, String>) -> *mut c_char {
    let set = |v: c_int| {
        if !status.is_null() {
            unsafe { *status = v };
        }
    };
    set(1);
    match catch_unwind(AssertUnwindSafe(work)) {
        Ok(Ok(text)) => {
            set(0);
            copy_out(&text)
        }
        Ok(Err(error)) => copy_out(&error),
        Err(_) => copy_out("Unknown native error"),
    }
}

fn catalog_json() -> Result<Json, String> {
    let mut tools = Vec::new();
    for tool in catalog() {
        let mut entry = Json::object();
        entry.set("id", tool.id);
        entry.set("title", tool.title);
        entry.set("category", tool.category);
        entry.set("help", tool.help);
        let mut parameters = Vec::new();
        for p in tool.parameters {
            let mut item = Json::object();
            item.set("name", p.name);
            item.set("label", p.label);
            item.set("kind", p.kind);
            item.set("required", p.required);
            item.set("default", if p.value.is_empty() { Json::Null } else { Json::parse(p.value)? });
            parameters.push(item);
        }
        entry.set("parameters", Json::Array(parameters));
        tools.push(entry);
    }
    Ok(Json::Array(tools))
}

#[no_mangle]
pub extern "C" fn afs_version() -> c_int {
    1
}

/// # Safety
/// `text` must come from this library (or be null).
#[no_mangle]
pub unsafe extern "C" fn afs_free(text: *mut c_char) {
    if !text.is_null() {
        drop(CString::from_raw(text));
    }
}

/// JSON array of {id, title, category, help, parameters[]}.
#[no_mangle]
pub extern "C" fn afs_catalog() -> *mut c_char {
    match catch_unwind(catalog_json) {
        Ok(Ok(json)) => copy_out(&json.dump(-1)),
        _ => std::ptr::null_mut(),
    }
}

fn modifiers_json() -> Result<Json, String> {
    let mut types = Vec::new();
    for t in modifier_types() {
        let mut entry = Json::object();
        entry.set("id", t.id.as_str());
        entry.set("title", t.title.as_str());
        entry.set("category", t.category.as_str());
        entry.set("help", t.help.as_str());
        let defaults = make_modifier(&t.id)?.parameters;
        let mut parameters = Vec::new();
        for p in &t.parameters {
            let mut item = Json::object();
            item.set("name", p.name.as_str());
            item.set("label", p.label.as_str());
            item.set("kind", p.kind.as_str());
            item.set("default", defaults.at(&p.name)?.clone());
            item.set("options", p.options.as_str());
            parameters.push(item);
        }
        entry.set("parameters", Json::Array(parameters));
        types.push(entry);
    }
    Ok(Json::Array(types))
}

/// JSON array of the structure-pipeline modifiers:
/// {id, title, category, help, parameters: [{name, label, kind, default, options}]}.
#[no_mangle]
pub extern "C" fn afs_pipeline_modifiers() -> *mut c_char {
    match catch_unwind(modifiers_json) {
        Ok(Ok(json)) => copy_out(&json.dump(-1)),
        _ => std::ptr::null_mut(),
    }
}

fn pipeline_of(definition: &Json) -> Result<Pipeline, String> {
    if definition.is_string() {
        Pipeline::parse(definition.string()?)
    } else {
        Pipeline::from_json(definition)
    }
}

fn run_pipeline(request: &str) -> Result<String, String> {
    let input = Json::parse(request)?;
    let pipeline = pipeline_of(input.at("pipeline")?)?;
    // An empty input (no atoms) is valid: builders start from nothing.
    let given = input.at("structure")?;
    let empty = given.is_object() && given.find("symbols").is_some_and(|s| s.is_array() && s.size() == 0);
    let mut structure = if empty { Structure::default() } else { structure_from_json(given)? };
    if let Some(property) = input.find("property") {
        if property.is_array() && property.size() > 0 {
            if property.size() != structure.atoms.len() {
                return Err("property needs one value per atom".into());
            }
            for value in property.items()? {
                structure.atom_property.push(if value.is_number() { value.number()? } else { f64::NAN });
            }
            structure.atom_property_name = match input.find("property_name") {
                Some(name) => name.string()?.to_string(),
                None => "property".into(),
            };
        }
    }
    let (data, stages) = pipeline.evaluate(&structure)?;
    let mut document = Json::object();
    document.set("structure", structure_json(&data.structure));
    document.set("selected", Json::Array(data.selected.iter().map(|&s| Json::from(s != 0)).collect()));
    let s = &data.structure;
    let property: Vec<Json> = if s.atom_property.len() == s.atoms.len() {
        s.atom_property.iter().map(|&v| if v.is_finite() { Json::from(v) } else { Json::Null }).collect()
    } else {
        Vec::new()
    };
    document.set("property", Json::Array(property));
    document.set("property_name", s.atom_property_name.as_str());
    let mut list = Vec::new();
    for (i, stage) in stages.iter().enumerate() {
        let mut item = Json::object();
        item.set("type", pipeline.modifiers[i].type_id.as_str());
        item.set("atoms", stage.atoms);
        item.set("selected", stage.selected);
        item.set("milliseconds", stage.milliseconds);
        item.set("skipped", stage.skipped);
        item.set("error", stage.error.as_str());
        item.set("notes", Json::Array(stage.notes.iter().map(|n| Json::from(n.as_str())).collect()));
        list.push(item);
    }
    document.set("stages", Json::Array(list));
    document.set("pipeline", pipeline.to_json());
    document.set("text", pipeline.to_text()?);
    Ok(document.dump(-1))
}

/// Runs a structure pipeline: request {"pipeline": JSON object or pipe text,
/// "structure": structure JSON, "property": [per-atom values] (optional)}.
/// Returns {"structure", "selected", "property", "property_name", "stages",
/// "pipeline", "text"} with *status = 0, or an error message with *status = 1.
///
/// # Safety
/// `request` is null or a NUL-terminated string; `status` is null or writable.
#[no_mangle]
pub unsafe extern "C" fn afs_pipeline_run(request: *const c_char, status: *mut c_int) -> *mut c_char {
    let request = input(request, "{}");
    with_status(status, || run_pipeline(&request))
}

/// Parses pipe text into pipeline JSON (validating every step), or *status = 1 with the message.
///
/// # Safety
/// As afs_pipeline_run.
#[no_mangle]
pub unsafe extern "C" fn afs_pipeline_parse(text: *const c_char, status: *mut c_int) -> *mut c_char {
    let text = input(text, "");
    with_status(status, || Ok(Pipeline::parse(&text)?.to_json().dump(-1)))
}

/// Pipe text for pipeline JSON (validating it), or *status = 1 with the message.
///
/// # Safety
/// As afs_pipeline_run.
#[no_mangle]
pub unsafe extern "C" fn afs_pipeline_text(json: *const c_char, status: *mut c_int) -> *mut c_char {
    let json = input(json, "{}");
    with_status(status, || Pipeline::from_json(&Json::parse(&json)?)?.to_text())
}

fn run(tool: &str, request: &str, base: &str) -> Result<String, String> {
    let output = run_tool(tool, &Json::parse(request)?, Path::new(if base.is_empty() { "." } else { base }), None)?;
    let mut document = Json::object();
    document.set("result", output.result);
    document.set("frames", Json::Array(output.frames.iter().map(structure_json).collect()));
    document.set(
        "velocities",
        Json::Array(output.velocities.iter().map(|frame| Json::Array(frame.iter().map(|v| Json::Array(vec![v[0].into(), v[1].into(), v[2].into()])).collect())).collect()),
    );
    document.set("times_fs", Json::Array(output.times.iter().map(|&t| Json::from(t)).collect()));
    Ok(document.dump(-1))
}

/// Runs a tool. Returns {"result": ..., "frames": [structures], "velocities": ...,
/// "times_fs": ...} with *status = 0, or an error message with *status = 1.
/// Relative file references resolve against base_dir.
///
/// # Safety
/// The strings are null or NUL-terminated; `status` is null or writable.
#[no_mangle]
pub unsafe extern "C" fn afs_run(tool: *const c_char, request: *const c_char, base_dir: *const c_char, status: *mut c_int) -> *mut c_char {
    let (tool, request, base) = (input(tool, ""), input(request, "{}"), input(base_dir, ""));
    with_status(status, || run(&tool, &request, &base))
}
