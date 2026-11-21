//! Batches of tool runs (the C++ science/Batch): the Cartesian product of
//! parameter sweeps and file matches, with a table of collected result fields.

use std::path::{Component, Path, PathBuf};

use af_core::json::Json;

use crate::data::StructureReader;
use crate::tools::run_tool;

fn matches(pattern: &[u8], text: &[u8]) -> bool {
    let (mut p, mut t) = (0, 0);
    while p < pattern.len() {
        if pattern[p] == b'*' {
            return (t..=text.len()).any(|skip| matches(&pattern[p + 1..], &text[skip..]));
        }
        if t >= text.len() || (pattern[p] != b'?' && pattern[p] != text[t]) {
            return false;
        }
        p += 1;
        t += 1;
    }
    t == text.len()
}

/// std::filesystem::path::lexically_normal: "." removed, "name/.." collapsed,
/// preferred separators.
pub fn lexically_normal(path: &Path) -> PathBuf {
    let mut parts: Vec<Component> = Vec::new();
    for c in path.components() {
        match c {
            Component::CurDir => {}
            Component::ParentDir => match parts.last() {
                Some(Component::Normal(_)) => {
                    parts.pop();
                }
                Some(Component::RootDir) | Some(Component::Prefix(_)) => {}
                _ => parts.push(c),
            },
            _ => parts.push(c),
        }
    }
    parts.iter().collect()
}

fn expand(pattern: &str, base: &Path) -> Result<Vec<PathBuf>, String> {
    let full = lexically_normal(&base.join(pattern));
    let name = full.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    if !name.contains(['*', '?']) {
        return Ok(vec![full]);
    }
    let mut result = Vec::new();
    if let Some(folder) = full.parent() {
        if let Ok(entries) = std::fs::read_dir(folder) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() && matches(name.as_bytes(), entry.file_name().to_string_lossy().as_bytes()) {
                    result.push(path);
                }
            }
        }
    }
    result.sort();
    if result.is_empty() {
        return Err(format!("No files match {}", full.display()));
    }
    Ok(result)
}

fn cell(value: &Json) -> String {
    match value {
        Json::Null => String::new(),
        Json::String(text) => {
            if text.contains([',', '"', '\n']) {
                format!("\"{}\"", text.replace('"', "\"\""))
            } else {
                text.clone()
            }
        }
        Json::Number(v) => af_core::format::format_g(*v, 12),
        Json::Bool(b) => if *b { "true" } else { "false" }.into(),
        other => cell(&Json::String(other.dump(-1))),
    }
}

/// std::stoul on an array index: leading spaces, digits; None without a number.
fn index_of(key: &str) -> Option<usize> {
    let t = key.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let t = t.strip_prefix('+').unwrap_or(t);
    let digits: String = t.chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse().ok()
}

/// Value at a dotted path ("a.b.0.c") in a result, or null.
pub fn value_at(value: &Json, path: &str) -> Json {
    let mut current = value;
    for key in path.split('.') {
        current = match current {
            Json::Array(items) => match index_of(key).and_then(|i| items.get(i)) {
                Some(v) => v,
                None => return Json::Null,
            },
            Json::Object(_) => match current.find(key) {
                Some(v) => v,
                None => return Json::Null,
            },
            _ => return Json::Null,
        };
    }
    current.clone()
}

/// Runs one tool many times. The batch object holds "tool", optional "base"
/// parameters, "sweep": {"param": [values]} (Cartesian product), "files":
/// {"param": "pattern"} (one run per sorted match of * and ? in the file
/// name), "collect": ["key.sub", ...] and "continue_on_error" (default true).
/// Relative paths resolve against `base`. Returns {"tool", "run_count",
/// "failures", "table": {"columns", "rows"}, "runs"}.
pub fn run_batch(batch: &Json, base: &Path, reader: StructureReader) -> Result<Json, String> {
    if !batch.is_object() || !batch.contains("tool") {
        return Err("A batch needs a \"tool\"".into());
    }
    let tool = batch.at("tool")?.string()?.to_string();
    let base_parameters = batch.find("base").cloned().unwrap_or_else(Json::object);
    if !base_parameters.is_object() {
        return Err("\"base\" must be an object".into());
    }
    let continue_on_error = match batch.find("continue_on_error") {
        Some(v) => v.boolean()?,
        None => true,
    };
    let mut axes: Vec<(String, Vec<Json>)> = Vec::new();
    if let Some(sweep) = batch.find("sweep") {
        for (name, values) in sweep.members()? {
            if !values.is_array() || values.size() == 0 {
                return Err(format!("sweep.{name} must be a nonempty list"));
            }
            axes.push((name.clone(), values.items()?.clone()));
        }
    }
    if let Some(files) = batch.find("files") {
        for (name, pattern) in files.members()? {
            let mut references = Vec::new();
            for path in expand(pattern.string()?, base)? {
                let text = path.to_string_lossy().into_owned();
                // *_file parameters take a path; data and structure inputs take a reference.
                if name.len() > 5 && name.ends_with("_file") {
                    references.push(Json::from(text));
                } else {
                    let mut reference = Json::object();
                    reference.set("file", text);
                    references.push(reference);
                }
            }
            axes.push((name.clone(), references));
        }
    }
    let mut collect = Vec::new();
    if let Some(fields) = batch.find("collect") {
        for field in fields.items()? {
            collect.push(field.string()?.to_string());
        }
    }
    let total: usize = axes.iter().map(|a| a.1.len()).product();
    if total > 100000 {
        return Err("The batch has more than 100000 runs".into());
    }
    let (mut runs, mut rows) = (Vec::new(), Vec::new());
    let mut columns: Vec<Json> = axes.iter().map(|a| Json::from(a.0.as_str())).collect();
    columns.extend(collect.iter().map(|c| Json::from(c.as_str())));
    columns.push("error".into());
    let mut failures = 0usize;
    let mut index = vec![0usize; axes.len()];
    for run in 0..total {
        let mut parameters = base_parameters.clone();
        let mut row = Vec::new();
        for (a, axis) in axes.iter().enumerate() {
            let value = &axis.1[index[a]];
            parameters.set(&axis.0, value.clone());
            let shown = match value.find("file") {
                Some(file) if value.is_object() => {
                    Json::from(Path::new(file.string()?).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default())
                }
                _ => value.clone(),
            };
            row.push(shown);
        }
        let mut record = Json::object();
        record.set("parameters", parameters.clone());
        match run_tool(&tool, &parameters, base, reader) {
            Ok(output) => {
                for field in &collect {
                    row.push(value_at(&output.result, field));
                }
                record.set("result", output.result);
                row.push(Json::Null);
            }
            Err(error) => {
                if error == "Calculation cancelled" {
                    return Err(error);
                }
                if !continue_on_error {
                    return Err(format!("Run {}: {error}", run + 1));
                }
                failures += 1;
                record.set("error", error.as_str());
                row.extend(collect.iter().map(|_| Json::Null));
                row.push(error.into());
            }
        }
        runs.push(record);
        rows.push(Json::Array(row));
        for a in (0..axes.len()).rev() {
            index[a] += 1;
            if index[a] < axes[a].1.len() {
                break;
            }
            index[a] = 0;
        }
    }
    let mut table = Json::object();
    table.set("columns", Json::Array(columns));
    table.set("rows", Json::Array(rows));
    let mut result = Json::object();
    result.set("tool", tool);
    result.set("run_count", total);
    result.set("failures", failures);
    result.set("table", table);
    result.set("runs", Json::Array(runs));
    Ok(result)
}

/// CSV text of a batch table (parameter columns then collected values).
pub fn batch_csv(batch_result: &Json) -> Result<String, String> {
    let table = batch_result.at("table")?;
    let mut text = table.at("columns")?.items()?.iter().map(cell).collect::<Vec<_>>().join(",");
    text.push('\n');
    for row in table.at("rows")?.items()? {
        text += &row.items()?.iter().map(cell).collect::<Vec<_>>().join(",");
        text.push('\n');
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wildcards_and_paths() {
        assert!(matches(b"a*.xyz", b"abc.xyz"));
        assert!(matches(b"a?c", b"abc"));
        assert!(!matches(b"a?c", b"abbc"));
        let v = Json::parse(r#"{"a": {"b": [1, {"c": 2}]}}"#).unwrap();
        assert_eq!(value_at(&v, "a.b.1.c"), Json::Number(2.0));
        assert_eq!(value_at(&v, "a.x"), Json::Null);
        assert_eq!(cell(&Json::from("x,y")), "\"x,y\"");
        assert_eq!(cell(&Json::Number(0.1)), "0.1");
    }
}
