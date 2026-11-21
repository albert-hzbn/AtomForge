//! Triangle meshes from Wavefront OBJ and STL files (the C++
//! algorithms/MeshLoader), read as the C++ streams read them.

use glam::Vec3;

/// Vertices and triangle indices.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    pub vertices: Vec<Vec3>,
    pub indices: Vec<u32>,
}

fn is_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')
}

/// The first word of a line and the rest after it (`ls >> tag`).
fn tag_and_rest(line: &str) -> (&str, &str) {
    let trimmed = line.trim_start_matches(is_space);
    let end = trimmed.find(is_space).unwrap_or(trimmed.len());
    (&trimmed[..end], &trimmed[end..])
}

fn three_floats(rest: &str) -> Option<Vec3> {
    let v = af_core::format::read_f32_words(rest, 3)?;
    Some(Vec3::new(v[0], v[1], v[2]))
}

/// Text lines as std::getline gives them from a text-mode stream on Windows
/// (the C++ opens OBJ files in text mode: CRLF reads as LF).
fn text_lines(bytes: &[u8]) -> Vec<String> {
    let text = String::from_utf8_lossy(bytes).replace("\r\n", "\n");
    let mut lines: Vec<String> = text.split('\n').map(str::to_string).collect();
    if text.ends_with('\n') {
        lines.pop();
    }
    lines
}

pub fn parse_obj_mesh(path: &std::path::Path) -> Result<Mesh, String> {
    let bytes = std::fs::read(path).map_err(|_| "Cannot open OBJ file.".to_string())?;
    let mut mesh = Mesh::default();
    for line in text_lines(&bytes) {
        if line.is_empty() {
            continue;
        }
        let (tag, rest) = tag_and_rest(&line);
        if tag == "v" {
            if let Some(v) = three_floats(rest) {
                mesh.vertices.push(v);
            }
        } else if tag == "f" {
            let mut face = Vec::new();
            for token in rest.split(is_space).filter(|t| !t.is_empty()) {
                let id = &token[..token.find('/').unwrap_or(token.len())];
                if id.is_empty() {
                    continue;
                }
                let idx = af_io::text::atoi(id);
                if idx <= 0 {
                    continue;
                }
                face.push((idx - 1) as u32);
            }
            if face.len() >= 3 {
                for i in 1..face.len() - 1 {
                    mesh.indices.extend_from_slice(&[face[0], face[i], face[i + 1]]);
                }
            }
        }
    }
    if mesh.vertices.is_empty() || mesh.indices.is_empty() {
        return Err("OBJ has no usable vertices/faces.".into());
    }
    if mesh.indices.iter().any(|&i| i as usize >= mesh.vertices.len()) {
        return Err("OBJ face index out of range.".into());
    }
    Ok(mesh)
}

pub fn parse_stl_mesh(path: &std::path::Path) -> Result<Mesh, String> {
    let bytes = std::fs::read(path).map_err(|_| "Cannot open STL file.".to_string())?;
    let mut mesh = Mesh::default();
    let tri_count = if bytes.len() >= 84 { u32::from_le_bytes(bytes[80..84].try_into().unwrap()) } else { 0 };
    let looks_binary = bytes.len() as u64 == 84 + tri_count as u64 * 50;
    if looks_binary {
        for t in 0..tri_count as usize {
            let record = &bytes[84 + t * 50..84 + t * 50 + 50];
            let f = |k: usize| f32::from_le_bytes(record[k * 4..k * 4 + 4].try_into().unwrap());
            let base = mesh.vertices.len() as u32;
            mesh.vertices.push(Vec3::new(f(3), f(4), f(5)));
            mesh.vertices.push(Vec3::new(f(6), f(7), f(8)));
            mesh.vertices.push(Vec3::new(f(9), f(10), f(11)));
            mesh.indices.extend_from_slice(&[base, base + 1, base + 2]);
        }
        if mesh.vertices.is_empty() {
            return Err(String::new());
        }
        return Ok(mesh);
    }
    // ASCII: "vertex x y z" lines, three per triangle (binary-mode lines keep '\r').
    let text = String::from_utf8_lossy(&bytes);
    let mut tri: Vec<Vec3> = Vec::with_capacity(3);
    for line in text.split('\n') {
        let (tag, rest) = tag_and_rest(line);
        if tag == "vertex" {
            if let Some(v) = three_floats(rest) {
                tri.push(v);
                if tri.len() == 3 {
                    let base = mesh.vertices.len() as u32;
                    mesh.vertices.append(&mut tri);
                    mesh.indices.extend_from_slice(&[base, base + 1, base + 2]);
                }
            }
        }
    }
    if mesh.vertices.is_empty() || mesh.indices.is_empty() {
        return Err("ASCII STL has no usable triangles.".into());
    }
    Ok(mesh)
}
