//! Science-tool data input and output (the C++ science/ScienceData): frames
//! of structure and trajectory files with dynamics data, streamed trajectory
//! access, numeric data files (NPY, text tables, EIGENVAL, JSON fields), the
//! structure JSON form and extXYZ output.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use af_core::elements;
use af_core::format::format_g;
use af_core::json::Json;
use af_core::task::checkpoint;
use af_core::{AtomSite, Structure};

use crate::core::{cell_matrix, cell_volume, mat3_json, parse_pbc, points, row_times, scale, to_array, Mat3, NdArray, Pbc, Vec3};

/// ASE time unit (CODATA 2014): 1 fs in sqrt(amu*A^2/eV).
pub fn ase_femtosecond() -> f64 {
    1e-5 * (1.6021766208e-19f64 / 1.660539040e-27).sqrt()
}

/// Optional reader for formats beyond the built-in ones.
pub type StructureReader<'a> = Option<&'a dyn Fn(&Path) -> Result<Structure, String>>;

/// The periodic axes of a structure (none without a cell).
pub fn periodic_axes(structure: &Structure) -> Pbc {
    let c = structure.has_unit_cell;
    [c && structure.pbc[0], c && structure.pbc[1], c && structure.pbc[2]]
}

/// One frame of a structure or trajectory file.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FrameData {
    pub structure: Structure,
    pub pbc: Pbc,
    /// Angstrom/fs; empty when absent.
    pub velocities: Vec<Vec3>,
    /// amu; empty when absent.
    pub masses: Vec<f64>,
    pub time_fs: Option<f64>,
}

fn frame_structure() -> Structure {
    Structure::default()
}

fn upper(text: &str) -> String {
    text.to_ascii_uppercase()
}
fn lower(text: &str) -> String {
    text.to_ascii_lowercase()
}

pub fn display(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn read_file(path: &Path) -> Result<Vec<u8>, String> {
    std::fs::read(path).map_err(|_| format!("Cannot open {}", display(path)))
}

fn split(line: &str) -> Vec<&str> {
    af_io::text::words(line)
}

fn blank(line: &str) -> bool {
    line.chars().all(|c| c == ' ' || c == '\t' || c == '\r')
}

/// std::stod over the whole token, finite.
pub fn parse_double(token: &str, what: &str) -> Result<f64, String> {
    match af_io::text::stod_used(token) {
        Ok((value, used)) if used == token.len() && value.is_finite() => Ok(value),
        _ => Err(format!("Invalid number in {what}: {token}")),
    }
}

pub fn atomic_number(symbol: &str) -> i32 {
    elements::atomic_number(symbol)
}

pub fn atomic_mass(symbol: &str) -> Result<f64, String> {
    let z = atomic_number(symbol);
    if z == 0 {
        return Err(format!("Unknown element: {symbol}"));
    }
    Ok(elements::atomic_mass(z))
}

fn make_atom(symbol: &str, position: Vec3) -> Result<AtomSite, String> {
    if atomic_number(symbol) == 0 {
        return Err(format!("Unknown element: {symbol}"));
    }
    Ok(AtomSite::new(symbol, position))
}

/// A line reader over text (lines without the newline; '\r' kept as the C++
/// binary streams keep it, which every parser treats as whitespace).
pub struct LineReader<'a> {
    text: &'a str,
    pos: usize,
}

impl<'a> LineReader<'a> {
    pub fn new(text: &'a str) -> Self {
        Self { text, pos: 0 }
    }
    pub fn at(text: &'a str, pos: usize) -> Self {
        Self { text, pos: pos.min(text.len()) }
    }
    pub fn tell(&self) -> usize {
        self.pos
    }
    pub fn seek(&mut self, pos: usize) {
        self.pos = pos.min(self.text.len());
    }
    pub fn next_line(&mut self) -> Option<&'a str> {
        if self.pos >= self.text.len() {
            return None;
        }
        let rest = &self.text[self.pos..];
        match rest.find('\n') {
            Some(n) => {
                self.pos += n + 1;
                Some(&rest[..n])
            }
            None => {
                self.pos = self.text.len();
                Some(rest)
            }
        }
    }
}

/// extXYZ comment-line key=value pairs, honouring quoted values.
fn parse_comment(line: &str) -> BTreeMap<String, String> {
    let bytes = line.as_bytes();
    let space = |c: u8| c.is_ascii_whitespace() || c == 0x0b;
    let mut result = BTreeMap::new();
    let mut i = 0;
    while i < bytes.len() {
        while i < bytes.len() && space(bytes[i]) {
            i += 1;
        }
        let key_start = i;
        while i < bytes.len() && bytes[i] != b'=' && !space(bytes[i]) {
            i += 1;
        }
        let key = &line[key_start..i];
        if key.is_empty() {
            i += 1;
            continue;
        }
        if i >= bytes.len() || bytes[i] != b'=' {
            result.insert(lower(key), "T".to_string());
            continue;
        }
        i += 1;
        let value;
        if i < bytes.len() && (bytes[i] == b'"' || bytes[i] == b'{') {
            let close = if bytes[i] == b'"' { '"' } else { '}' };
            let end = line[i + 1..].find(close).map(|e| e + i + 1);
            value = line[i + 1..end.unwrap_or(line.len())].to_string();
            i = end.map_or(line.len(), |e| e + 1);
        } else {
            let start = i;
            while i < bytes.len() && !space(bytes[i]) {
                i += 1;
            }
            value = line[start..i].to_string();
        }
        result.insert(lower(key), value);
    }
    result
}

/// Up to `max_frames` extended-XYZ frames from the reader position.
pub fn read_extxyz_frames(input: &mut LineReader, max_frames: usize) -> Result<Vec<FrameData>, String> {
    let mut frames = Vec::new();
    let mut total = 0usize;
    while frames.len() < max_frames {
        let Some(line) = input.next_line() else { break };
        checkpoint()?;
        if blank(line) {
            continue;
        }
        let count_tokens = split(line);
        if count_tokens.len() != 1 {
            return Err("Invalid XYZ frame atom count".into());
        }
        let count = parse_double(count_tokens[0], "XYZ atom count")? as i64;
        if !(0..=10_000_000).contains(&count) {
            return Err("Invalid XYZ frame atom count".into());
        }
        total += count as usize;
        if total > 50_000_000 {
            return Err("Trajectory exceeds memory limit".into());
        }
        let line = input.next_line().ok_or("Missing XYZ frame comment")?;
        let info = parse_comment(line);
        let mut frame = FrameData { structure: frame_structure(), ..Default::default() };
        if let Some(lattice) = info.get("lattice") {
            let tokens = split(lattice);
            if tokens.len() != 9 {
                return Err("Invalid XYZ lattice".into());
            }
            for k in 0..9 {
                frame.structure.cell_vectors[k / 3][k % 3] = parse_double(tokens[k], "XYZ lattice")?;
            }
            frame.structure.has_unit_cell = true;
            frame.pbc = [true; 3];
        }
        if let Some(pbc) = info.get("pbc") {
            let tokens = split(pbc);
            if tokens.len() == 3 {
                for k in 0..3 {
                    frame.pbc[k] = upper(tokens[k]).starts_with('T');
                }
            }
        }
        if let Some(time) = info.get("time_fs") {
            frame.time_fs = Some(parse_double(time, "time_fs")?);
        }
        struct Column {
            name: String,
            kind: char,
            count: usize,
        }
        let mut columns = vec![
            Column { name: "species".into(), kind: 'S', count: 1 },
            Column { name: "pos".into(), kind: 'R', count: 3 },
        ];
        if let Some(properties) = info.get("properties") {
            columns.clear();
            // std::getline on ':' drops a trailing empty part.
            let mut parts: Vec<&str> = properties.split(':').collect();
            if parts.last() == Some(&"") {
                parts.pop();
            }
            if parts.len() % 3 != 0 {
                return Err("Invalid extXYZ Properties".into());
            }
            for p in (0..parts.len()).step_by(3) {
                let kind = parts[p + 1].chars().next().unwrap_or('\0').to_ascii_uppercase();
                columns.push(Column { name: lower(parts[p]), kind, count: parse_double(parts[p + 2], "Properties")? as i64 as usize });
            }
        }
        let mut momenta = Vec::new();
        let (mut has_momenta, mut has_velocities, mut has_masses) = (false, false, false);
        for _ in 0..count {
            let line = input.next_line().ok_or("Truncated XYZ trajectory")?;
            let tokens = split(line);
            let mut cursor = 0usize;
            let mut symbol = String::new();
            let (mut position, mut velocity, mut momentum) = ([0.0; 3], [0.0; 3], [0.0; 3]);
            let mut mass = 0.0;
            for column in &columns {
                if cursor + column.count > tokens.len() {
                    return Err("Invalid XYZ atom row".into());
                }
                let name = column.name.as_str();
                if name == "species" {
                    symbol = tokens[cursor].to_string();
                } else if name == "z" && column.kind == 'I' && symbol.is_empty() {
                    symbol = elements::symbol(parse_double(tokens[cursor], "atomic number")? as i32).to_string();
                } else if (name == "pos" || name == "positions") && column.count == 3 {
                    for k in 0..3 {
                        position[k] = parse_double(tokens[cursor + k], "XYZ positions")?;
                    }
                } else if name == "momenta" && column.count == 3 {
                    has_momenta = true;
                    for k in 0..3 {
                        momentum[k] = parse_double(tokens[cursor + k], "XYZ momenta")?;
                    }
                } else if (name == "velocities" || name == "vel") && column.count == 3 {
                    has_velocities = true;
                    for k in 0..3 {
                        velocity[k] = parse_double(tokens[cursor + k], "XYZ velocities")?;
                    }
                } else if name == "masses" && column.count == 1 {
                    has_masses = true;
                    mass = parse_double(tokens[cursor], "XYZ masses")?;
                } else if name == "property" && column.count == 1 {
                    frame.structure.atom_property.push(parse_double(tokens[cursor], "XYZ property")?);
                }
                cursor += column.count;
            }
            if symbol.is_empty() {
                return Err("extXYZ frame has no species column".into());
            }
            frame.structure.atoms.push(make_atom(&symbol, position)?);
            frame.velocities.push(velocity);
            momenta.push(momentum);
            frame.masses.push(if has_masses { mass } else { atomic_mass(&symbol)? });
        }
        if has_momenta && !has_velocities {
            for a in 0..momenta.len() {
                frame.velocities[a] = scale(momenta[a], ase_femtosecond() / frame.masses[a]);
            }
        }
        if !has_momenta && !has_velocities {
            frame.velocities.clear();
        }
        if !has_masses {
            frame.masses.clear();
        }
        frames.push(frame);
    }
    if frames.is_empty() {
        return Err("Structure/trajectory file contains no frames".into());
    }
    Ok(frames)
}

/// VASP 5 POSCAR/XDATCAR blocks: a header sets the cell and species.
#[derive(Default)]
pub struct VaspReader {
    pub cell: Mat3,
    pub symbols: Vec<String>,
}

impl VaspReader {
    /// false at end of file (when reading the comment).
    pub fn read_header(&mut self, input: &mut LineReader, skip_comment: bool) -> Result<bool, String> {
        if !skip_comment && input.next_line().is_none() {
            return Ok(false);
        }
        let line = input.next_line().ok_or("Truncated VASP header")?;
        let scale_tokens = split(line);
        if scale_tokens.is_empty() {
            return Err("Missing VASP scale factor".into());
        }
        let mut factor = parse_double(scale_tokens[0], "VASP scale")?;
        for row in 0..3 {
            let line = input.next_line().ok_or("Truncated VASP lattice")?;
            let tokens = split(line);
            if tokens.len() < 3 {
                return Err("Invalid VASP lattice".into());
            }
            for k in 0..3 {
                self.cell[row][k] = parse_double(tokens[k], "VASP lattice")?;
            }
        }
        if factor < 0.0 {
            factor = (-factor / cell_volume(&self.cell)).cbrt();
        }
        for row in self.cell.iter_mut() {
            *row = scale(*row, factor);
        }
        let line = input.next_line().ok_or("Truncated VASP species")?;
        let names = split(line);
        if names.is_empty() || names[0].as_bytes()[0].is_ascii_digit() {
            return Err("VASP files must list element symbols (VASP 5 format)".into());
        }
        let line = input.next_line().ok_or("Truncated VASP counts")?;
        let counts = split(line);
        if counts.len() != names.len() {
            return Err("VASP species and counts differ".into());
        }
        self.symbols.clear();
        for (s, name) in names.iter().enumerate() {
            let name = &name[..name.find(['/', '_']).unwrap_or(name.len())];
            let count = parse_double(counts[s], "VASP counts")? as i32;
            for _ in 0..count.max(0) {
                self.symbols.push(name.to_string());
            }
        }
        Ok(true)
    }

    pub fn read_positions(&self, input: &mut LineReader, direct: bool) -> Result<FrameData, String> {
        let mut frame = FrameData { pbc: [true; 3], ..Default::default() };
        frame.structure.has_unit_cell = true;
        frame.structure.cell_vectors = self.cell;
        for symbol in &self.symbols {
            let line = input.next_line().ok_or("Truncated VASP coordinates")?;
            let tokens = split(line);
            if tokens.len() < 3 {
                return Err("Invalid VASP coordinates".into());
            }
            let value = [
                parse_double(tokens[0], "VASP coordinates")?,
                parse_double(tokens[1], "VASP coordinates")?,
                parse_double(tokens[2], "VASP coordinates")?,
            ];
            frame.structure.atoms.push(make_atom(symbol, if direct { row_times(value, &self.cell) } else { value })?);
        }
        Ok(frame)
    }
}

fn text_of(bytes: Vec<u8>) -> String {
    String::from_utf8(bytes).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned())
}

fn first_non_blank(line: &str) -> Option<usize> {
    line.find(|c| c != ' ' && c != '\t')
}

fn read_vasp(path: &Path, trajectory: bool) -> Result<Vec<FrameData>, String> {
    let text = text_of(read_file(path)?);
    let mut input = LineReader::new(&text);
    let mut vasp = VaspReader::default();
    let mut frames = Vec::new();
    if !vasp.read_header(&mut input, false)? {
        return Err("Empty VASP file".into());
    }
    if !trajectory {
        let mut line = input.next_line().ok_or("Truncated POSCAR")?;
        line = line.strip_suffix('\r').unwrap_or(line);
        if line.starts_with(['S', 's']) {
            line = input.next_line().ok_or("Truncated POSCAR")?;
            line = line.strip_suffix('\r').unwrap_or(line);
        }
        let mode = match first_non_blank(line) {
            None => 'D',
            Some(i) => (line.as_bytes()[i] as char).to_ascii_uppercase(),
        };
        frames.push(vasp.read_positions(&mut input, mode != 'C' && mode != 'K')?);
        return Ok(frames);
    }
    while let Some(line) = input.next_line() {
        checkpoint()?;
        if blank(line) {
            continue;
        }
        let head = lower(&line[first_non_blank(line).unwrap_or(0)..]);
        if head.starts_with("direct") || head.starts_with("cartesian") {
            frames.push(vasp.read_positions(&mut input, head.starts_with('d'))?);
        } else {
            vasp.read_header(&mut input, true)?;
        }
    }
    if frames.is_empty() {
        return Err("XDATCAR contains no configurations".into());
    }
    Ok(frames)
}

/// Structure from {"symbols", "positions", "cell", "pbc", per-atom data}.
pub fn structure_from_json(value: &Json) -> Result<Structure, String> {
    let Json::Object(members) = value else {
        return Err("Structure requires a file reference or symbols, positions and optional cell".into());
    };
    const KNOWN: [&str; 8] = ["symbols", "positions", "cell", "pbc", "magmoms", "charges", "velocities", "fixed"];
    if members.iter().any(|(k, _)| !KNOWN.contains(&k.as_str())) {
        return Err("Structure requires a file reference or symbols, positions, optional cell, pbc, magmoms, charges, velocities and fixed".into());
    }
    let mut structure = Structure::default();
    let symbols = value.at("symbols")?.items()?;
    let positions = points(to_array(value.at("positions")?, "positions")?, "positions")?;
    if symbols.len() != positions.len() {
        return Err("Structure symbols and positions differ in length".into());
    }
    for (symbol, position) in symbols.iter().zip(&positions) {
        structure.atoms.push(make_atom(symbol.string()?, *position)?);
    }
    if let Some(cell) = value.find("cell").filter(|c| !c.is_null()) {
        structure.cell_vectors = cell_matrix(to_array(cell, "cell")?)?;
        structure.has_unit_cell = true;
    }
    if let Some(pbc) = value.find("pbc").filter(|c| !c.is_null()) {
        structure.pbc = parse_pbc(pbc)?;
    }
    let count = structure.atoms.len();
    let per_atom = |key: &str| -> Result<Option<&Vec<Json>>, String> {
        match value.find(key) {
            None | Some(Json::Null) => Ok(None),
            Some(Json::Array(items)) if items.len() == count => Ok(Some(items)),
            Some(_) => Err(format!("Structure {key} needs one entry per atom")),
        }
    };
    let component = |item: &Json, k: usize| -> Result<Json, String> {
        item.items()?.get(k).cloned().ok_or_else(|| "vector::_M_range_check".to_string())
    };
    if let Some(list) = per_atom("magmoms")? {
        for (atom, item) in structure.atoms.iter_mut().zip(list) {
            atom.magmom = item.number()?;
        }
    }
    if let Some(list) = per_atom("charges")? {
        for (atom, item) in structure.atoms.iter_mut().zip(list) {
            atom.charge = item.number()?;
        }
    }
    if let Some(list) = per_atom("velocities")? {
        for (atom, item) in structure.atoms.iter_mut().zip(list) {
            for k in 0..3 {
                atom.velocity[k] = component(item, k)?.number()?;
            }
        }
    }
    if let Some(list) = per_atom("fixed")? {
        for (atom, item) in structure.atoms.iter_mut().zip(list) {
            for k in 0..3 {
                if component(item, k)?.boolean()? {
                    atom.fixed |= 1 << k;
                }
            }
        }
    }
    Ok(structure)
}

/// A .npy array (little-endian numeric; Fortran order reordered).
pub fn load_npy(path: &Path) -> Result<NdArray, String> {
    let data = read_file(path)?;
    if data.len() < 10 || &data[..6] != b"\x93NUMPY" {
        return Err("Invalid NPY file".into());
    }
    let major = data[6];
    let (header_length, offset) = if major == 1 {
        (data[8] as usize | (data[9] as usize) << 8, 10usize)
    } else {
        if data.len() < 12 {
            return Err("Invalid NPY file".into());
        }
        ((0..4).fold(0usize, |acc, b| acc | (data[8 + b] as usize) << (8 * b)), 12usize)
    };
    if offset + header_length > data.len() {
        return Err("Invalid NPY header".into());
    }
    let header = String::from_utf8_lossy(&data[offset..offset + header_length]).into_owned();
    let field = |key: &str| -> Result<String, String> {
        let position = header.find(&format!("'{key}'")).ok_or_else(|| format!("NPY header lacks {key}"))?;
        let colon = header[position..].find(':').map(|c| c + position + 1).unwrap_or(header.len());
        Ok(header[colon.min(header.len())..].to_string())
    };
    let mut descr = field("descr")?;
    descr = descr[descr.find('\'').map_or(0, |q| q + 1)..].to_string();
    descr = descr[..descr.find('\'').unwrap_or(descr.len())].to_string();
    let fortran_text = field("fortran_order")?;
    let npos = usize::MAX;
    let fortran = fortran_text.find("True").unwrap_or(npos) < fortran_text.find(',').unwrap_or(npos);
    let mut shape_text = field("shape")?;
    shape_text = shape_text[shape_text.find('(').map_or(0, |p| p + 1)..].to_string();
    shape_text = shape_text[..shape_text.find(')').unwrap_or(shape_text.len())].to_string();
    let mut result = NdArray::default();
    let mut parts: Vec<&str> = shape_text.split(',').collect();
    if parts.last() == Some(&"") {
        parts.pop();
    }
    for token in parts {
        if token.chars().any(|c| c != ' ') {
            let t = token.trim_start();
            let digits: String = t.chars().take_while(|c| c.is_ascii_digit()).collect();
            result.shape.push(digits.parse::<usize>().map_err(|_| "stoull".to_string())?);
        }
    }
    let count: usize = result.shape.iter().product();
    if descr.len() < 3 || descr.starts_with('>') {
        return Err("Only little-endian numeric NPY arrays are supported".into());
    }
    let kind = descr.as_bytes()[1];
    let width: usize = descr[2..].chars().take_while(|c| c.is_ascii_digit()).collect::<String>().parse().map_err(|_| "stoul".to_string())?;
    let start = offset + header_length;
    if start + count * width > data.len() {
        return Err("Truncated NPY data".into());
    }
    let mut values = vec![0.0; count];
    for (i, value) in values.iter_mut().enumerate() {
        let b = &data[start + i * width..start + (i + 1) * width];
        *value = match (kind, width) {
            (b'f', 8) => f64::from_le_bytes(b.try_into().unwrap()),
            (b'f', 4) => f32::from_le_bytes(b.try_into().unwrap()) as f64,
            (b'i', 8) => i64::from_le_bytes(b.try_into().unwrap()) as f64,
            (b'i', 4) => i32::from_le_bytes(b.try_into().unwrap()) as f64,
            (b'i', 2) => i16::from_le_bytes(b.try_into().unwrap()) as f64,
            (b'u', 1) => b[0] as f64,
            (b'b', 1) => {
                if b[0] != 0 {
                    1.0
                } else {
                    0.0
                }
            }
            _ => return Err(format!("Unsupported NPY data type {descr}")),
        };
    }
    if fortran && result.shape.len() > 1 {
        result.values = vec![0.0; count];
        let mut index = vec![0usize; result.shape.len()];
        for value in values {
            let (mut c, mut stride) = (0usize, 1usize);
            for axis in (0..result.shape.len()).rev() {
                c += index[axis] * stride;
                stride *= result.shape[axis];
            }
            result.values[c] = value;
            for axis in 0..index.len() {
                index[axis] += 1;
                if index[axis] < result.shape[axis] {
                    break;
                }
                index[axis] = 0;
            }
        }
    } else {
        result.values = values;
    }
    Ok(result)
}

/// A whitespace- (or comma-) separated numeric table; '#' starts a comment.
pub fn load_table(path: &Path, comma: bool) -> Result<NdArray, String> {
    let text = text_of(read_file(path)?);
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for line in af_io::text::lines(&text, true) {
        let mut line = line[..line.find('#').unwrap_or(line.len())].to_string();
        if comma {
            line = line.replace(',', " ");
        }
        let tokens = split(&line);
        if tokens.is_empty() {
            continue;
        }
        let row = tokens.iter().map(|t| parse_double(t, "numeric table")).collect::<Result<Vec<_>, _>>()?;
        if !rows.is_empty() && row.len() != rows[0].len() {
            return Err("Numeric table rows have different lengths".into());
        }
        rows.push(row);
    }
    if rows.is_empty() {
        return Err("Numeric table is empty".into());
    }
    let columns = rows[0].len();
    Ok(NdArray { shape: vec![rows.len(), columns], values: rows.concat() })
}

/// A whole text file, decompressing .gz (.bz2 and .xz need decompressing first).
pub fn read_text_file(path: &Path) -> Result<String, String> {
    let extension = path.extension().map(|e| format!(".{}", lower(&e.to_string_lossy()))).unwrap_or_default();
    if extension != ".gz" && extension != ".bz2" && extension != ".xz" {
        return Ok(text_of(read_file(path)?));
    }
    let packed = read_file(path)?;
    if extension == ".gz" {
        use std::io::Read;
        let mut text = Vec::new();
        flate2::read::MultiGzDecoder::new(&packed[..])
            .read_to_end(&mut text)
            .map_err(|_| format!("Invalid gzip file {}", display(path)))?;
        return Ok(text_of(text));
    }
    if extension == ".bz2" {
        use std::io::Read;
        let mut text = Vec::new();
        bzip2::read::MultiBzDecoder::new(&packed[..])
            .read_to_end(&mut text)
            .map_err(|_| format!("Invalid bzip2 file {}", display(path)))?;
        return Ok(text_of(text));
    }
    let mut text = Vec::new();
    lzma_rs::xz_decompress(&mut &packed[..], &mut text).map_err(|_| format!("Invalid xz file {}", display(path)))?;
    Ok(text_of(text))
}

/// VASP EIGENVAL -> (spin, k, band) energies.
pub fn load_eigenval(path: &Path) -> Result<NdArray, String> {
    let text = read_text_file(path)?;
    let mut input = LineReader::new(&text);
    let mut line = input.next_line().unwrap_or("");
    let first = split(line);
    if first.len() < 4 {
        return Err("Invalid EIGENVAL header".into());
    }
    let spins = parse_double(first[3], "EIGENVAL")? as usize;
    for _ in 0..5 {
        line = input.next_line().unwrap_or("");
    }
    let sizes = split(line);
    if sizes.len() < 3 || (spins != 1 && spins != 2) {
        return Err("Invalid EIGENVAL header".into());
    }
    let kpoints = parse_double(sizes[1], "EIGENVAL")? as usize;
    let bands = parse_double(sizes[2], "EIGENVAL")? as usize;
    let mut result = NdArray::new(vec![spins, kpoints, bands], 0.0);
    for k in 0..kpoints {
        loop {
            let line = input.next_line().ok_or("Truncated EIGENVAL")?;
            if !split(line).is_empty() {
                break;
            }
        }
        for b in 0..bands {
            let line = input.next_line().ok_or("Truncated EIGENVAL")?;
            let tokens = split(line);
            if tokens.len() < 1 + spins {
                return Err("Invalid EIGENVAL band row".into());
            }
            for s in 0..spins {
                result.values[(s * kpoints + k) * bands + b] = parse_double(tokens[1 + s], "EIGENVAL")?;
            }
        }
    }
    Ok(result)
}

/// A dotted JSON path ("a.0.b").
pub fn select_field(value: &Json, field: &str) -> Result<Json, String> {
    let mut value = value.clone();
    let mut parts: Vec<&str> = field.split('.').collect();
    if parts.last() == Some(&"") {
        parts.pop();
    }
    for component in parts {
        value = if value.is_array() {
            let index = af_io::text::stoll(component)?;
            if index < 0 || index as usize >= value.size() {
                return Err(format!("JSON field index out of range: {component}"));
            }
            value.items()?[index as usize].clone()
        } else {
            value.at(component)?.clone()
        };
    }
    Ok(value)
}

pub fn select_column(table: NdArray, column: Option<usize>) -> Result<NdArray, String> {
    let Some(column) = column else { return Ok(table) };
    if table.ndim() != 2 || column >= table.shape[1] {
        return Err("Column is outside the two-dimensional data table".into());
    }
    Ok(NdArray { shape: vec![table.shape[0]], values: (0..table.shape[0]).map(|i| table.at2(i, column)).collect() })
}

fn file_name(path: &Path) -> String {
    path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default()
}

pub fn extension(path: &Path) -> String {
    path.extension().map(|e| format!(".{}", lower(&e.to_string_lossy()))).unwrap_or_default()
}

fn is_vasp_structure(path: &Path) -> bool {
    let name = upper(&file_name(path));
    let ext = extension(path);
    name.contains("POSCAR") || name.contains("CONTCAR") || ext == ".vasp" || ext == ".poscar"
}

fn is_lammps_dump(path: &Path) -> bool {
    let ext = extension(path);
    if ext == ".lammpstrj" || ext == ".dump" {
        return true;
    }
    let Ok(bytes) = std::fs::read(path) else { return false };
    let text = String::from_utf8_lossy(&bytes);
    for line in text.split('\n') {
        if !blank(line) {
            return line.starts_with("ITEM:");
        }
    }
    false
}

fn read_extxyz(path: &Path) -> Result<Vec<FrameData>, String> {
    let text = text_of(read_file(path)?);
    read_extxyz_frames(&mut LineReader::new(&text), usize::MAX)
}

fn lammps_frame(structure: Structure) -> FrameData {
    let p = structure.has_unit_cell;
    FrameData { pbc: [p, p, p], structure, ..Default::default() }
}

/// All frames of a structure or trajectory file.
pub fn read_frames(path: &Path, reader: StructureReader) -> Result<Vec<FrameData>, String> {
    if !path.is_file() {
        return Err(format!("File not found: {}", display(path)));
    }
    let name = upper(&file_name(path));
    let ext = extension(path);
    if name.starts_with("XDATCAR") {
        return read_vasp(path, true);
    }
    if is_vasp_structure(path) {
        return read_vasp(path, false);
    }
    if ext == ".xyz" || ext == ".extxyz" {
        return read_extxyz(path);
    }
    if ext == ".json" {
        let structure = structure_from_json(&Json::parse(&text_of(read_file(path)?))?)?;
        let p = structure.has_unit_cell;
        return Ok(vec![FrameData { pbc: if p { [true; 3] } else { [false; 3] }, structure, ..Default::default() }]);
    }
    if is_lammps_dump(path) {
        return Ok(af_io::trajectory::load_lammps_trajectory(&display(path))?.into_iter().map(lammps_frame).collect());
    }
    let Some(reader) = reader else {
        return Err(format!("Unsupported structure file: {}", display(path)));
    };
    let structure = reader(path)?;
    if structure.atoms.is_empty() {
        return Err("Structure file contains no atoms".into());
    }
    let p = structure.has_unit_cell;
    Ok(vec![FrameData { pbc: if p { [true; 3] } else { [false; 3] }, structure, ..Default::default() }])
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum StreamFormat {
    Extxyz,
    Lammps,
    Xdatcar,
    Kept,
}

#[derive(Clone, Copy, Debug)]
struct Entry {
    start: usize,
    header: usize,
    direct: bool,
}

/// Random access to the frames of a trajectory file: one indexing pass
/// records where each frame starts; frame(i) parses only that frame.
pub struct TrajectoryStream {
    path: PathBuf,
    format: StreamFormat,
    entries: Vec<Entry>,
    kept: Vec<FrameData>,
}

impl TrajectoryStream {
    pub fn new(path: &Path, reader: StructureReader) -> Result<Self, String> {
        if !path.is_file() {
            return Err(format!("File not found: {}", display(path)));
        }
        let name = upper(&file_name(path));
        let ext = extension(path);
        let mut stream = Self { path: path.to_path_buf(), format: StreamFormat::Kept, entries: Vec::new(), kept: Vec::new() };
        let truncated = || format!("Truncated trajectory frame in {}", display(path));
        if name.starts_with("XDATCAR") {
            stream.format = StreamFormat::Xdatcar;
            let text = text_of(read_file(path)?);
            let mut input = LineReader::new(&text);
            let mut vasp = VaspReader::default();
            let mut header = 0usize;
            if !vasp.read_header(&mut input, false)? {
                return Err("Empty VASP file".into());
            }
            loop {
                let offset = input.tell();
                let Some(line) = input.next_line() else { break };
                checkpoint()?;
                if blank(line) {
                    continue;
                }
                let head = lower(&line[first_non_blank(line).unwrap_or(0)..]);
                if head.starts_with("direct") || head.starts_with("cartesian") {
                    stream.entries.push(Entry { start: input.tell(), header, direct: head.starts_with('d') });
                    for _ in 0..vasp.symbols.len() {
                        input.next_line().ok_or_else(truncated)?;
                    }
                } else {
                    header = offset;
                    vasp.read_header(&mut input, true)?;
                }
            }
        } else if ext == ".xyz" || ext == ".extxyz" {
            stream.format = StreamFormat::Extxyz;
            let text = text_of(read_file(path)?);
            let mut input = LineReader::new(&text);
            loop {
                let offset = input.tell();
                let Some(line) = input.next_line() else { break };
                checkpoint()?;
                if blank(line) {
                    continue;
                }
                let tokens = split(line);
                if tokens.len() != 1 {
                    return Err("Invalid XYZ frame atom count".into());
                }
                let count = parse_double(tokens[0], "XYZ atom count")?;
                if !(0.0..=10_000_000.0).contains(&count) {
                    return Err("Invalid XYZ frame atom count".into());
                }
                stream.entries.push(Entry { start: offset, header: 0, direct: true });
                for _ in 0..(count as i64 + 1) {
                    input.next_line().ok_or_else(truncated)?;
                }
            }
        } else if !is_vasp_structure(path) && ext != ".json" && is_lammps_dump(path) {
            stream.format = StreamFormat::Lammps;
            let text = text_of(read_file(path)?);
            let mut input = LineReader::new(&text);
            loop {
                let offset = input.tell();
                let Some(line) = input.next_line() else { break };
                checkpoint()?;
                if line.starts_with("ITEM: TIMESTEP") {
                    stream.entries.push(Entry { start: offset, header: 0, direct: true });
                }
            }
        } else {
            stream.kept = read_frames(path, reader)?;
            return Ok(stream);
        }
        if stream.entries.is_empty() {
            return Err("Structure/trajectory file contains no frames".into());
        }
        Ok(stream)
    }

    pub fn len(&self) -> usize {
        if self.format == StreamFormat::Kept {
            self.kept.len()
        } else {
            self.entries.len()
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn streamed(&self) -> bool {
        self.kept.is_empty()
    }

    pub fn frame(&self, index: usize) -> Result<FrameData, String> {
        if index >= self.len() {
            return Err(format!("Trajectory frame {index} of {}", self.len()));
        }
        if self.format == StreamFormat::Kept {
            return Ok(self.kept[index].clone());
        }
        let text = text_of(read_file(&self.path)?);
        let entry = self.entries[index];
        let mut input = LineReader::at(&text, entry.header);
        match self.format {
            StreamFormat::Xdatcar => {
                let mut vasp = VaspReader::default();
                vasp.read_header(&mut input, false)?;
                input.seek(entry.start);
                vasp.read_positions(&mut input, entry.direct)
            }
            StreamFormat::Extxyz => {
                input.seek(entry.start);
                Ok(read_extxyz_frames(&mut input, 1)?.remove(0))
            }
            _ => {
                input.seek(entry.start);
                let mut lines = std::iter::from_fn(|| input.next_line());
                let structure = af_io::trajectory::load_lammps_frames(&mut lines, 1)?.remove(0);
                Ok(lammps_frame(structure))
            }
        }
    }
}

/// positions/velocities/cell/masses of all frames as an array.
pub fn frames_field(frames: &[FrameData], field: &str, timestep: Option<f64>) -> Result<NdArray, String> {
    if frames.is_empty() {
        return Err("Structure/trajectory file contains no frames".into());
    }
    if let Some(timestep) = timestep.filter(|_| field == "positions" || field == "velocities") {
        let any = frames.iter().any(|f| f.time_fs.is_some());
        let all = frames.iter().all(|f| f.time_fs.is_some());
        if any {
            if !all {
                return Err("Trajectory has incomplete time_fs metadata".into());
            }
            for f in 1..frames.len() {
                let interval = frames[f].time_fs.unwrap() - frames[f - 1].time_fs.unwrap();
                if !(interval > 0.0) || (interval - timestep).abs() > 1e-10 + 1e-8 * timestep.abs() {
                    return Err("Trajectory sampling must be uniform and match timestep_fs; resample or trim the final partial interval".into());
                }
            }
        }
    }
    let first = &frames[0].structure.atoms;
    for frame in frames {
        let atoms = &frame.structure.atoms;
        if atoms.len() != first.len() || atoms.iter().zip(first).any(|(a, b)| a.symbol != b.symbol) {
            return Err("Trajectory atom species/order must remain consistent".into());
        }
    }
    let atoms = first.len();
    match field {
        "velocities" => {
            let mut result = NdArray::new(vec![frames.len(), atoms, 3], 0.0);
            for (f, frame) in frames.iter().enumerate() {
                if frame.velocities.len() != atoms {
                    return Err("Trajectory has no velocities; supply explicit A/fs data".into());
                }
                for a in 0..atoms {
                    for k in 0..3 {
                        result.values[(f * atoms + a) * 3 + k] = frame.velocities[a][k];
                    }
                }
            }
            Ok(result)
        }
        "cell" => {
            if !frames[0].structure.has_unit_cell {
                return Err("Structure file has no cell".into());
            }
            let reference = frames[0].structure.cell_vectors;
            for frame in frames {
                for r in 0..3 {
                    for c in 0..3 {
                        if !frame.structure.has_unit_cell || (frame.structure.cell_vectors[r][c] - reference[r][c]).abs() > 1e-10 {
                            return Err("This analysis requires a fixed cell".into());
                        }
                    }
                }
            }
            Ok(NdArray { shape: vec![3, 3], values: reference.iter().flatten().copied().collect() })
        }
        "masses" => {
            let values = (0..atoms)
                .map(|a| if frames[0].masses.len() == atoms { Ok(frames[0].masses[a]) } else { atomic_mass(&first[a].symbol) })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(NdArray { shape: vec![atoms], values })
        }
        "positions" => {
            let mut result = NdArray::new(vec![frames.len(), atoms, 3], 0.0);
            for (f, frame) in frames.iter().enumerate() {
                for (a, atom) in frame.structure.atoms.iter().enumerate() {
                    let base = (f * atoms + a) * 3;
                    result.values[base..base + 3].copy_from_slice(&[atom.x, atom.y, atom.z]);
                }
            }
            if frames.len() == 1 {
                result.shape.remove(0);
            }
            Ok(result)
        }
        _ => Err("Structure field must be positions, velocities, cell or masses".into()),
    }
}

/// A structure input together with its periodicity.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StructureInput {
    pub structure: Structure,
    pub pbc: Pbc,
}

/// The LAMMPS types of one file reference, restored afterwards.
pub struct ScopedLammpsTypes {
    previous: String,
    active: bool,
}

impl ScopedLammpsTypes {
    pub fn new(reference: &Json) -> Result<Self, String> {
        let previous = af_io::trajectory::lammps_type_elements();
        let mut active = false;
        if let Some(types) = reference.find("types") {
            af_io::trajectory::set_lammps_type_elements(types.string()?)?;
            active = true;
        }
        Ok(Self { previous, active })
    }
}

impl Drop for ScopedLammpsTypes {
    fn drop(&mut self) {
        if self.active {
            let _ = af_io::trajectory::set_lammps_type_elements(&self.previous);
        }
    }
}

/// Lexically normal base/relative (as std::filesystem::lexically_normal).
pub fn join_normal(base: &Path, relative: &str) -> PathBuf {
    let joined = base.join(relative);
    let mut parts: Vec<std::path::Component> = Vec::new();
    for component in joined.components() {
        match component {
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                if matches!(parts.last(), Some(std::path::Component::Normal(_))) {
                    parts.pop();
                } else if !matches!(parts.last(), Some(std::path::Component::RootDir | std::path::Component::Prefix(_))) {
                    parts.push(component);
                }
            }
            other => parts.push(other),
        }
    }
    let result: PathBuf = parts.iter().collect();
    if result.as_os_str().is_empty() {
        PathBuf::from(".")
    } else {
        result
    }
}

/// A structure parameter: a file reference ({"file", optional "types"}) or
/// the JSON form.
pub fn resolve_structure(value: &Json, base: &Path, reader: StructureReader) -> Result<StructureInput, String> {
    if value.is_object() && value.contains("file") && (value.size() == 1 || (value.size() == 2 && value.contains("types"))) {
        let _types = ScopedLammpsTypes::new(value)?;
        let path = join_normal(base, value.at("file")?.string()?);
        if extension(&path) == ".json" {
            let parsed = Json::parse(&text_of(read_file(&path)?))?;
            return resolve_structure(&parsed, path.parent().unwrap_or(Path::new("")), reader);
        }
        let mut frames = read_frames(&path, reader)?;
        let last = frames.pop().expect("read_frames returns frames");
        return Ok(StructureInput { structure: last.structure, pbc: last.pbc });
    }
    let structure = structure_from_json(value)?;
    let p = structure.has_unit_cell;
    Ok(StructureInput { structure, pbc: [p, p, p] })
}

/// Writes frames as extXYZ (ASE-compatible momenta and time_fs metadata).
pub fn write_extxyz(path: &Path, frames: &[Structure], velocities: &[Vec<Vec3>], times: &[f64]) -> Result<(), String> {
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent).map_err(|_| format!("Cannot write {}", display(path)))?;
    }
    let mut output = String::new();
    let g = |v: f64| format_g(v, 12);
    for (f, frame) in frames.iter().enumerate() {
        let with_momenta = f < velocities.len() && velocities[f].len() == frame.atoms.len();
        output.push_str(&format!("{}\n", frame.atoms.len()));
        if frame.has_unit_cell {
            output.push_str("Lattice=\"");
            for k in 0..9 {
                if k > 0 {
                    output.push(' ');
                }
                output.push_str(&g(frame.cell_vectors[k / 3][k % 3]));
            }
            output.push_str("\" ");
        }
        output.push_str("Properties=species:S:1:pos:R:3");
        if with_momenta {
            output.push_str(":momenta:R:3");
        }
        if f < times.len() {
            output.push_str(&format!(" time_fs={}", g(times[f])));
        }
        output.push_str(&format!(" pbc=\"{}\"\n", if frame.has_unit_cell { "T T T" } else { "F F F" }));
        for (a, atom) in frame.atoms.iter().enumerate() {
            output.push_str(&format!("{} {} {} {}", atom.symbol, g(atom.x), g(atom.y), g(atom.z)));
            if with_momenta {
                let mass = atomic_mass(&atom.symbol)?;
                for k in 0..3 {
                    output.push_str(&format!(" {}", g(velocities[f][a][k] * mass / ase_femtosecond())));
                }
            }
            output.push('\n');
        }
    }
    // The C++ writes in text mode (CRLF on Windows).
    #[cfg(windows)]
    let output = output.replace('\n', "\r\n");
    std::fs::write(path, output).map_err(|_| format!("Cannot write {}", display(path)))
}

/// The JSON form of a structure.
pub fn structure_json(structure: &Structure) -> Json {
    let mut result = Json::object();
    let mut symbols = Json::array(Vec::new());
    let mut positions = Json::array(Vec::new());
    for atom in &structure.atoms {
        symbols.push(Json::String(atom.symbol.clone()));
        positions.push(Json::Array(vec![Json::Number(atom.x), Json::Number(atom.y), Json::Number(atom.z)]));
    }
    result.set("symbols", symbols);
    result.set("positions", positions);
    if structure.has_unit_cell {
        result.set("cell", mat3_json(&structure.cell_vectors));
        result.set("pbc", Json::Array(structure.pbc.iter().map(|&b| Json::Bool(b)).collect()));
    } else {
        result.set("cell", Json::Null);
    }
    let atoms = &structure.atoms;
    if atoms.iter().any(|a| a.magmom != 0.0) {
        result.set("magmoms", Json::Array(atoms.iter().map(|a| Json::Number(a.magmom)).collect()));
    }
    if atoms.iter().any(|a| a.charge != 0.0) {
        result.set("charges", Json::Array(atoms.iter().map(|a| Json::Number(a.charge)).collect()));
    }
    if atoms.iter().any(|a| a.velocity.iter().any(|&v| v != 0.0)) {
        result.set("velocities", Json::Array(atoms.iter().map(|a| Json::Array(a.velocity.iter().map(|&v| Json::Number(v)).collect())).collect()));
    }
    if atoms.iter().any(|a| a.fixed != 0) {
        result.set(
            "fixed",
            Json::Array(atoms.iter().map(|a| Json::Array((0..3).map(|k| Json::Bool(a.fixed & (1 << k) != 0)).collect())).collect()),
        );
    }
    result
}

/// Structure inputs are the catalog parameters of kind "structure" (without
/// a catalog entry: structure, initial and final).
fn is_structure_parameter(tool: &str, name: &str) -> bool {
    match crate::catalog::find_tool(tool) {
        Some(definition) => definition.parameters.iter().find(|p| p.name == name).is_some_and(|p| p.kind == "structure"),
        None => name == "structure" || name == "initial" || name == "final",
    }
}

/// Resolved tool parameters. Data references ({"file": ...}) are loaded when
/// constructed, relative to the base directory; files are never evaluated.
#[derive(Clone, Debug, Default)]
pub struct Parameters {
    json: BTreeMap<String, Json>,
    arrays: BTreeMap<String, NdArray>,
    structures: BTreeMap<String, StructureInput>,
    frames: BTreeMap<String, Vec<FrameData>>,
    cell_from_structure_file: bool,
}

impl Parameters {
    pub fn new(tool: &str, request: &Json, base: &Path, reader: StructureReader) -> Result<Self, String> {
        let Json::Object(members) = request else { return Err("parameters must be an object".into()) };
        let mut p = Parameters::default();
        let mut timestep = None;
        if tool == "msd" || tool == "vacf" || tool == "vibrational-spectrum" {
            if let Some(Json::Number(value)) = request.find("timestep_fs") {
                timestep = Some(*value);
            }
        }
        for (name, value) in members {
            checkpoint()?;
            let name = name.as_str();
            if is_structure_parameter(tool, name) {
                p.structures.insert(name.into(), resolve_structure(value, base, reader)?);
                continue;
            }
            if name == "restart_images" {
                if !value.is_object() || !value.contains("file") {
                    return Err("restart_images must reference a structure/trajectory file".into());
                }
                p.frames.insert(name.into(), read_frames(&join_normal(base, value.at("file")?.string()?), reader)?);
                continue;
            }
            if name.len() > 5 && name.ends_with("_file") {
                let path = if let Json::String(s) = value { s.as_str() } else { value.at("file")?.string()? };
                p.json.insert(name.into(), Json::String(display(&join_normal(base, path))));
                continue;
            }
            if !value.is_object() || !value.contains("file") || name == "calculator" || name == "calculator_factory" {
                p.json.insert(name.into(), value.clone());
                continue;
            }
            for (key, _) in value.members()? {
                if !["file", "field", "column", "types"].contains(&key.as_str()) {
                    return Err(format!("Unknown data-reference options for {name}"));
                }
            }
            let _types = ScopedLammpsTypes::new(value)?;
            let path = join_normal(base, value.at("file")?.string()?);
            let field = match value.find("field") {
                Some(f) if !f.is_null() => Some(f.string()?.to_string()),
                _ => None,
            };
            let column = match value.find("column") {
                Some(c) if !c.is_null() => match c {
                    Json::Number(n) if *n >= 0.0 && *n == n.floor() => Some(*n as usize),
                    _ => return Err("Column must be a nonnegative integer".into()),
                },
                _ => None,
            };
            let ext = extension(&path);
            let filename = upper(&file_name(&path));
            if !path.is_file() {
                return Err(format!("File not found: {}", display(&path)));
            }
            if ext == ".json" {
                let mut data = Json::parse(&text_of(read_file(&path)?))?;
                if let Some(field) = &field {
                    data = select_field(&data, field)?;
                }
                if field.is_none() && data.is_object() && data.contains("positions") && data.contains("symbols") {
                    data = data.at("positions")?.clone();
                }
                if column.is_some() {
                    p.arrays.insert(name.into(), select_column(to_array(&data, name)?, column)?);
                } else {
                    p.json.insert(name.into(), data);
                }
            } else if filename.starts_with("EIGENVAL") && name == "energies_eV" {
                if field.is_some() || column.is_some() {
                    return Err("EIGENVAL supplies all spin/k-point/band energies; field/column selection is unsupported".into());
                }
                p.arrays.insert(name.into(), load_eigenval(&path)?);
            } else if ext == ".npy" {
                if field.is_some() {
                    return Err("NPY arrays do not have named fields".into());
                }
                p.arrays.insert(name.into(), select_column(load_npy(&path)?, column)?);
            } else if ext == ".csv" || ext == ".txt" || ext == ".dat" {
                if field.is_some() {
                    return Err("Numeric text tables do not have named fields; select a column instead".into());
                }
                let mut table = load_table(&path, ext == ".csv")?;
                if column.is_none() && table.shape[1] == 1 {
                    table.shape = vec![table.shape[0]];
                }
                p.arrays.insert(name.into(), select_column(table, column)?);
            } else {
                if column.is_some() {
                    return Err("Column selection is supported only for JSON/NPY/numeric tables".into());
                }
                let selected = match (&field, name) {
                    (Some(f), _) => f.clone(),
                    (None, "velocities") => "velocities".into(),
                    (None, "masses") => "masses".into(),
                    (None, "cell" | "reference_cell" | "current_cell") => "cell".into(),
                    _ => "positions".into(),
                };
                p.arrays.insert(name.into(), frames_field(&read_frames(&path, reader)?, &selected, timestep)?);
                if name == "cell" && selected == "cell" {
                    p.cell_from_structure_file = true;
                }
            }
        }
        Ok(p)
    }

    /// An explicit null selects the documented default, as None does in Python.
    pub fn has(&self, name: &str) -> bool {
        if let Some(value) = self.json.get(name) {
            return !value.is_null();
        }
        self.arrays.contains_key(name) || self.structures.contains_key(name) || self.frames.contains_key(name)
    }

    pub fn json(&self, name: &str) -> Result<&Json, String> {
        self.json.get(name).ok_or_else(|| format!("Provide {name}"))
    }

    pub fn array(&self, name: &str) -> Result<NdArray, String> {
        if let Some(found) = self.arrays.get(name) {
            return Ok(found.clone());
        }
        to_array(self.json(name)?, name)
    }

    pub fn number(&self, name: &str) -> Result<f64, String> {
        if let Some(Json::Number(n)) = self.json.get(name) {
            if !n.is_finite() {
                return Err(format!("{name} must be finite"));
            }
            return Ok(*n);
        }
        let value = self.array(name)?;
        if value.size() != 1 || !value.values[0].is_finite() {
            return Err(format!("{name} must be one finite number"));
        }
        Ok(value.values[0])
    }

    pub fn number_or(&self, name: &str, fallback: f64) -> Result<f64, String> {
        if !self.has(name) {
            return Ok(fallback);
        }
        self.number(name)
    }

    pub fn boolean(&self, name: &str, fallback: bool) -> Result<bool, String> {
        if !self.has(name) {
            return Ok(fallback);
        }
        match self.json(name)? {
            Json::Bool(b) => Ok(*b),
            _ => Err(format!("{name} must be a boolean, not a string or number")),
        }
    }

    /// A cell read from a structure file comes with its periodicity.
    pub fn pbc(&self, name: &str, fallback: Pbc) -> Result<Pbc, String> {
        if !self.has(name) && name == "pbc" && self.cell_from_structure_file {
            return Ok([true; 3]);
        }
        if !self.has(name) {
            return Ok(fallback);
        }
        parse_pbc(self.json(name)?)
    }

    pub fn structure(&self, name: &str) -> Result<&StructureInput, String> {
        self.structures.get(name).ok_or_else(|| format!("Provide {name}"))
    }

    pub fn frames(&self, name: &str) -> Result<&Vec<FrameData>, String> {
        self.frames.get(name).ok_or_else(|| format!("Provide {name}"))
    }
}
