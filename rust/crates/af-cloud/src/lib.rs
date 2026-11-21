//! Out-of-core atom clouds (.afcloud) for structures far larger than memory
//! (the C++ cloud/AtomCloud): billions of atoms, as in large molecular-dynamics
//! snapshots.
//!
//! Atoms are stored by space in chunks of up to ~0.5 M atoms; each atom takes
//! 8 bytes (x, y, z quantized to 16 bits within its chunk's box and a 16-bit
//! species). Inside a chunk the atoms are in random order, so any prefix of a
//! chunk is a uniform random subsample of it: a renderer reads only as much of
//! each chunk as its size on screen needs.
//!
//! Layout: "AFCLOUD1", a little-endian u64 header offset and length, the chunk
//! data, then the header as JSON (species, chunks with offsets and boxes).
//! With one build thread the files are byte-identical to the C++ ones; with
//! more, chunk data is laid out in completion order (the header is the same).

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use af_core::cxx::nth_element_by;
use af_core::cxx_random::{shuffle, Mt19937_64};
use af_core::format::{format_g, stod};
use af_core::json::Json;

const MAGIC: &[u8; 8] = b"AFCLOUD1";
/// Records a bucket may hold in memory.
const IN_MEMORY_LIMIT: u64 = 8 << 20;

#[derive(Clone, Debug, PartialEq)]
pub struct Species {
    /// Element symbol or "type N".
    pub name: String,
    pub color: [f32; 3],
    /// Display radius (Angstrom).
    pub radius: f32,
}

impl Default for Species {
    fn default() -> Self {
        Self { name: String::new(), color: [0.7; 3], radius: 1.3 }
    }
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PackedAtom {
    pub x: u16,
    pub y: u16,
    pub z: u16,
    pub species: u16,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Chunk {
    /// File offset of its first atom.
    pub offset: u64,
    pub count: u64,
    /// Box of its atoms (quantization frame).
    pub lower: [f64; 3],
    pub upper: [f64; 3],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CloudInfo {
    pub atoms: u64,
    pub lower: [f64; 3],
    pub upper: [f64; 3],
    pub has_cell: bool,
    pub cell: [[f64; 3]; 3],
    pub cell_origin: [f64; 3],
    pub species: Vec<Species>,
    pub chunks: Vec<Chunk>,
    /// What it was built from.
    pub source: String,
}

/// A bucket record on disk (16 bytes).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
struct Record {
    x: f32,
    y: f32,
    z: f32,
    species: u16,
    pad: u16,
}

fn as_bytes<T: Copy>(items: &[T]) -> &[u8] {
    // Plain repr(C) records of integers and floats.
    unsafe { std::slice::from_raw_parts(items.as_ptr() as *const u8, std::mem::size_of_val(items)) }
}

fn as_bytes_mut<T: Copy>(items: &mut [T]) -> &mut [u8] {
    unsafe { std::slice::from_raw_parts_mut(items.as_mut_ptr() as *mut u8, std::mem::size_of_val(items)) }
}

fn display(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn element_number(symbol: &str) -> i32 {
    if let Some(z) = (1..=118).find(|&z| af_core::elements::symbol(z) == symbol) {
        return z;
    }
    // Lower/upper-case variants ("CU", "cu").
    if symbol.is_empty() || symbol.len() > 3 {
        return 0;
    }
    let mut fixed: Vec<u8> = symbol.bytes().collect();
    fixed[0] = fixed[0].to_ascii_uppercase();
    for c in fixed.iter_mut().skip(1) {
        *c = c.to_ascii_lowercase();
    }
    let fixed = String::from_utf8_lossy(&fixed).into_owned();
    (1..=118).find(|&z| af_core::elements::symbol(z) == fixed).unwrap_or(0)
}

fn species_for(name: &str, index: usize) -> Species {
    let z = element_number(name);
    if z != 0 {
        let radius = af_core::elements::literature_covalent_radii()[z as usize];
        return Species { name: name.to_string(), color: af_core::elements::default_colour(z), radius: if 0.3 < radius { radius } else { 0.3 } };
    }
    // Distinct colours for LAMMPS types without element names.
    const PALETTE: [[f32; 3]; 8] = [
        [0.85, 0.53, 0.20],
        [0.25, 0.55, 0.90],
        [0.35, 0.75, 0.35],
        [0.85, 0.30, 0.35],
        [0.65, 0.45, 0.85],
        [0.90, 0.80, 0.30],
        [0.40, 0.80, 0.80],
        [0.60, 0.60, 0.60],
    ];
    Species { name: name.to_string(), color: PALETTE[index % 8], radius: 1.3 }
}

fn vec_json(v: &[f64; 3]) -> Json {
    Json::Array(vec![v[0].into(), v[1].into(), v[2].into()])
}

fn vec_of(j: &Json) -> Result<[f64; 3], String> {
    let items = j.items()?;
    if items.len() < 3 {
        return Err("vector::_M_range_check".into());
    }
    Ok([items[0].number()?, items[1].number()?, items[2].number()?])
}

fn read_u64(file: &mut File) -> u64 {
    let mut bytes = [0u8; 8];
    let _ = file.read_exact(&mut bytes);
    u64::from_le_bytes(bytes)
}

/// std::from_chars for a double over a whole field (no leading '+').
fn parse_double(field: &[u8]) -> Option<f64> {
    if field.is_empty() || field[0] == b'+' {
        return None;
    }
    std::str::from_utf8(field).ok()?.parse::<f64>().ok()
}

/// std::stoull (leading space, optional sign, digits).
fn stoull(text: &str) -> Result<u64, String> {
    let t = text.trim_start_matches([' ', '\t', '\n', '\x0B', '\x0C', '\r']);
    let (negative, digits) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let end = digits.bytes().take_while(u8::is_ascii_digit).count();
    if end == 0 {
        return Err("stoull".into());
    }
    let value: u64 = digits[..end].parse().map_err(|_| "stoull".to_string())?;
    Ok(if negative { value.wrapping_neg() } else { value })
}

fn stod_prefix(text: &str) -> Result<f64, String> {
    stod(text).map(|(v, _)| v).map_err(|_| "stod".to_string())
}

fn stoi_prefix(text: &str) -> Result<i32, String> {
    let t = text.trim_start_matches([' ', '\t', '\n', '\x0B', '\x0C', '\r']);
    let sign = usize::from(t.starts_with(['+', '-']));
    let end = sign + t[sign..].bytes().take_while(u8::is_ascii_digit).count();
    if end == sign {
        return Err("stoi".into());
    }
    t[..end].parse::<i64>().ok().and_then(|v| i32::try_from(v).ok()).ok_or_else(|| "stoi".to_string())
}

fn words(line: &str) -> Vec<&str> {
    line.split([' ', '\t', '\n', '\x0B', '\x0C', '\r']).filter(|w| !w.is_empty()).collect()
}

fn thread_count(requested: u32) -> u32 {
    if requested != 0 {
        return requested;
    }
    let hardware = std::thread::available_parallelism().map_or(1, |n| n.get() as u32).max(1);
    hardware.min(16).max(1)
}

// ---------------------------------------------------------------- reading

/// A read-only .afcloud file; reads are thread-safe.
pub struct CloudFile {
    path: PathBuf,
    info: CloudInfo,
    file: Mutex<File>,
}

impl CloudFile {
    pub fn open(path: &Path) -> Result<Self, String> {
        let mut file = File::open(path).map_err(|_| format!("Cannot open {}", display(path)))?;
        let mut magic = [0u8; 8];
        if file.read_exact(&mut magic).is_err() || &magic != MAGIC {
            return Err(format!("{} is not an AtomForge atom cloud", display(path)));
        }
        let header_offset = read_u64(&mut file);
        let header_length = read_u64(&mut file);
        if header_offset == 0 || header_length > (1u64 << 32) {
            return Err(format!("{} is incomplete (its build did not finish)", display(path)));
        }
        let mut text = vec![0u8; header_length as usize];
        file.seek(SeekFrom::Start(header_offset)).map_err(|_| format!("Cannot read the header of {}", display(path)))?;
        file.read_exact(&mut text).map_err(|_| format!("Cannot read the header of {}", display(path)))?;
        let header = Json::parse(&String::from_utf8_lossy(&text))?;
        let mut info = CloudInfo {
            atoms: header.at("atoms")?.number()? as u64,
            lower: vec_of(header.at("lower")?)?,
            upper: vec_of(header.at("upper")?)?,
            ..Default::default()
        };
        if let Some(source) = header.find("source") {
            info.source = source.string()?.to_string();
        }
        if let Some(cell) = header.find("cell") {
            info.has_cell = true;
            let rows = cell.items()?;
            for r in 0..3 {
                info.cell[r] = vec_of(rows.get(r).ok_or("vector::_M_range_check")?)?;
            }
            info.cell_origin = vec_of(header.at("cell_origin")?)?;
        }
        for s in header.at("species")?.items()? {
            let c = vec_of(s.at("color")?)?;
            info.species.push(Species {
                name: s.at("name")?.string()?.to_string(),
                color: [c[0] as f32, c[1] as f32, c[2] as f32],
                radius: s.at("radius")?.number()? as f32,
            });
        }
        for c in header.at("chunks")?.items()? {
            info.chunks.push(Chunk {
                offset: c.at("offset")?.number()? as u64,
                count: c.at("count")?.number()? as u64,
                lower: vec_of(c.at("lower")?)?,
                upper: vec_of(c.at("upper")?)?,
            });
        }
        Ok(Self { path: path.to_path_buf(), info, file: Mutex::new(file) })
    }

    pub fn info(&self) -> &CloudInfo {
        &self.info
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Atoms [first, first + count) of a chunk.
    pub fn read(&self, chunk: usize, first: u64, count: u64) -> Result<Vec<PackedAtom>, String> {
        let c = self.info.chunks.get(chunk).ok_or("vector::_M_range_check")?;
        if first + count > c.count {
            return Err("Atom cloud read past the end of a chunk".into());
        }
        let mut out = vec![PackedAtom::default(); count as usize];
        let mut file = self.file.lock().unwrap();
        file.seek(SeekFrom::Start(c.offset + first * 8))
            .and_then(|_| file.read_exact(as_bytes_mut(&mut out)))
            .map_err(|_| format!("Cannot read atoms from {}", display(&self.path)))?;
        Ok(out)
    }

    /// A packed atom of a chunk decoded to a position.
    pub fn position(&self, chunk: usize, atom: &PackedAtom) -> Result<[f64; 3], String> {
        let c = self.info.chunks.get(chunk).ok_or("vector::_M_range_check")?;
        let q = [atom.x, atom.y, atom.z];
        Ok(std::array::from_fn(|k| c.lower[k] + (c.upper[k] - c.lower[k]) * q[k] as f64 / 65535.0))
    }
}

/// Whether a file is (by its first bytes) an .afcloud.
pub fn is_cloud_file(path: &Path) -> bool {
    let mut magic = [0u8; 8];
    File::open(path).and_then(|mut f| f.read_exact(&mut magic)).is_ok() && &magic == MAGIC
}

/// Atom count of a LAMMPS dump or XYZ file from its first header (0 if unknown).
pub fn peek_atom_count(path: &Path) -> u64 {
    let Ok(file) = File::open(path) else { return 0 };
    let mut reader = BufReader::new(file);
    let mut line = Vec::new();
    let mut getline = |line: &mut Vec<u8>| -> bool {
        line.clear();
        match reader.read_until(b'\n', line) {
            Ok(0) | Err(_) => false,
            Ok(_) => {
                if line.last() == Some(&b'\n') {
                    line.pop();
                }
                true
            }
        }
    };
    for i in 0..64 {
        if !getline(&mut line) {
            break;
        }
        let text = String::from_utf8_lossy(&line).into_owned();
        if i == 0 {
            // XYZ: the first line is the atom count.
            let w = words(&text);
            if w.len() == 1 {
                if let Ok(n) = stoull(w[0]) {
                    if w[0].bytes().skip(usize::from(w[0].starts_with(['+', '-']))).all(|c| c.is_ascii_digit()) {
                        return n;
                    }
                }
            }
        }
        if text.starts_with("ITEM: NUMBER OF ATOMS") {
            if !getline(&mut line) {
                return 0;
            }
            return stoull(&String::from_utf8_lossy(&line)).unwrap_or(0);
        }
    }
    0
}

// ---------------------------------------------------------------- building

/// Progress of long operations: fraction (0..1, or < 0 when unknown) and a
/// message; return false to cancel.
pub type Progress<'a> = Option<&'a (dyn Fn(f64, &str) -> bool + Sync)>;

#[derive(Clone, Debug)]
pub struct BuildOptions {
    /// Target atoms per chunk.
    pub chunk_atoms: u64,
    /// Folder for the bucket files (default: next to the output).
    pub temporary: Option<PathBuf>,
    /// 0 = hardware concurrency (at most 16).
    pub threads: u32,
    /// Shuffle seed (reproducible files).
    pub seed: u64,
}

impl Default for BuildOptions {
    fn default() -> Self {
        Self { chunk_atoms: 1 << 19, temporary: None, threads: 0, seed: 12345 }
    }
}

/// One atom handed to the builder.
#[derive(Clone, Copy, Debug, Default)]
pub struct SourceAtom {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub species: u16,
}

struct Bucket {
    path: PathBuf,
    pending: Vec<Record>,
    count: u64,
}

/// Streams atoms into a .afcloud without holding them in memory: atoms are
/// spread over spatial buckets on disk, then each bucket is split into chunks,
/// shuffled, quantized and written.
pub struct CloudBuilder {
    output: PathBuf,
    folder: PathBuf,
    options: BuildOptions,
    lower: [f64; 3],
    upper: [f64; 3],
    grid: [i32; 3],
    buckets: Vec<Bucket>,
    species: Vec<Species>,
    has_cell: bool,
    cell: [[f64; 3]; 3],
    cell_origin: [f64; 3],
    source: String,
    added: u64,
    finished: bool,
}

impl CloudBuilder {
    /// `lower`/`upper` bound the atoms (atoms outside are clamped into the
    /// grid); `expected_atoms` sizes the bucket grid (an estimate is fine).
    pub fn new(output: &Path, lower: [f64; 3], upper: [f64; 3], expected_atoms: u64, mut options: BuildOptions) -> Result<Self, String> {
        if options.chunk_atoms < 1024 {
            options.chunk_atoms = 1024;
        }
        let folder = options.temporary.clone().unwrap_or_else(|| output.with_extension("afcloud-build"));
        std::fs::create_dir_all(&folder).map_err(|e| e.to_string())?;
        // Buckets of a few chunks each, laid out in proportion to the box.
        let buckets = (expected_atoms as f64 / (4.0 * options.chunk_atoms as f64)).ceil().clamp(1.0, 4096.0);
        let mut extent = [0.0f64; 3];
        let mut volume = 1.0;
        let mut used = 0;
        for k in 0..3 {
            let e = upper[k] - lower[k];
            extent[k] = if 0.0 < e { e } else { 0.0 };
            if extent[k] > 1e-9 {
                volume *= extent[k];
                used += 1;
            }
        }
        let per_axis = if used > 0 { (buckets / volume).powf(1.0 / used as f64) } else { 0.0 };
        let mut grid = [1i32; 3];
        let mut total = 1usize;
        for k in 0..3 {
            grid[k] = if extent[k] > 1e-9 { ((extent[k] * per_axis).round() as i64).clamp(1, 64) as i32 } else { 1 };
            total *= grid[k] as usize;
        }
        let buckets = (0..total)
            .map(|i| {
                let path = folder.join(format!("bucket-{i}.bin"));
                let _ = std::fs::remove_file(&path);
                Bucket { path, pending: Vec::new(), count: 0 }
            })
            .collect();
        Ok(Self {
            output: output.to_path_buf(),
            folder,
            options,
            lower,
            upper,
            grid,
            buckets,
            species: Vec::new(),
            has_cell: false,
            cell: [[0.0; 3]; 3],
            cell_origin: [0.0; 3],
            source: String::new(),
            added: 0,
            finished: false,
        })
    }

    /// The species index of a name (added on first use; colours and radii from the element table).
    pub fn species(&mut self, name: &str) -> Result<u16, String> {
        if let Some(i) = self.species.iter().position(|s| s.name == name) {
            return Ok(i as u16);
        }
        if self.species.len() >= 65535 {
            return Err("Too many species in an atom cloud".into());
        }
        self.species.push(species_for(name, self.species.len()));
        Ok((self.species.len() - 1) as u16)
    }

    pub fn set_species(&mut self, species: Vec<Species>) {
        self.species = species;
    }

    pub fn set_cell(&mut self, cell: [[f64; 3]; 3], origin: [f64; 3]) {
        self.has_cell = true;
        self.cell = cell;
        self.cell_origin = origin;
    }

    pub fn set_source(&mut self, source: &str) {
        self.source = source.to_string();
    }

    pub fn added(&self) -> u64 {
        self.added
    }

    fn flush(bucket: &mut Bucket) -> Result<(), String> {
        if bucket.pending.is_empty() {
            return Ok(());
        }
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&bucket.path)
            .map_err(|_| format!("Cannot write {}", display(&bucket.path)))?;
        file.write_all(as_bytes(&bucket.pending)).map_err(|_| format!("Cannot write {} (disk full?)", display(&bucket.path)))?;
        bucket.pending.clear();
        Ok(())
    }

    pub fn add(&mut self, atom: &SourceAtom) -> Result<(), String> {
        let (mut index, mut stride) = (0usize, 1usize);
        let p = [atom.x, atom.y, atom.z];
        for k in 0..3 {
            let e = self.upper[k] - self.lower[k];
            let mut cell = if e > 1e-9 { ((p[k] - self.lower[k]) / e * self.grid[k] as f64) as i32 } else { 0 };
            cell = cell.clamp(0, self.grid[k] - 1);
            index += cell as usize * stride;
            stride *= self.grid[k] as usize;
        }
        let bucket = &mut self.buckets[index];
        bucket.pending.push(Record { x: atom.x as f32, y: atom.y as f32, z: atom.z as f32, species: atom.species, pad: 0 });
        bucket.count += 1;
        self.added += 1;
        if bucket.pending.len() >= 16384 {
            Self::flush(bucket)?;
        }
        Ok(())
    }

    pub fn add_all(&mut self, atoms: &[SourceAtom]) -> Result<(), String> {
        for atom in atoms {
            self.add(atom)?;
        }
        Ok(())
    }

    /// Writes the file; returns its description.
    pub fn finish(&mut self, progress: Progress) -> Result<CloudInfo, String> {
        if self.finished {
            return Err("CloudBuilder::finish called twice".into());
        }
        self.finished = true;
        for bucket in self.buckets.iter_mut() {
            Self::flush(bucket)?;
        }
        let mut out = File::create(&self.output).map_err(|_| format!("Cannot write {}", display(&self.output)))?;
        let mut header_bytes = MAGIC.to_vec();
        header_bytes.extend_from_slice(&[0u8; 16]); // header offset and length, written at the end
        out.write_all(&header_bytes).map_err(|_| format!("Cannot write {}", display(&self.output)))?;

        // Buckets (and parts of oversized ones) to turn into chunks.
        struct Job {
            path: PathBuf,
            count: u64,
            order: u64,
        }
        struct Shared {
            jobs: Vec<Job>,
            next: usize,
            failure: Option<String>,
            cancelled: bool,
        }
        struct Written {
            out: File,
            offset: u64,
            chunks: Vec<Chunk>,
            chunk_order: Vec<(u64, usize)>,
            done: u64,
        }
        let jobs = self
            .buckets
            .iter()
            .enumerate()
            .filter(|(_, b)| b.count > 0)
            .map(|(i, b)| Job { path: b.path.clone(), count: b.count, order: i as u64 })
            .collect();
        let shared = Mutex::new(Shared { jobs, next: 0, failure: None, cancelled: false });
        let written = Mutex::new(Written { out, offset: 24, chunks: Vec::new(), chunk_order: Vec::new(), done: 0 });
        let (seed, chunk_atoms, added, output) = (self.options.seed, self.options.chunk_atoms, self.added, self.output.clone());
        let worker = || {
            let run = || -> Result<(), String> {
                loop {
                    let job = {
                        let mut s = shared.lock().unwrap();
                        if s.next >= s.jobs.len() || s.cancelled || s.failure.is_some() {
                            return Ok(());
                        }
                        let j = &s.jobs[s.next];
                        let job = Job { path: j.path.clone(), count: j.count, order: j.order };
                        s.next += 1;
                        job
                    };
                    if job.count > IN_MEMORY_LIMIT {
                        // Too large for memory: split on disk and queue the parts.
                        let parts = split_file(&job.path, job.count)?;
                        let mut s = shared.lock().unwrap();
                        for (p, (path, count)) in parts.into_iter().enumerate() {
                            if count > 0 {
                                s.jobs.push(Job { path, count, order: job.order * 8 + p as u64 + (1 << 20) });
                            }
                        }
                        continue;
                    }
                    let mut records = read_records(&job.path, job.count)?;
                    let _ = std::fs::remove_file(&job.path);
                    let mut random = Mt19937_64::new(seed ^ 0x9E3779B97F4A7C15u64.wrapping_mul(job.order + 1));
                    let mut built = Vec::new();
                    split(&mut records, chunk_atoms, &mut random, &mut built);
                    drop(records);
                    let mut w = written.lock().unwrap();
                    for (c, (atoms, mut chunk)) in built.into_iter().enumerate() {
                        chunk.offset = w.offset;
                        w.out.write_all(as_bytes(&atoms)).map_err(|_| format!("Cannot write {} (disk full?)", display(&output)))?;
                        w.offset += atoms.len() as u64 * 8;
                        w.done += chunk.count;
                        w.chunks.push(chunk);
                        w.chunk_order.push((job.order, c));
                    }
                    let fraction = if added > 0 { w.done as f64 / added as f64 } else { 1.0 };
                    drop(w);
                    if let Some(progress) = progress {
                        if !progress(fraction, "Writing chunks") {
                            shared.lock().unwrap().cancelled = true;
                        }
                    }
                }
            };
            if let Err(e) = run() {
                let mut s = shared.lock().unwrap();
                if s.failure.is_none() {
                    s.failure = Some(e);
                }
            }
        };
        // Workers pull jobs until none are left (a worker that splits a bucket queues and takes its parts).
        loop {
            let threads = thread_count(self.options.threads);
            if threads == 1 {
                worker();
            } else {
                std::thread::scope(|scope| {
                    for _ in 0..threads {
                        scope.spawn(worker);
                    }
                });
            }
            let s = shared.lock().unwrap();
            if let Some(e) = &s.failure {
                return Err(e.clone());
            }
            if s.cancelled {
                return Err("Cancelled".into());
            }
            if s.next >= s.jobs.len() {
                break;
            }
        }
        let Written { mut out, offset, chunks, chunk_order, done } = written.into_inner().unwrap();

        // A stable chunk order (by bucket) regardless of thread timing.
        let mut order: Vec<usize> = (0..chunks.len()).collect();
        order.sort_by(|&a, &b| chunk_order[a].cmp(&chunk_order[b]));
        let mut info = CloudInfo {
            atoms: done,
            species: self.species.clone(),
            has_cell: self.has_cell,
            cell: self.cell,
            cell_origin: self.cell_origin,
            source: self.source.clone(),
            lower: [f64::MAX; 3],
            upper: [f64::MIN; 3],
            chunks: Vec::new(),
        };
        let mut chunk_list = Vec::new();
        for &i in &order {
            let c = chunks[i];
            info.chunks.push(c);
            for k in 0..3 {
                if c.lower[k] < info.lower[k] {
                    info.lower[k] = c.lower[k];
                }
                if info.upper[k] < c.upper[k] {
                    info.upper[k] = c.upper[k];
                }
            }
            let mut entry = Json::object();
            entry.set("offset", c.offset as f64);
            entry.set("count", c.count as f64);
            entry.set("lower", vec_json(&c.lower));
            entry.set("upper", vec_json(&c.upper));
            chunk_list.push(entry);
        }
        if info.chunks.is_empty() {
            info.lower = [0.0; 3];
            info.upper = [0.0; 3];
        }
        let mut header = Json::object();
        header.set("format", "atomforge-cloud");
        header.set("version", 1);
        header.set("atoms", info.atoms as f64);
        header.set("lower", vec_json(&info.lower));
        header.set("upper", vec_json(&info.upper));
        header.set("source", info.source.as_str());
        if self.has_cell {
            header.set("cell", Json::Array(self.cell.iter().map(vec_json).collect()));
            header.set("cell_origin", vec_json(&self.cell_origin));
        }
        let species_list: Vec<Json> = self
            .species
            .iter()
            .map(|s| {
                let mut entry = Json::object();
                entry.set("name", s.name.as_str());
                entry.set("color", Json::Array(s.color.iter().map(|&c| Json::from(c)).collect()));
                entry.set("radius", s.radius);
                entry
            })
            .collect();
        header.set("species", Json::Array(species_list));
        header.set("chunks", Json::Array(chunk_list));
        let text = header.dump(-1);
        let fail = |_| format!("Cannot write {}", display(&self.output));
        out.write_all(text.as_bytes()).map_err(fail)?;
        out.seek(SeekFrom::Start(8)).map_err(fail)?;
        out.write_all(&offset.to_le_bytes()).map_err(fail)?;
        out.write_all(&(text.len() as u64).to_le_bytes()).map_err(fail)?;
        drop(out);
        if let Some(progress) = progress {
            progress(1.0, "Done");
        }
        Ok(info)
    }
}

impl Drop for CloudBuilder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.folder);
    }
}

fn bounds_of(records: &[Record]) -> ([f64; 3], [f64; 3]) {
    let (mut lower, mut upper) = ([f64::MAX; 3], [f64::MIN; 3]);
    for r in records {
        let p = [r.x as f64, r.y as f64, r.z as f64];
        for k in 0..3 {
            if p[k] < lower[k] {
                lower[k] = p[k];
            }
            if upper[k] < p[k] {
                upper[k] = p[k];
            }
        }
    }
    (lower, upper)
}

fn coordinate(r: &Record, axis: usize) -> f32 {
    match axis {
        0 => r.x,
        1 => r.y,
        _ => r.z,
    }
}

/// Splits records into chunks of at most `limit` atoms (median cuts along the longest axis).
fn split(r: &mut [Record], limit: u64, random: &mut Mt19937_64, out: &mut Vec<(Vec<PackedAtom>, Chunk)>) {
    let n = r.len();
    let (lower, upper) = bounds_of(r);
    if n as u64 > limit {
        let mut axis = 0;
        for k in 1..3 {
            if upper[k] - lower[k] > upper[axis] - lower[axis] {
                axis = k;
            }
        }
        let half = n / 2;
        nth_element_by(r, half, |a, b| coordinate(a, axis) < coordinate(b, axis));
        let (left, right) = r.split_at_mut(half);
        split(left, limit, random, out);
        split(right, limit, random, out);
        return;
    }
    // Random order: any prefix is a uniform subsample.
    shuffle(r, random);
    let mut scale = [0.0f64; 3];
    for k in 0..3 {
        let e = upper[k] - lower[k];
        scale[k] = if e > 0.0 { 65535.0 / e } else { 0.0 };
    }
    let atoms = r
        .iter()
        .map(|record| {
            let p = [record.x as f64, record.y as f64, record.z as f64];
            let q: [u16; 3] = std::array::from_fn(|k| (((p[k] - lower[k]) * scale[k]).round() as i64).clamp(0, 65535) as u16);
            PackedAtom { x: q[0], y: q[1], z: q[2], species: record.species }
        })
        .collect();
    out.push((atoms, Chunk { offset: 0, count: n as u64, lower, upper }));
}

fn read_records(path: &Path, count: u64) -> Result<Vec<Record>, String> {
    let mut records = vec![Record::default(); count as usize];
    if count == 0 {
        return Ok(records);
    }
    let mut file = File::open(path).map_err(|_| format!("Cannot read {}", display(path)))?;
    file.read_exact(as_bytes_mut(&mut records)).map_err(|_| format!("Short read from {}", display(path)))?;
    Ok(records)
}

/// Splits a bucket file too large for memory into eight files by octant, streaming.
fn split_file(path: &Path, count: u64) -> Result<Vec<(PathBuf, u64)>, String> {
    let (mut lower, mut upper) = ([f64::MAX; 3], [f64::MIN; 3]);
    let mut block = vec![Record::default(); 1 << 20];
    let stream = |block: &mut Vec<Record>, mut using: Box<dyn FnMut(&[Record]) -> Result<(), String> + '_>| -> Result<(), String> {
        let mut file = File::open(path).map_err(|_| format!("Cannot read {}", display(path)))?;
        loop {
            let bytes = as_bytes_mut(block.as_mut_slice());
            let mut got = 0;
            while got < bytes.len() {
                match file.read(&mut bytes[got..]) {
                    Ok(0) => break,
                    Ok(n) => got += n,
                    Err(_) => break,
                }
            }
            let records = got / 16;
            if records == 0 {
                return Ok(());
            }
            using(&block[..records])?;
        }
    };
    stream(
        &mut block,
        Box::new(|r: &[Record]| {
            let (lo, hi) = bounds_of(r);
            for k in 0..3 {
                if lo[k] < lower[k] {
                    lower[k] = lo[k];
                }
                if upper[k] < hi[k] {
                    upper[k] = hi[k];
                }
            }
            Ok(())
        }),
    )?;
    let mid = [0.5 * (lower[0] + upper[0]), 0.5 * (lower[1] + upper[1]), 0.5 * (lower[2] + upper[2])];
    let mut parts: Vec<(PathBuf, u64)> = (0..8)
        .map(|i| {
            let p = PathBuf::from(format!("{}.{i}", path.to_string_lossy()));
            let _ = std::fs::remove_file(&p);
            (p, 0)
        })
        .collect();
    let mut pending: Vec<Vec<Record>> = vec![Vec::new(); 8];
    let flush_part = |parts: &Vec<(PathBuf, u64)>, pending: &mut Vec<Vec<Record>>, i: usize| -> Result<(), String> {
        if pending[i].is_empty() {
            return Ok(());
        }
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&parts[i].0)
            .map_err(|_| "Cannot write a bucket part".to_string())?;
        let _ = file.write_all(as_bytes(&pending[i]));
        pending[i].clear();
        Ok(())
    };
    stream(
        &mut block,
        Box::new(|r: &[Record]| {
            for record in r {
                let octant = usize::from(record.x as f64 > mid[0]) | (usize::from(record.y as f64 > mid[1]) << 1) | (usize::from(record.z as f64 > mid[2]) << 2);
                pending[octant].push(*record);
                parts[octant].1 += 1;
                if pending[octant].len() >= 65536 {
                    flush_part(&parts, &mut pending, octant)?;
                }
            }
            Ok(())
        }),
    )?;
    for i in 0..8 {
        flush_part(&parts, &mut pending, i)?;
    }
    let _ = std::fs::remove_file(path);
    // All atoms at one point cannot be split by space.
    if parts.iter().any(|p| p.1 == count) {
        return Err("Too many atoms at one position to split into chunks".into());
    }
    Ok(parts)
}

// ---------------------------------------------------------------- converters

/// Reads a text file line by line, tracking the bytes consumed.
struct LineReader {
    reader: BufReader<File>,
    size: u64,
    consumed: u64,
    line: Vec<u8>,
}

impl LineReader {
    fn new(path: &Path) -> Result<Self, String> {
        let file = File::open(path).map_err(|_| format!("Cannot read {}", display(path)))?;
        let size = file.metadata().map(|m| m.len()).unwrap_or(0);
        Ok(Self { reader: BufReader::with_capacity(1 << 20, file), size, consumed: 0, line: Vec::new() })
    }

    /// The next line without its newline (and a final '\r'); None at the end.
    fn next(&mut self) -> Option<&[u8]> {
        self.line.clear();
        let n = self.reader.read_until(b'\n', &mut self.line).ok()?;
        if n == 0 {
            return None;
        }
        self.consumed += n as u64;
        if self.line.last() == Some(&b'\n') {
            self.line.pop();
        }
        if self.line.last() == Some(&b'\r') {
            self.line.pop();
        }
        Some(&self.line)
    }

    fn line(&mut self) -> Result<String, String> {
        self.next().map(|l| String::from_utf8_lossy(l).into_owned()).ok_or_else(|| "Unexpected end of file".to_string())
    }

    fn fraction(&self) -> f64 {
        if self.size > 0 {
            self.consumed as f64 / self.size as f64
        } else {
            -1.0
        }
    }
}

fn is_field_space(c: u8) -> bool {
    c == b' ' || c == b'\t' || c == b'\r'
}

/// Up to `count` whitespace-separated fields of a line (empty ones past its end).
fn fields(line: &[u8], count: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::with_capacity(count);
    let mut p = 0;
    for _ in 0..count {
        while p < line.len() && is_field_space(line[p]) {
            p += 1;
        }
        let start = p;
        while p < line.len() && !is_field_space(line[p]) && line[p] != b'\n' {
            p += 1;
        }
        out.push((start, p));
    }
    out
}

fn min3(values: [f64; 4]) -> f64 {
    let mut m = values[0];
    for &v in &values[1..] {
        if v < m {
            m = v;
        }
    }
    m
}

fn max3(values: [f64; 4]) -> f64 {
    let mut m = values[0];
    for &v in &values[1..] {
        if m < v {
            m = v;
        }
    }
    m
}

fn file_name(path: &Path) -> String {
    path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default()
}

fn report(progress: Progress, fraction: f64, message: &str) -> Result<(), String> {
    match progress {
        Some(p) if !p(fraction, message) => Err("Cancelled".into()),
        _ => Ok(()),
    }
}

/// A LAMMPS text dump (first frame unless `frame` is given): columns x y z,
/// xs ys zs, xu yu zu or xsu ysu zsu, with element or type; `type_names` names
/// LAMMPS types ("Cu", "Ni" for types 1 and 2).
pub fn build_from_lammps_dump(input: &Path, output: &Path, type_names: &[String], frame: usize, progress: Progress, options: BuildOptions) -> Result<CloudInfo, String> {
    let mut reader = LineReader::new(input)?;
    let mut current = 0usize;
    loop {
        // Frame header.
        loop {
            if reader.line()?.starts_with("ITEM: TIMESTEP") {
                break;
            }
        }
        reader.line()?;
        if !reader.line()?.starts_with("ITEM: NUMBER OF ATOMS") {
            return Err(format!("Not a LAMMPS dump: {}", display(input)));
        }
        let atoms = stoull(&reader.line()?)?;
        let box_line = reader.line()?;
        if !box_line.starts_with("ITEM: BOX BOUNDS") {
            return Err("LAMMPS dump without BOX BOUNDS".into());
        }
        let triclinic = box_line.contains("xy");
        let mut bounds = [[0.0f64; 3]; 3];
        for row in bounds.iter_mut() {
            let line = reader.line()?;
            let w = words(&line);
            if w.len() < 2 {
                return Err("Bad BOX BOUNDS line".into());
            }
            row[0] = stod_prefix(w[0])?;
            row[1] = stod_prefix(w[1])?;
            row[2] = if w.len() > 2 { stod_prefix(w[2])? } else { 0.0 };
        }
        let columns_line = reader.line()?;
        let columns = words(&columns_line);
        if columns.len() < 3 || columns[0] != "ITEM:" || columns[1] != "ATOMS" {
            return Err("LAMMPS dump without an ATOMS section".into());
        }
        if current < frame {
            for _ in 0..atoms {
                reader.line()?;
            }
            current += 1;
            continue;
        }
        // Box: orthogonal or triclinic (LAMMPS bounding-box convention).
        let (xy, xz, yz) = if triclinic { (bounds[0][2], bounds[1][2], bounds[2][2]) } else { (0.0, 0.0, 0.0) };
        let xlo = bounds[0][0] - min3([0.0, xy, xz, xy + xz]);
        let xhi = bounds[0][1] - max3([0.0, xy, xz, xy + xz]);
        let ylo = bounds[1][0] - if yz < 0.0 { yz } else { 0.0 };
        let yhi = bounds[1][1] - if 0.0 < yz { yz } else { 0.0 };
        let (zlo, zhi) = (bounds[2][0], bounds[2][1]);
        // Columns of the positions, element and type.
        let (mut cx, mut cy, mut cz, mut element, mut type_) = (-1i32, -1i32, -1i32, -1i32, -1i32);
        let mut scaled = false;
        for (i, c) in columns.iter().enumerate().skip(2) {
            let index = i as i32 - 2;
            match *c {
                "x" | "xu" => cx = index,
                "y" | "yu" => cy = index,
                "z" | "zu" => cz = index,
                "element" => element = index,
                "type" => type_ = index,
                _ => {}
            }
        }
        if cx < 0 {
            for (i, c) in columns.iter().enumerate().skip(2) {
                let index = i as i32 - 2;
                match *c {
                    "xs" | "xsu" => cx = index,
                    "ys" | "ysu" => cy = index,
                    "zs" | "zsu" => cz = index,
                    _ => {}
                }
            }
            scaled = true;
        }
        if cx < 0 || cy < 0 || cz < 0 {
            return Err("The dump has no x y z (or xs ys zs) columns".into());
        }
        let mut builder = CloudBuilder::new(output, [bounds[0][0], bounds[1][0], bounds[2][0]], [bounds[0][1], bounds[1][1], bounds[2][1]], atoms, options)?;
        builder.set_cell([[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]], [xlo, ylo, zlo]);
        builder.set_source(&format!("LAMMPS dump {}, frame {frame}", file_name(input)));
        // Type and element names map to species indices through small caches.
        let mut type_species: Vec<i32> = Vec::new();
        let mut element_species: Vec<(Vec<u8>, u16)> = Vec::new();
        let needed = [cx, cy, cz, element, type_].into_iter().max().unwrap() as usize;
        let mut batch: Vec<SourceAtom> = Vec::with_capacity(1 << 16);
        for i in 0..atoms {
            let input_name = file_name(input);
            let Some(line) = reader.next() else {
                return Err(format!("The dump ends after {i} of {atoms} atoms"));
            };
            let line = line.to_vec();
            let f = fields(&line, needed + 1);
            if f.iter().any(|(a, b)| a == b) {
                return Err(format!("Short atom line in the dump: {}", String::from_utf8_lossy(&line)));
            }
            let mut r = [0.0f64; 3];
            for (k, col) in [cx, cy, cz].into_iter().enumerate() {
                let (a, b) = f[col as usize];
                r[k] = parse_double(&line[a..b]).ok_or_else(|| format!("Bad coordinate in the dump: {}", String::from_utf8_lossy(&line[a..b])))?;
            }
            if scaled {
                let (s0, s1, s2) = (r[0], r[1], r[2]);
                r[0] = xlo + s0 * (xhi - xlo) + s1 * xy + s2 * xz;
                r[1] = ylo + s1 * (yhi - ylo) + s2 * yz;
                r[2] = zlo + s2 * (zhi - zlo);
            }
            let s = if element >= 0 {
                let (a, b) = f[element as usize];
                let name = &line[a..b];
                match element_species.iter().find(|(n, _)| n.as_slice() == name) {
                    Some(&(_, index)) => index,
                    None => {
                        let index = builder.species(&String::from_utf8_lossy(name))?;
                        element_species.push((name.to_vec(), index));
                        index
                    }
                }
            } else if type_ >= 0 {
                let (a, b) = f[type_ as usize];
                // std::from_chars: the leading integer, else 0.
                let field = &line[a..b];
                let sign = usize::from(field.first() == Some(&b'-'));
                let digits = field[sign..].iter().take_while(|c| c.is_ascii_digit()).count();
                let mut t: i64 = if digits > 0 {
                    std::str::from_utf8(&field[..sign + digits]).ok().and_then(|s| s.parse::<i32>().ok()).map_or(0, i64::from)
                } else {
                    0
                };
                if !(0..=65535).contains(&t) {
                    t = 0;
                }
                let t = t as usize;
                if t >= type_species.len() {
                    type_species.resize(t + 1, -1);
                }
                if type_species[t] < 0 {
                    let name = if t >= 1 && t <= type_names.len() { type_names[t - 1].clone() } else { format!("type {t}") };
                    type_species[t] = i32::from(builder.species(&name)?);
                }
                type_species[t] as u16
            } else {
                builder.species("X")?
            };
            batch.push(SourceAtom { x: r[0], y: r[1], z: r[2], species: s });
            if batch.len() == batch.capacity() {
                builder.add_all(&batch)?;
                batch.clear();
                report(progress, 0.5 * reader.fraction(), &format!("Reading {input_name}"))?;
            }
        }
        builder.add_all(&batch)?;
        let scaled_progress = |f: f64, m: &str| progress.is_none_or(|p| p(0.5 + 0.5 * f, m));
        return builder.finish(Some(&scaled_progress));
    }
}

/// XYZ and extended XYZ (first frame unless `frame` is given).
pub fn build_from_xyz(input: &Path, output: &Path, frame: usize, progress: Progress, options: BuildOptions) -> Result<CloudInfo, String> {
    // Pass 1: the frame's bounds (XYZ has no box); pass 2: the atoms.
    let (mut lower, mut upper) = ([f64::MAX; 3], [f64::MIN; 3]);
    let mut atoms = 0u64;
    let mut comment = String::new();
    let input_name = file_name(input);
    let scan = |using: &mut dyn FnMut(&[u8], &[f64; 3]) -> Result<(), String>, share: f64, base: f64, atoms: &mut u64, comment: &mut String| -> Result<(), String> {
        let mut reader = LineReader::new(input)?;
        let mut current = 0usize;
        loop {
            let count_line = reader.line()?;
            *atoms = stoull(&count_line)?;
            *comment = reader.line()?;
            // Extended XYZ: species and position columns from Properties=.
            let (mut species, mut position) = (0usize, 1usize);
            if let Some(at) = comment.find("Properties=") {
                let rest = &comment[at + 11..];
                let spec = &rest[..rest.find([' ', '\t']).unwrap_or(rest.len())];
                let mut parts: Vec<&str> = spec.split(':').collect();
                if spec.is_empty() || spec.ends_with(':') {
                    parts.pop();
                }
                let mut column = 0usize;
                let mut i = 0;
                while i + 2 < parts.len() {
                    let width = stoi_prefix(parts[i + 2])?;
                    if parts[i] == "species" {
                        species = column;
                    }
                    if parts[i] == "pos" {
                        position = column;
                    }
                    column = (column as i64 + i64::from(width)) as usize;
                    i += 3;
                }
            }
            if current < frame {
                for _ in 0..*atoms {
                    reader.line()?;
                }
                current += 1;
                continue;
            }
            let needed = species.max(position + 2);
            for i in 0..*atoms {
                let Some(line) = reader.next() else {
                    return Err("The XYZ file ends early".into());
                };
                let line = line.to_vec();
                let f = fields(&line, needed + 1);
                let mut r = [0.0f64; 3];
                for k in 0..3 {
                    let (a, b) = f[position + k];
                    r[k] = parse_double(&line[a..b]).ok_or_else(|| format!("Bad coordinate in the XYZ file: {}", String::from_utf8_lossy(&line)))?;
                }
                let (a, b) = f[species];
                using(&line[a..b], &r)?;
                if i & 0xFFFF == 0 {
                    report(progress, base + share * reader.fraction(), &format!("Reading {input_name}"))?;
                }
            }
            return Ok(());
        }
    };
    scan(
        &mut |_, r| {
            for k in 0..3 {
                if r[k] < lower[k] {
                    lower[k] = r[k];
                }
                if upper[k] < r[k] {
                    upper[k] = r[k];
                }
            }
            Ok(())
        },
        0.25,
        0.0,
        &mut atoms,
        &mut comment,
    )?;
    let mut builder = CloudBuilder::new(output, lower, upper, atoms, options)?;
    builder.set_source(&format!("XYZ {input_name}, frame {frame}"));
    // Extended XYZ lattice.
    if let Some(at) = comment.find("Lattice=\"") {
        let rest = &comment[at + 9..];
        let values = &rest[..rest.find('"').unwrap_or(rest.len())];
        let numbers: Vec<f64> = words(values).iter().map_while(|w| stod(w).ok().filter(|&(_, used)| used == w.len()).map(|(v, _)| v)).collect();
        if numbers.len() >= 9 {
            builder.set_cell([[numbers[0], numbers[1], numbers[2]], [numbers[3], numbers[4], numbers[5]], [numbers[6], numbers[7], numbers[8]]], [0.0; 3]);
        }
    }
    let mut names: Vec<(Vec<u8>, u16)> = Vec::new();
    let mut second_atoms = 0u64;
    let mut second_comment = String::new();
    scan(
        &mut |symbol, r| {
            let s = match names.iter().find(|(n, _)| n.as_slice() == symbol) {
                Some(&(_, index)) => index,
                None => {
                    let index = builder.species(&String::from_utf8_lossy(symbol))?;
                    names.push((symbol.to_vec(), index));
                    index
                }
            };
            builder.add(&SourceAtom { x: r[0], y: r[1], z: r[2], species: s })
        },
        0.25,
        0.25,
        &mut second_atoms,
        &mut second_comment,
    )?;
    let scaled_progress = |f: f64, m: &str| progress.is_none_or(|p| p(0.5 + 0.5 * f, m));
    builder.finish(Some(&scaled_progress))
}

/// A crystal of nx x ny x nz cells: lattice fcc, bcc, sc, diamond or hcp with
/// lattice constant a (c for hcp, ideal when 0) and one element (or two,
/// alternating on the basis, with `second`).
#[allow(clippy::too_many_arguments)]
pub fn generate_crystal(output: &Path, lattice: &str, a: f64, c: f64, cells: [u64; 3], element: &str, second: &str, progress: Progress, options: BuildOptions) -> Result<CloudInfo, String> {
    if a <= 0.0 {
        return Err("The lattice constant must be positive".into());
    }
    if cells.iter().any(|&n| n < 1) {
        return Err("Cell counts must be at least 1".into());
    }
    // Cell vectors and basis (fractional).
    let mut cell = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]];
    let basis: Vec<[f64; 3]> = match lattice {
        "fcc" => vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]],
        "bcc" => vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        "sc" => vec![[0.0, 0.0, 0.0]],
        "diamond" => vec![
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.75, 0.25, 0.75],
            [0.25, 0.75, 0.75],
        ],
        "hcp" => {
            // Orthohexagonal cell (a, sqrt(3) a, c) with four atoms.
            let cc = if c > 0.0 { c } else { a * (8.0f64 / 3.0).sqrt() };
            cell = [[a, 0.0, 0.0], [0.0, 3f64.sqrt() * a, 0.0], [0.0, 0.0, cc]];
            vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 1.0 / 6.0, 0.5], [0.0, 2.0 / 3.0, 0.5]]
        }
        _ => return Err(format!("Unknown lattice {lattice} (fcc, bcc, sc, diamond or hcp)")),
    };
    let total = cells[0] as f64 * cells[1] as f64 * cells[2] as f64 * basis.len() as f64;
    if total > 2.0e10 {
        return Err("More than 20 billion atoms".into());
    }
    let sup: [[f64; 3]; 3] = std::array::from_fn(|r| std::array::from_fn(|k| cell[r][k] * cells[r] as f64));
    let upper = [sup[0][0], sup[1][1], sup[2][2]];
    let mut builder = CloudBuilder::new(output, [0.0; 3], upper, total as u64, options)?;
    builder.set_cell(sup, [0.0; 3]);
    let pair = if second.is_empty() { String::new() } else { format!("/{second}") };
    builder.set_source(&format!("{lattice} {element}{pair}, a = {}, {} x {} x {} cells", format_g(a, 6), cells[0], cells[1], cells[2]));
    let first = builder.species(element)?;
    let other = if second.is_empty() { first } else { builder.species(second)? };
    let mut batch: Vec<SourceAtom> = Vec::with_capacity(1 << 16);
    for k in 0..cells[2] {
        for j in 0..cells[1] {
            for i in 0..cells[0] {
                for (b, base) in basis.iter().enumerate() {
                    let f = [i as f64 + base[0], j as f64 + base[1], k as f64 + base[2]];
                    batch.push(SourceAtom {
                        x: f[0] * cell[0][0] + f[1] * cell[1][0] + f[2] * cell[2][0],
                        y: f[0] * cell[0][1] + f[1] * cell[1][1] + f[2] * cell[2][1],
                        z: f[0] * cell[0][2] + f[1] * cell[1][2] + f[2] * cell[2][2],
                        species: if b % 2 == 1 { other } else { first },
                    });
                    if batch.len() == batch.capacity() {
                        builder.add_all(&batch)?;
                        batch.clear();
                    }
                }
            }
        }
        report(progress, 0.5 * (k + 1) as f64 / cells[2] as f64, "Generating atoms")?;
    }
    builder.add_all(&batch)?;
    let scaled_progress = |f: f64, m: &str| progress.is_none_or(|p| p(0.5 + 0.5 * f, m));
    builder.finish(Some(&scaled_progress))
}
