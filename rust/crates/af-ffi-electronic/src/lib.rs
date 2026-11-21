//! C interface of the electronic-structure post-processing for the Python
//! package (atomforge.electronic, via ctypes): the C++ electronic/PythonAPI,
//! BaderAPI, LobsterAPI and WannierAPI under the same names. Errors never
//! cross the boundary: a failing call returns null (or 0) and the message is
//! read with the matching *_error function. Each returned handle is owned by
//! the caller; getters expose borrowed data until the handle is freed.

#![allow(clippy::missing_safety_doc)]

use std::cell::RefCell;
use std::ffi::{c_char, c_int, c_void, CStr, CString};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;

use af_electronic::bader::{bader_on_grid, BaderPartition};
use af_electronic::charge::*;
use af_electronic::density::*;
use af_electronic::field::*;
use af_electronic::fourier::*;
use af_electronic::lobster::{read_cohpcar, read_icohplist, CohpData, IcohpList};
use af_electronic::topology::betti_curve;
use af_electronic::volume::{arithmetic, minimum_image, periodic_endpoints, resample, scale};
use af_electronic::wannier::{berry_curvature, bloch_hamiltonian, chern_number, read_wannier_hr, wannier_bands, WannierHamiltonian};
use af_electronic::{isosurface, load_volume, save_mesh, save_volume, Grid, Mesh, Site, Volume};
use glam::{DMat3, DVec3, IVec3};

// ------------------------------------------------------------------ errors

thread_local! {
    static ELECTRONIC_ERROR: RefCell<CString> = RefCell::new(CString::default());
    static BADER_ERROR: RefCell<CString> = RefCell::new(CString::default());
    static LOBSTER_ERROR: RefCell<CString> = RefCell::new(CString::default());
    static WANNIER_ERROR: RefCell<CString> = RefCell::new(CString::default());
}

type ErrorSlot = std::thread::LocalKey<RefCell<CString>>;

fn c_string(text: &str) -> CString {
    CString::new(text.replace('\0', "?")).unwrap_or_default()
}

/// Runs `work` with the slot cleared; on failure (error or panic) stores the message.
fn guarded<T>(slot: &'static ErrorSlot, unknown: &str, work: impl FnOnce() -> Result<T, String>) -> Option<T> {
    slot.with(|e| *e.borrow_mut() = CString::default());
    let message = match catch_unwind(AssertUnwindSafe(work)) {
        Ok(Ok(value)) => return Some(value),
        Ok(Err(message)) => message,
        Err(_) => unknown.to_string(),
    };
    slot.with(|e| *e.borrow_mut() = c_string(&message));
    None
}

fn error_pointer(slot: &'static ErrorSlot) -> *const c_char {
    slot.with(|e| e.borrow().as_ptr())
}

fn text_arg(p: *const c_char) -> Option<String> {
    if p.is_null() {
        None
    } else {
        Some(unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned())
    }
}

fn into_handle<T>(value: T) -> *mut c_void {
    Box::into_raw(Box::new(value)) as *mut c_void
}

fn borrow<'a, T>(h: *mut c_void, what: &str) -> Result<&'a mut T, String> {
    if h.is_null() {
        Err(what.to_string())
    } else {
        Ok(unsafe { &mut *(h as *mut T) })
    }
}

/// glm cell from a row-major buffer whose row a is cell vector a (column a of the glm matrix).
unsafe fn read_grid(shape: *const c_int, cell: *const f64, origin: *const f64, periodic: c_int, values: *const f64, count: usize, unit: *const c_char) -> Grid {
    let shape = std::slice::from_raw_parts(shape, 3);
    let c = std::slice::from_raw_parts(cell, 9);
    let o = std::slice::from_raw_parts(origin, 3);
    Grid {
        shape: [shape[0], shape[1], shape[2]],
        cell: DMat3::from_cols(DVec3::new(c[0], c[1], c[2]), DVec3::new(c[3], c[4], c[5]), DVec3::new(c[6], c[7], c[8])),
        origin: DVec3::new(o[0], o[1], o[2]),
        periodic: periodic != 0,
        unit: text_arg(unit).unwrap_or_default(),
        values: std::slice::from_raw_parts(values, count).to_vec(),
        ..Default::default()
    }
}

/// std::vector::at (libstdc++'s out_of_range message; a negative int index
/// converts to a huge size_t in C++).
fn at<T>(v: &[T], index: i64) -> Result<&T, String> {
    if index < 0 || index as usize >= v.len() {
        return Err(format!("vector::_M_range_check: __n (which is {}) >= this->size() (which is {})", index as u64, v.len()));
    }
    Ok(&v[index as usize])
}

// ------------------------------------------------------------------ volumes (PythonAPI)

#[derive(Default)]
struct Handle {
    volume: Volume,
    table: Vec<f64>,
    columns: c_int,
    /// (name, unit) of each field, NUL-terminated for af_electronic_label.
    labels: Vec<(CString, CString)>,
}

impl Handle {
    fn finish(mut self) -> *mut c_void {
        self.labels = self.volume.fields.iter().map(|g| (c_string(&g.name), c_string(&g.unit))).collect();
        into_handle(self)
    }
    fn clone_data(&self) -> Handle {
        Handle { volume: self.volume.clone(), table: self.table.clone(), columns: self.columns, labels: Vec::new() }
    }
}

fn result<'a>(h: *mut c_void) -> Result<&'a mut Handle, String> {
    borrow(h, "Null result handle")
}

fn grid<'a>(h: *mut c_void) -> Result<&'a Grid, String> {
    result(h)?.volume.fields.first().ok_or_else(|| "Operation requires a grid".to_string())
}

fn exact_int(value: f64) -> Result<i32, String> {
    if !value.is_finite() || value != value.floor() || value.abs() > 100000000.0 {
        return Err("Expected a bounded integer parameter".into());
    }
    Ok(value as i32)
}

fn electronic<T>(work: impl FnOnce() -> Result<T, String>) -> Option<T> {
    guarded(&ELECTRONIC_ERROR, "Unknown native post-processing error", work)
}

#[no_mangle]
pub extern "C" fn af_electronic_version() -> c_int {
    1
}

#[no_mangle]
pub extern "C" fn af_electronic_error() -> *const c_char {
    error_pointer(&ELECTRONIC_ERROR)
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_free(h: *mut c_void) {
    if !h.is_null() {
        drop(Box::from_raw(h as *mut Handle));
    }
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_load(path: *const c_char, quantity: *const c_char, coordinate_unit: *const c_char) -> *mut c_void {
    electronic(|| {
        let (Some(path), Some(quantity), Some(unit)) = (text_arg(path), text_arg(quantity), text_arg(coordinate_unit)) else {
            return Err("Null load argument".into());
        };
        Ok(Handle { volume: load_volume(&path, &quantity, &unit)?, ..Default::default() }.finish())
    })
    .unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_create(
    shape: *const c_int,
    cell: *const f64,
    origin: *const f64,
    periodic: c_int,
    values: *const f64,
    count: usize,
    unit: *const c_char,
) -> *mut c_void {
    electronic(|| {
        if shape.is_null() || cell.is_null() || origin.is_null() || values.is_null() || unit.is_null() || count > 100000000 {
            return Err("Invalid grid buffers".into());
        }
        let g = read_grid(shape, cell, origin, periodic, values, count, unit);
        g.validate()?;
        let mut out = Handle::default();
        out.volume.fields.push(g);
        Ok(out.finish())
    })
    .unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub extern "C" fn af_electronic_select(h: *mut c_void, index: c_int) -> *mut c_void {
    electronic(|| {
        let r = result(h)?;
        let mut out = Handle::default();
        out.volume.sites = r.volume.sites.clone();
        out.volume.fields.push(at(&r.volume.fields, i64::from(index))?.clone());
        Ok(out.finish())
    })
    .unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub extern "C" fn af_electronic_fields(h: *mut c_void) -> c_int {
    if h.is_null() {
        0
    } else {
        unsafe { &*(h as *const Handle) }.volume.fields.len() as c_int
    }
}

#[no_mangle]
pub extern "C" fn af_electronic_label(h: *mut c_void, field: c_int, unit: c_int) -> *const c_char {
    if h.is_null() || field < 0 || field >= af_electronic_fields(h) {
        return c"".as_ptr();
    }
    let (name, u) = &unsafe { &*(h as *const Handle) }.labels[field as usize];
    if unit != 0 {
        u.as_ptr()
    } else {
        name.as_ptr()
    }
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_data(h: *mut c_void, destination: *mut f64) -> usize {
    if h.is_null() {
        return 0;
    }
    let r = &*(h as *const Handle);
    let data = match r.volume.fields.first() {
        Some(g) => &g.values,
        None => &r.table,
    };
    if !destination.is_null() {
        std::ptr::copy_nonoverlapping(data.as_ptr(), destination, data.len());
    }
    data.len()
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_geometry(h: *mut c_void, shape: *mut c_int, cell: *mut f64, origin: *mut f64) -> c_int {
    if h.is_null() {
        return -1;
    }
    let Some(g) = (*(h as *const Handle)).volume.fields.first() else { return -1 };
    for a in 0..3 {
        if !shape.is_null() {
            *shape.add(a) = g.shape[a];
        }
        if !origin.is_null() {
            *origin.add(a) = g.origin[a];
        }
        for b in 0..3 {
            if !cell.is_null() {
                *cell.add(a * 3 + b) = g.cell.col(a)[b];
            }
        }
    }
    c_int::from(g.periodic)
}

#[no_mangle]
pub extern "C" fn af_electronic_columns(h: *mut c_void) -> c_int {
    if h.is_null() {
        0
    } else {
        unsafe { &*(h as *const Handle) }.columns
    }
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_sites(h: *mut c_void, destination: *mut f64) -> usize {
    if h.is_null() {
        return 0;
    }
    let sites = &(*(h as *const Handle)).volume.sites;
    if !destination.is_null() {
        for (i, s) in sites.iter().enumerate() {
            *destination.add(i * 4) = f64::from(s.number);
            for a in 0..3 {
                *destination.add(i * 4 + a + 1) = s.position[a];
            }
        }
    }
    sites.len()
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_with_sites(h: *mut c_void, sites: *const f64, count: usize) -> *mut c_void {
    electronic(|| {
        if count > 10000000 || (count > 0 && sites.is_null()) {
            return Err("Invalid site buffer".into());
        }
        let mut out = result(h)?.clone_data();
        out.volume.sites.clear();
        let data = if count > 0 { std::slice::from_raw_parts(sites, count * 4) } else { &[] };
        for i in 0..count {
            let z = exact_int(data[i * 4])?;
            if !(0..=118).contains(&z) {
                return Err("Atomic number outside 0..118".into());
            }
            let p = DVec3::new(data[i * 4 + 1], data[i * 4 + 2], data[i * 4 + 3]);
            if !p.is_finite() {
                return Err("Non-finite site".into());
            }
            out.volume.sites.push(Site { number: z, position: p });
        }
        Ok(out.finish())
    })
    .unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_save(h: *mut c_void, path: *const c_char, format: *const c_char) -> c_int {
    let done = electronic(|| {
        let (Some(path), Some(format)) = (text_arg(path), text_arg(format)) else {
            return Err("Null output argument".into());
        };
        let r = result(h)?;
        if format == "obj" {
            if r.columns != 4 {
                return Err("OBJ export requires an isosurface result".into());
            }
            let mut m = Mesh::default();
            for row in r.table.chunks_exact(4) {
                m.vertices.push(DVec3::new(row[0], row[1], row[2]));
                m.colors.push(row[3]);
            }
            save_mesh(&m, Path::new(&path))
        } else {
            save_volume(&r.volume, &path, &format)
        }
    });
    c_int::from(done.is_some())
}

fn calculate(h: *mut c_void, other: *mut c_void, op: &str, parameters: &[f64]) -> Result<Handle, String> {
    if parameters.iter().any(|v| !v.is_finite()) {
        return Err("Non-finite calculation parameter".into());
    }
    let g = grid(h)?;
    let sites = result(h)?.volume.sites.clone();
    let mut out = Handle::default();
    out.volume.sites = sites.clone();
    let count = parameters.len();
    let need = |n: usize| if count != n { Err(format!("Incorrect parameter count for {op}")) } else { Ok(()) };
    let p = |i: usize| parameters.get(i).copied().ok_or_else(|| format!("Missing parameter for {op}"));
    let vec = |i: usize| -> Result<DVec3, String> { Ok(DVec3::new(p(i)?, p(i + 1)?, p(i + 2)?)) };
    let ivec = |i: usize| -> Result<IVec3, String> { Ok(IVec3::new(exact_int(p(i)?)?, exact_int(p(i + 1)?)?, exact_int(p(i + 2)?)?)) };
    let mut fields: Vec<Grid> = Vec::new();
    let profile = |out: &mut Handle, rows: Vec<glam::DVec2>| {
        out.columns = 2;
        for row in rows {
            out.table.push(row.x);
            out.table.push(row.y);
        }
    };
    let reflections = |out: &mut Handle, list: Vec<Reflection>| {
        out.columns = 5;
        for r in list {
            out.table.extend([f64::from(r.hkl.x), f64::from(r.hkl.y), f64::from(r.hkl.z), r.value.re, r.value.im]);
        }
    };
    match op {
        "add" | "subtract" | "multiply" | "divide" => {
            need(0)?;
            fields.push(arithmetic(g, grid(other)?, op)?);
        }
        "scale" => {
            need(1)?;
            fields.push(scale(g, p(0)?)?);
        }
        "density_difference" => {
            need(1)?;
            fields.push(density_difference(g, grid(other)?, p(0)?)?);
        }
        "threshold_mask" => {
            need(2)?;
            fields.push(threshold_mask(g, p(0)?, p(1)?)?);
        }
        "union" | "intersection" | "difference" | "xor" => {
            need(0)?;
            fields.push(boolean_mask(g, grid(other)?, op)?);
        }
        "apply_mask" => {
            need(0)?;
            fields.push(apply_mask(g, grid(other)?)?);
        }
        "split_density" => {
            need(0)?;
            fields.extend(split_density(g)?);
        }
        "charge_summary" => {
            need(0)?;
            out.columns = 3;
            out.table = charge_summary(g)?.to_vec();
        }
        "cumulative_charge" => {
            need(1)?;
            profile(&mut out, cumulative_charge(g, exact_int(p(0)?)?)?);
        }
        "resample" => {
            need(0)?;
            fields.push(resample(g, grid(other)?)?);
        }
        "periodic" => {
            need(1)?;
            fields.push(periodic_endpoints(g, p(0)?)?);
        }
        "smooth" => {
            need(2)?;
            fields.push(smooth(g, p(0)?, exact_int(p(1)?)?)?);
        }
        "gradient" => {
            need(0)?;
            fields.extend(gradient(g)?);
        }
        "laplacian" => {
            need(0)?;
            fields.push(laplacian(g)?);
        }
        "energy" => {
            need(1)?;
            fields.extend(energy_density(g, p(0)?)?);
        }
        "reduced_density_gradient" => {
            need(1)?;
            fields.push(reduced_density_gradient(g, p(0)?)?);
        }
        "signed_density" => {
            need(1)?;
            fields.push(signed_density(g, p(0)?)?);
        }
        "dori" => {
            need(1)?;
            fields.push(dori(g, p(0)?)?);
        }
        "line" => {
            need(7)?;
            profile(&mut out, line_profile(g, vec(0)?, vec(3)?, exact_int(p(6)?)?)?);
        }
        "planar" => {
            need(1)?;
            profile(&mut out, planar_average(g, exact_int(p(0)?)?)?);
        }
        "macro" => {
            need(2)?;
            profile(&mut out, macroscopic_average(g, exact_int(p(0)?)?, exact_int(p(1)?)?)?);
        }
        "section" => {
            need(11)?;
            fields.push(section(g, vec(0)?, vec(3)?, vec(6)?, exact_int(p(9)?)?, exact_int(p(10)?)?)?);
        }
        "contours" => {
            need(1)?;
            out.columns = 3;
            for v in contours(g, p(0)?)? {
                out.table.extend([v.x, v.y, v.z]);
            }
        }
        "integrate" => {
            need(0)?;
            out.columns = 1;
            out.table = vec![integrate(g)?];
        }
        "sphere" => {
            need(4)?;
            out.columns = 1;
            out.table = vec![integrate_sphere(g, vec(0)?, p(3)?)?];
        }
        "voronoi" => {
            if count % 3 != 0 {
                return Err("Expected XYZ sites".into());
            }
            let points: Vec<DVec3> = (0..count).step_by(3).map(vec).collect::<Result<_, _>>()?;
            out.columns = 2;
            for b in voronoi_integrate(g, &points)? {
                out.table.extend([b.integral, b.volume]);
            }
        }
        "peaks" => {
            need(1)?;
            out.columns = 4;
            for peak in peaks(g, exact_int(p(0)?)?)? {
                out.table.extend([peak.position.x, peak.position.y, peak.position.z, peak.value]);
            }
        }
        "patterson" => {
            need(0)?;
            fields.push(patterson(g)?);
        }
        "betti_curve" => {
            out.columns = 4;
            for row in betti_curve(g, parameters)? {
                out.table.extend([row.x, row.y, row.z, row.w]);
            }
        }
        "factors" => {
            if count % 3 != 0 {
                return Err("Expected hkl triples".into());
            }
            let indices: Vec<IVec3> = (0..count).step_by(3).map(ivec).collect::<Result<_, _>>()?;
            reflections(&mut out, structure_factors(g, &indices)?);
        }
        "synthesis" => {
            if count % 5 != 0 {
                return Err("Expected h k l real imag rows".into());
            }
            let mut list = Vec::new();
            for i in (0..count).step_by(5) {
                list.push(Reflection { hkl: ivec(i)?, value: Complex::new(p(i + 3)?, p(i + 4)?) });
            }
            fields.push(fourier_synthesis(g, &list)?);
        }
        "atomic_factors" => {
            let n = sites.len();
            if n == 0 || count < n * 2 || (count - n * 2) % (3 + n * 2) != 0 {
                return Err("Invalid atomic-factor parameter layout".into());
            }
            let positions: Vec<DVec3> = sites.iter().map(|s| s.position).collect();
            let (mut occupancy, mut b) = (Vec::new(), Vec::new());
            for i in 0..n {
                occupancy.push(p(i * 2)?);
                b.push(p(i * 2 + 1)?);
            }
            let (mut indices, mut factors) = (Vec::new(), Vec::new());
            let mut i = n * 2;
            while i < count {
                indices.push(ivec(i)?);
                for j in 0..n {
                    factors.push(Complex::new(p(i + 3 + j * 2)?, p(i + 4 + j * 2)?));
                }
                i += 3 + n * 2;
            }
            reflections(&mut out, atomic_structure_factors(g, &positions, &occupancy, &b, &indices, &factors)?);
        }
        "ewald" => {
            need(3 + sites.len())?;
            let positions: Vec<DVec3> = sites.iter().map(|s| s.position).collect();
            let charges: Vec<f64> = (0..sites.len()).map(|i| p(i + 3)).collect::<Result<_, _>>()?;
            let e = ewald(&g.cell, &positions, &charges, p(0)?, p(1)?, p(2)?)?;
            out.columns = 1;
            out.table = vec![e.energy];
            out.table.extend(e.potentials);
        }
        "isosurface" => {
            need(1)?;
            out.columns = 4;
            let colour = if other.is_null() { None } else { Some(grid(other)?) };
            let mesh = isosurface(g, p(0)?, colour)?;
            for (v, c) in mesh.vertices.iter().zip(&mesh.colors) {
                out.table.extend([v.x, v.y, v.z, *c]);
            }
        }
        _ => return Err(format!("Unknown electronic operation: {op}")),
    }
    out.volume.fields = fields;
    Ok(out)
}

#[no_mangle]
pub unsafe extern "C" fn af_electronic_calculate(h: *mut c_void, other: *mut c_void, operation: *const c_char, parameters: *const f64, count: usize) -> *mut c_void {
    electronic(|| {
        let Some(op) = text_arg(operation) else { return Err("Null calculation arguments".into()) };
        if count > 0 && parameters.is_null() {
            return Err("Null calculation arguments".into());
        }
        let values = if count > 0 { std::slice::from_raw_parts(parameters, count) } else { &[] };
        Ok(calculate(h, other, &op, values)?.finish())
    })
    .unwrap_or(std::ptr::null_mut())
}

// ------------------------------------------------------------------ Bader (BaderAPI)

struct BaderHandle {
    partition: BaderPartition,
    cell: DMat3,
    periodic: bool,
}

fn bader<T>(work: impl FnOnce() -> Result<T, String>) -> Option<T> {
    guarded(&BADER_ERROR, "Unknown native Bader error", work)
}

fn bader_handle<'a>(h: *mut c_void) -> Result<&'a mut BaderHandle, String> {
    borrow(h, "Null Bader partition handle")
}

#[no_mangle]
pub extern "C" fn af_bader_error() -> *const c_char {
    error_pointer(&BADER_ERROR)
}

#[no_mangle]
pub unsafe extern "C" fn af_bader_free(h: *mut c_void) {
    if !h.is_null() {
        drop(Box::from_raw(h as *mut BaderHandle));
    }
}

#[no_mangle]
pub unsafe extern "C" fn af_bader_partition(
    shape: *const c_int,
    cell: *const f64,
    origin: *const f64,
    periodic: c_int,
    values: *const f64,
    count: usize,
    unit: *const c_char,
) -> *mut c_void {
    bader(|| {
        if shape.is_null() || cell.is_null() || origin.is_null() || values.is_null() || unit.is_null() {
            return Err("Invalid grid buffers".into());
        }
        let g = read_grid(shape, cell, origin, periodic, values, count, unit);
        Ok(into_handle(BaderHandle { cell: g.cell, periodic: g.periodic, partition: bader_on_grid(&g)? }))
    })
    .unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub extern "C" fn af_bader_num_basins(h: *mut c_void) -> c_int {
    if h.is_null() {
        0
    } else {
        unsafe { &*(h as *const BaderHandle) }.partition.basins.len() as c_int
    }
}

/// basin_ids: one double per voxel, x-fastest order, matching the source grid.
#[no_mangle]
pub unsafe extern "C" fn af_bader_basin_ids(h: *mut c_void, basin_ids: *mut f64) -> c_int {
    c_int::from(
        bader(|| {
            if basin_ids.is_null() {
                return Err("Null buffer".into());
            }
            for (i, &id) in bader_handle(h)?.partition.basin.iter().enumerate() {
                *basin_ids.add(i) = f64::from(id);
            }
            Ok(())
        })
        .is_some(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn af_bader_basin(h: *mut c_void, index: c_int, charge: *mut f64, volume: *mut f64, maximum: *mut f64) -> c_int {
    c_int::from(
        bader(|| {
            if charge.is_null() || volume.is_null() || maximum.is_null() {
                return Err("Null buffer".into());
            }
            let basin = at(&bader_handle(h)?.partition.basins, i64::from(index))?;
            *charge = basin.charge;
            *volume = basin.volume;
            for a in 0..3 {
                *maximum.add(a) = basin.maximum[a];
            }
            Ok(())
        })
        .is_some(),
    )
}

/// sites: site_count*3 Cartesian doubles; populations: site_count doubles out.
#[no_mangle]
pub unsafe extern "C" fn af_bader_populations(h: *mut c_void, sites: *const f64, site_count: usize, populations: *mut f64) -> c_int {
    c_int::from(
        bader(|| {
            if sites.is_null() || populations.is_null() || site_count == 0 {
                return Err("Invalid site buffer".into());
            }
            let info = bader_handle(h)?;
            let s = std::slice::from_raw_parts(sites, site_count * 3);
            let points: Vec<DVec3> = s.chunks_exact(3).map(|c| DVec3::new(c[0], c[1], c[2])).collect();
            let out = std::slice::from_raw_parts_mut(populations, site_count);
            out.fill(0.0);
            for basin in &info.partition.basins {
                let mut nearest = 0;
                let mut best = f64::INFINITY;
                for (k, &point) in points.iter().enumerate() {
                    let mut delta = basin.maximum - point;
                    if info.periodic {
                        delta = minimum_image(delta, &info.cell)?;
                    }
                    let d2 = delta.dot(delta);
                    if d2 < best {
                        best = d2;
                        nearest = k;
                    }
                }
                out[nearest] += basin.charge;
            }
            Ok(())
        })
        .is_some(),
    )
}

// ------------------------------------------------------------------ LOBSTER (LobsterAPI)

fn lobster<T>(work: impl FnOnce() -> Result<T, String>) -> Option<T> {
    guarded(&LOBSTER_ERROR, "Unknown native LOBSTER error", work)
}

fn lobster_ok(work: impl FnOnce() -> Result<(), String>) -> c_int {
    c_int::from(lobster(work).is_some())
}

fn cohpcar<'a>(h: *mut c_void) -> Result<&'a mut CohpData, String> {
    borrow(h, "Null COHPCAR handle")
}

fn icohplist<'a>(h: *mut c_void) -> Result<&'a mut IcohpList, String> {
    borrow(h, "Null ICOHPLIST handle")
}

unsafe fn copy_to(values: &[f64], buffer: *mut f64) {
    std::ptr::copy_nonoverlapping(values.as_ptr(), buffer, values.len());
}

#[no_mangle]
pub extern "C" fn af_lobster_error() -> *const c_char {
    error_pointer(&LOBSTER_ERROR)
}

#[no_mangle]
pub unsafe extern "C" fn af_cohpcar_load(path: *const c_char) -> *mut c_void {
    lobster(|| Ok(into_handle(read_cohpcar(&text_arg(path).ok_or("Null path")?)?))).unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub unsafe extern "C" fn af_cohpcar_free(h: *mut c_void) {
    if !h.is_null() {
        drop(Box::from_raw(h as *mut CohpData));
    }
}

#[no_mangle]
pub unsafe extern "C" fn af_cohpcar_info(h: *mut c_void, spin_polarized: *mut c_int, num_energies: *mut usize, num_bonds: *mut c_int) -> c_int {
    lobster_ok(|| {
        if spin_polarized.is_null() || num_energies.is_null() || num_bonds.is_null() {
            return Err("Null buffer".into());
        }
        let data = cohpcar(h)?;
        *spin_polarized = c_int::from(data.spin_polarized);
        *num_energies = data.energies.len();
        *num_bonds = data.bonds.len() as c_int;
        Ok(())
    })
}

#[no_mangle]
pub extern "C" fn af_cohpcar_fermi_energy(h: *mut c_void) -> f64 {
    if h.is_null() {
        0.0
    } else {
        unsafe { &*(h as *const CohpData) }.fermi_energy
    }
}

#[no_mangle]
pub unsafe extern "C" fn af_cohpcar_energies(h: *mut c_void, buffer: *mut f64) -> c_int {
    lobster_ok(|| {
        if buffer.is_null() {
            return Err("Null buffer".into());
        }
        copy_to(&cohpcar(h)?.energies, buffer);
        Ok(())
    })
}

#[no_mangle]
pub unsafe extern "C" fn af_cohpcar_bond(h: *mut c_void, bond: c_int, atom1: *mut c_int, atom2: *mut c_int, length: *mut f64) -> c_int {
    lobster_ok(|| {
        if atom1.is_null() || atom2.is_null() || length.is_null() {
            return Err("Null buffer".into());
        }
        let info = at(&cohpcar(h)?.bonds, i64::from(bond))?;
        *atom1 = info.atom1;
        *atom2 = info.atom2;
        *length = info.length;
        Ok(())
    })
}

/// cohp/icohp: num_energies doubles out each. spin: 0 (up/only channel) or 1 (down).
#[no_mangle]
pub unsafe extern "C" fn af_cohpcar_average(h: *mut c_void, spin: c_int, cohp: *mut f64, icohp: *mut f64) -> c_int {
    lobster_ok(|| {
        if cohp.is_null() || icohp.is_null() {
            return Err("Null buffer".into());
        }
        let data = cohpcar(h)?;
        copy_to(at(&data.average_cohp, i64::from(spin))?, cohp);
        copy_to(at(&data.average_icohp, i64::from(spin))?, icohp);
        Ok(())
    })
}

#[no_mangle]
pub unsafe extern "C" fn af_cohpcar_bond_curve(h: *mut c_void, bond: c_int, spin: c_int, cohp: *mut f64, icohp: *mut f64) -> c_int {
    lobster_ok(|| {
        if cohp.is_null() || icohp.is_null() {
            return Err("Null buffer".into());
        }
        let data = cohpcar(h)?;
        copy_to(at(at(&data.cohp, i64::from(spin))?, i64::from(bond))?, cohp);
        copy_to(at(at(&data.icohp, i64::from(spin))?, i64::from(bond))?, icohp);
        Ok(())
    })
}

#[no_mangle]
pub unsafe extern "C" fn af_icohplist_load(path: *const c_char) -> *mut c_void {
    lobster(|| Ok(into_handle(read_icohplist(&text_arg(path).ok_or("Null path")?)?))).unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub unsafe extern "C" fn af_icohplist_free(h: *mut c_void) {
    if !h.is_null() {
        drop(Box::from_raw(h as *mut IcohpList));
    }
}

#[no_mangle]
pub extern "C" fn af_icohplist_count(h: *mut c_void) -> c_int {
    if h.is_null() {
        0
    } else {
        unsafe { &*(h as *const IcohpList) }.entries.len() as c_int
    }
}

#[no_mangle]
pub extern "C" fn af_icohplist_spin_polarized(h: *mut c_void) -> c_int {
    c_int::from(!h.is_null() && unsafe { &*(h as *const IcohpList) }.spin_polarized)
}

#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn af_icohplist_entry(
    h: *mut c_void,
    index: c_int,
    atom1: *mut c_int,
    atom2: *mut c_int,
    length: *mut f64,
    num_bonds: *mut c_int,
    icohp_up: *mut f64,
    icohp_down: *mut f64,
) -> c_int {
    lobster_ok(|| {
        if atom1.is_null() || atom2.is_null() || length.is_null() || num_bonds.is_null() || icohp_up.is_null() || icohp_down.is_null() {
            return Err("Null buffer".into());
        }
        let e = at(&icohplist(h)?.entries, i64::from(index))?;
        *atom1 = e.atom1;
        *atom2 = e.atom2;
        *length = e.length;
        *num_bonds = e.num_bonds;
        *icohp_up = e.icohp[0];
        *icohp_down = e.icohp[1];
        Ok(())
    })
}

// ------------------------------------------------------------------ Wannier (WannierAPI)

fn wannier_ok(work: impl FnOnce() -> Result<(), String>) -> c_int {
    c_int::from(guarded(&WANNIER_ERROR, "Unknown native Wannier error", work).is_some())
}

fn model<'a>(h: *mut c_void) -> Result<&'a mut WannierHamiltonian, String> {
    borrow(h, "Null Wannier model handle")
}

unsafe fn kpoint_of(k: *const f64) -> DVec3 {
    DVec3::new(*k, *k.add(1), *k.add(2))
}

#[no_mangle]
pub extern "C" fn af_wannier_error() -> *const c_char {
    error_pointer(&WANNIER_ERROR)
}

#[no_mangle]
pub unsafe extern "C" fn af_wannier_free(h: *mut c_void) {
    if !h.is_null() {
        drop(Box::from_raw(h as *mut WannierHamiltonian));
    }
}

#[no_mangle]
pub unsafe extern "C" fn af_wannier_load(path: *const c_char) -> *mut c_void {
    guarded(&WANNIER_ERROR, "Unknown native Wannier error", || Ok(into_handle(read_wannier_hr(&text_arg(path).ok_or("Null path")?)?)))
        .unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub extern "C" fn af_wannier_num_wann(h: *mut c_void) -> c_int {
    if h.is_null() {
        0
    } else {
        unsafe { &*(h as *const WannierHamiltonian) }.num_wann
    }
}

/// matrix: num_wann*num_wann interleaved (real, imag) pairs out, row-major H(k).
#[no_mangle]
pub unsafe extern "C" fn af_wannier_hamiltonian(h: *mut c_void, kpoint: *const f64, matrix: *mut f64) -> c_int {
    wannier_ok(|| {
        if kpoint.is_null() || matrix.is_null() {
            return Err("Null buffer".into());
        }
        for (i, v) in bloch_hamiltonian(model(h)?, kpoint_of(kpoint))?.iter().enumerate() {
            *matrix.add(2 * i) = v.re;
            *matrix.add(2 * i + 1) = v.im;
        }
        Ok(())
    })
}

/// kpoints: count*3 fractional-reciprocal doubles; energies: count*num_wann doubles out, ascending per k point.
#[no_mangle]
pub unsafe extern "C" fn af_wannier_bands(h: *mut c_void, kpoints: *const f64, count: usize, energies: *mut f64) -> c_int {
    wannier_ok(|| {
        if kpoints.is_null() || energies.is_null() {
            return Err("Null buffer".into());
        }
        let points: Vec<DVec3> = (0..count).map(|i| kpoint_of(kpoints.add(i * 3))).collect();
        copy_to(&wannier_bands(model(h)?, &points)?, energies);
        Ok(())
    })
}

/// curvature: num_wann doubles out, ascending band order, Hamiltonian-gauge (Kubo) term.
#[no_mangle]
pub unsafe extern "C" fn af_wannier_berry_curvature(h: *mut c_void, kpoint: *const f64, plane_a: c_int, plane_b: c_int, curvature: *mut f64) -> c_int {
    wannier_ok(|| {
        if kpoint.is_null() || curvature.is_null() {
            return Err("Null buffer".into());
        }
        copy_to(&berry_curvature(model(h)?, kpoint_of(kpoint), plane_a, plane_b)?, curvature);
        Ok(())
    })
}

#[no_mangle]
pub unsafe extern "C" fn af_wannier_chern_number(h: *mut c_void, band: c_int, plane_a: c_int, plane_b: c_int, grid: c_int, fixed: f64, result: *mut f64) -> c_int {
    wannier_ok(|| {
        if result.is_null() {
            return Err("Null buffer".into());
        }
        *result = chern_number(model(h)?, band, plane_a, plane_b, grid, fixed)?;
        Ok(())
    })
}
