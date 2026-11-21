//! Analysis, property and colouring steps (the C++ pipeline/AnalysisSteps, in
//! the spirit of OVITO's modifiers): per-atom analyses through the scientific
//! tools, statistics, binning and selections derived from them.

use std::path::Path;

use af_core::format::format_g;
use af_core::json::Json;
use af_core::Structure;
use af_science::atom_properties::{property_name, AtomProperty};
use af_science::core::{fractional, identity, neighbor_list};
use af_science::tools::{analyse_per_atom, atoms_request as science_atoms_request, periodic_axes, run_tool};

use crate::expression::Expression;
use crate::pipeline::{params, ModifierParameter, ModifierType, PipelineData};
use crate::support::*;

type R = Result<(), String>;

/// Positions, cell and periodic axes of the structure as a tool request.
fn atoms_request(s: &Structure) -> Json {
    science_atoms_request(s, &periodic_axes(s))
}

fn tool(name: &str, request: &Json) -> Result<Json, String> {
    Ok(run_tool(name, request, Path::new("."), None)?.result)
}

/// snprintf into a buffer of `size` bytes: at most size - 1 bytes.
fn truncated(mut text: String, size: usize) -> String {
    if text.len() >= size {
        let mut end = size - 1;
        while !text.is_char_boundary(end) {
            end -= 1;
        }
        text.truncate(end);
    }
    text
}

fn mask_of(list: &Json) -> Result<Vec<u8>, String> {
    list.items()?.iter().map(|f| Ok(f.boolean()? as u8)).collect()
}

/// Atoms on the surface of the solid: those a probe sphere touches.
fn select_surface(p: &Json, data: &mut PipelineData) -> R {
    let mut request = atoms_request(&data.structure);
    request.set("probe_radius_A", number(p, "select-surface", "probe_radius")?);
    request.set("export_mesh", false);
    let result = tool("surface-mesh", &request)?;
    let mask = mask_of(result.at("surface_atom")?)?;
    data.notes.push(truncated(
        format!(
            "Surface area {} A^2, solid volume {} A^3",
            format_g(result.at("surface_area_A2")?.number()?, 6),
            format_g(result.at("solid_volume_A3")?.number()?, 6)
        ),
        160,
    ));
    select(data, &mask, &text(p, "select-surface", "mode")?);
    Ok(())
}

/// Atoms inside LAMMPS regions (region commands, ';' between them).
fn select_region(p: &Json, data: &mut PipelineData) -> R {
    let mut request = atoms_request(&data.structure);
    request.set("regions", text(p, "select-region", "regions")?);
    request.set("lattice_scale", number(p, "select-region", "lattice_scale")?);
    request.set("export_mesh", false);
    let result = tool("lammps-regions", &request)?;
    let mask = mask_of(result.at("inside_any")?)?;
    select(data, &mask, &text(p, "select-region", "mode")?);
    Ok(())
}

/// Pairs of atoms closer than a distance: both atoms, or only the later one
/// of each pair (so deleting the selection keeps one atom per overlap).
fn select_overlapping(p: &Json, data: &mut PipelineData) -> R {
    let distance = number(p, "select-overlapping", "distance")?;
    if !(distance > 0.0) {
        return Err("select-overlapping needs a positive distance".into());
    }
    let both = text(p, "select-overlapping", "atoms")? == "both";
    let s = &data.structure;
    let positions: Vec<_> = s.atoms.iter().map(position_of).collect();
    let pbc = periodic_axes(s);
    let mut mask = vec![0u8; s.atoms.len()];
    let mut pairs = 0usize;
    let cell = if s.has_unit_cell { cell_of(s) } else { identity() };
    for n in neighbor_list(&positions, &cell, pbc, distance)? {
        if n.i >= n.j {
            continue; // each pair once; skips an atom's own periodic images
        }
        pairs += 1;
        mask[n.j] = 1;
        if both {
            mask[n.i] = 1;
        }
    }
    data.notes.push(format!("{pairs} overlapping pairs"));
    select(data, &mask, &text(p, "select-overlapping", "mode")?);
    Ok(())
}

fn set_property(data: &mut PipelineData, values: Vec<f64>, name: &str) {
    let valid = values.iter().filter(|v| v.is_finite()).count();
    let count = values.len();
    data.structure.atom_property = values;
    data.structure.atom_property_name = name.to_string();
    data.notes.push(format!("{name} computed for {valid} of {count} atoms"));
}

fn numbers_of(list: &Json) -> Result<Vec<f64>, String> {
    Ok(list.items()?.iter().map(|item| if item.is_number() { item.number().unwrap() } else { f64::NAN }).collect())
}

/// Elastic strain relative to an ideal lattice (polyhedral template matching).
fn elastic_strain_step(p: &Json, data: &mut PipelineData) -> R {
    let id = "elastic-strain";
    let mut request = atoms_request(&data.structure);
    request.set("lattice", text(p, id, "lattice")?);
    request.set("lattice_constant_A", number(p, id, "a")?);
    request.set("c_over_a", number(p, id, "c_over_a")?);
    let result = tool("elastic-strain", &request)?;
    let component = text(p, id, "component")?;
    if component == "volumetric" {
        set_property(data, numbers_of(result.at("volumetric_strain")?)?, "Volumetric elastic strain");
        return Ok(());
    }
    if component == "shear" {
        set_property(data, numbers_of(result.at("shear_strain")?)?, "Von Mises elastic shear strain");
        return Ok(());
    }
    let names = ["xx", "yy", "zz", "yz", "xz", "xy"];
    let c = names.iter().position(|n| *n == component).ok_or_else(|| format!("Unknown strain component {component}"))?;
    let mut values = Vec::new();
    for row in result.at("strain_voigt")?.items()? {
        let item = &row.items()?[c];
        values.push(if item.is_number() { item.number()? } else { f64::NAN });
    }
    set_property(data, values, &format!("Elastic strain {component}"));
    Ok(())
}

/// Displacements from a reference structure file (same atoms, same order).
fn displacement_step(p: &Json, data: &mut PipelineData) -> R {
    let id = "displacements";
    let file = text(p, id, "reference")?;
    if file.is_empty() {
        return Err("Choose the reference structure file".into());
    }
    let reference = af_science::data::read_frames(Path::new(&file), None)?.swap_remove(0).structure;
    let s = &data.structure;
    if reference.atoms.len() != s.atoms.len() {
        return Err(format!("The reference has {} atoms, the structure {}", reference.atoms.len(), s.atoms.len()));
    }
    let mut request = Json::object();
    let current = atoms_request(s).at("positions")?.clone();
    let reference_request = atoms_request(&reference);
    request.set("reference", reference_request.at("positions")?.clone());
    request.set("current", current);
    let mapping = text(p, id, "affine_mapping")?;
    if reference.has_unit_cell {
        request.set("reference_cell", reference_request.at("cell")?.clone());
        request.set("pbc", reference_request.at("pbc")?.clone());
        if s.has_unit_cell {
            request.set("current_cell", atoms_request(s).at("cell")?.clone());
        }
    }
    if mapping != "off" {
        request.set("affine_mapping", mapping);
    }
    let result = tool("displacements", &request)?;
    let component = text(p, id, "component")?;
    if component == "magnitude" {
        set_property(data, numbers_of(result.at("magnitude_A")?)?, "Displacement magnitude (A)");
        return Ok(());
    }
    let c = match component.as_str() {
        "x" => 0,
        "y" => 1,
        _ => 2,
    };
    let mut values = Vec::new();
    for row in result.at("displacement_A")?.items()? {
        values.push(row.items()?[c].number()?);
    }
    set_property(data, values, &format!("Displacement {component} (A)"));
    Ok(())
}

/// Statistics of the per-atom property (OVITO's Reduce property), as notes.
fn reduce_property(p: &Json, data: &mut PipelineData) -> R {
    let s = &data.structure;
    if s.atom_property.len() != s.atoms.len() || s.atoms.is_empty() {
        return Err("reduce-property needs a per-atom property; add compute-property above".into());
    }
    let selected_only = text(p, "reduce-property", "target")? == "selected";
    let values: Vec<f64> =
        (0..s.atoms.len()).filter(|&i| (!selected_only || data.selected[i] != 0) && s.atom_property[i].is_finite()).map(|i| s.atom_property[i]).collect();
    if values.is_empty() {
        return Err("No atom has a value".into());
    }
    let mut request = Json::object();
    request.set("values", Json::Array(values.iter().map(|&v| Json::from(v)).collect()));
    request.set("bins", number(p, "reduce-property", "bins")?);
    let r = tool("property-histogram", &request)?;
    let g = |key: &str| -> Result<String, String> { Ok(format_g(r.at(key)?.number()?, 6)) };
    let name = if s.atom_property_name.is_empty() { "property".to_string() } else { s.atom_property_name.clone() };
    let line = format!(
        "{name}: n={:.0} sum={} mean={} std={} min={} median={} max={}",
        r.at("count")?.number()?,
        g("sum")?,
        g("mean")?,
        g("std")?,
        g("min")?,
        g("median")?,
        g("max")?
    );
    let line = truncated(line, 256);
    let mut histogram = "histogram:".to_string();
    let (centers, counts) = (numbers_of(r.at("centers")?)?, numbers_of(r.at("counts")?)?);
    for b in 0..centers.len() {
        histogram += &truncated(format!(" {}:{:.0}", format_g(centers[b], 4), counts[b]), 256);
    }
    data.notes.push(line);
    data.notes.push(histogram);
    Ok(())
}

/// Spatial binning of the per-atom property (or atom counts); each atom takes
/// the value of its bin, so the profile or map can colour the structure.
fn bin_property(p: &Json, data: &mut PipelineData) -> R {
    let id = "bin-property";
    require_cell(data, id)?;
    let s = &data.structure;
    let mut request = atoms_request(s);
    let has_property = s.atom_property.len() == s.atoms.len() && !s.atoms.is_empty();
    if has_property {
        request.set("values", Json::Array(s.atom_property.iter().map(|&v| Json::from(v)).collect()));
    }
    let bins = vector(p, id, "bins")?;
    request.set("bins", Json::Array(bins.iter().map(|&v| Json::from(v)).collect()));
    request.set("operation", text(p, id, "operation")?);
    let result = tool("spatial-binning", &request)?;
    let grid = numbers_of(result.at("values")?)?;
    let cell = cell_of(s);
    let n: [i64; 3] = std::array::from_fn(|d| bins[d].round() as i64);
    let mut per_atom = vec![0.0; s.atoms.len()];
    for (i, atom) in s.atoms.iter().enumerate() {
        let mut f = fractional(position_of(atom), &cell)?;
        let mut index = [0i64; 3];
        for d in 0..3 {
            f[d] -= f[d].floor();
            index[d] = (n[d] - 1).min((f[d] * n[d] as f64) as i64);
        }
        per_atom[i] = grid[((index[0] * n[1] + index[1]) * n[2] + index[2]) as usize];
    }
    let operation = text(p, id, "operation")?;
    let name = format!("Binned {operation}{}", if has_property { format!(" of {}", s.atom_property_name) } else { " (atoms)".to_string() });
    set_property(data, per_atom, &name);
    Ok(())
}

/// Colours the selected (or all) atoms; the colour stays with the atoms.
fn assign_color(p: &Json, data: &mut PipelineData) -> R {
    let colour = vector(p, "assign-color", "color")?;
    for c in colour {
        if !(0.0..=1.0).contains(&c) {
            return Err("Colour components are between 0 and 1 (red green blue)".into());
        }
    }
    let indices = targets(data, text(p, "assign-color", "target")? == "selected", "assign-color")?;
    let s = &mut data.structure;
    let rgb = [colour[0] as f32, colour[1] as f32, colour[2] as f32];
    for &i in &indices {
        [s.atoms[i].r, s.atoms[i].g, s.atoms[i].b] = rgb;
    }
    if s.grain_colors.len() == s.atoms.len() {
        for &i in &indices {
            s.grain_colors[i] = rgb;
        }
    }
    data.notes.push(format!("{} atoms coloured", indices.len()));
    Ok(())
}

/// Darkens atom colours by ambient occlusion (covalent radii).
fn ambient_occlusion_step(p: &Json, data: &mut PipelineData) -> R {
    let id = "ambient-occlusion";
    if data.structure.atoms.is_empty() {
        return Ok(());
    }
    let covalent = af_core::elements::literature_covalent_radii();
    let s = &data.structure;
    let positions: Vec<Json> = s.atoms.iter().map(|a| Json::Array(vec![a.x.into(), a.y.into(), a.z.into()])).collect();
    let radii: Vec<Json> = s
        .atoms
        .iter()
        .map(|a| {
            let z = a.atomic_number.max(0) as usize;
            Json::from(if z < covalent.len() && covalent[z] > 0.0 { covalent[z] as f64 } else { 1.0 })
        })
        .collect();
    let mut request = Json::object();
    request.set("positions", Json::Array(positions));
    request.set("radii", Json::Array(radii));
    request.set("directions", number(p, id, "directions")?);
    let result = tool("ambient-occlusion", &request)?;
    let intensity = number(p, id, "intensity")?;
    if !(0.0..=1.0).contains(&intensity) {
        return Err("intensity is between 0 and 1".into());
    }
    let brightness = numbers_of(result.at("brightness")?)?;
    let s = &mut data.structure;
    let grains = s.grain_colors.len() == s.atoms.len();
    for i in 0..s.atoms.len() {
        let factor = (1.0 - intensity + intensity * brightness[i]) as f32;
        s.atoms[i].r *= factor;
        s.atoms[i].g *= factor;
        s.atoms[i].b *= factor;
        if grains {
            for c in s.grain_colors[i].iter_mut() {
                *c *= factor;
            }
        }
    }
    s.atom_property = brightness;
    s.atom_property_name = "Ambient brightness".into();
    data.notes.push(format!("Ambient occlusion from {} directions", number(p, id, "directions")? as i64));
    Ok(())
}

/// Keeps a copy of the per-atom property for later comparison.
fn freeze_property(_: &Json, data: &mut PipelineData) -> R {
    let s = &mut data.structure;
    if s.atom_property.len() != s.atoms.len() || s.atoms.is_empty() {
        return Err("freeze-property needs a per-atom property; add compute-property above".into());
    }
    s.frozen_property = s.atom_property.clone();
    s.frozen_property_name = s.atom_property_name.clone();
    let name = if s.atom_property_name.is_empty() { "the property".to_string() } else { s.atom_property_name.clone() };
    data.notes.push(format!("Kept {name} as frozen"));
    Ok(())
}

fn clear_property(_: &Json, data: &mut PipelineData) -> R {
    data.structure.atom_property.clear();
    data.structure.atom_property_name.clear();
    Ok(())
}

fn compute_expression(p: &Json, data: &mut PipelineData) -> R {
    let id = "compute-expression";
    let formula = text(p, id, "expression")?;
    let values = evaluate_for_atoms(&Expression::new(&formula)?, data)?;
    let selected_only = text(p, id, "target")? == "selected";
    let n = data.structure.atoms.len();
    let mut property = if data.structure.atom_property.len() == n { data.structure.atom_property.clone() } else { vec![f64::NAN; n] };
    for (i, v) in values.iter().enumerate() {
        if !selected_only || data.selected[i] != 0 {
            property[i] = *v;
        }
    }
    data.structure.atom_property = property;
    let name = text(p, id, "name")?;
    data.structure.atom_property_name = if name.is_empty() { formula } else { name };
    let count = if selected_only { data.selected_count() } else { values.len() };
    data.notes.push(format!("{} computed for {count} atoms", data.structure.atom_property_name));
    Ok(())
}

/// Grains from local lattice orientations; atoms are coloured by grain and the
/// grain id becomes the per-atom property.
fn segment_grains(p: &Json, data: &mut PipelineData) -> R {
    let id = "segment-grains";
    let mut request = atoms_request(&data.structure);
    request.set("merge_threshold_deg", number(p, id, "threshold")?);
    request.set("min_grain_size", number(p, id, "min_size")?);
    request.set("adopt_orphans", flag(p, id, "adopt")?);
    let result = tool("grain-segmentation", &request)?;
    let s = &mut data.structure;
    s.atom_property.clear();
    s.grain_colors.clear();
    for id in result.at("grain_id")?.items()? {
        let g = id.number()?;
        s.atom_property.push(g);
        // Golden-ratio hues give neighbouring grain ids distinct colours; grey for none.
        if g <= 0.0 {
            s.grain_colors.push([0.6, 0.6, 0.6]);
            continue;
        }
        let hue = (g * 0.6180339887) % 1.0 * 6.0;
        let sector = hue as i32;
        let f = (hue - sector as f64) as f32;
        let (v, lo) = (0.9f32, 0.25f32);
        let rise = lo + (v - lo) * f;
        let fall = v - (v - lo) * f;
        let colours = [[v, rise, lo], [fall, v, lo], [lo, v, rise], [lo, fall, v], [rise, lo, v], [v, lo, fall]];
        s.grain_colors.push(colours[(sector % 6) as usize]);
    }
    s.atom_property_name = "Grain id (0 none)".into();
    data.notes.push(format!("{} grains", result.at("grain_count")?.number()? as i64));
    Ok(())
}

/// The analyses of compute-property: the choice, the scientific tool, the
/// property it reports (empty: its first one) and, for analyses of one
/// species (Chill+ looks at water oxygens), that element.
const PROPERTY_CHOICES: [(&str, &str, &str, &str); 18] = [
    ("structure-type", "structure-type", "", ""),
    ("cna", "cna", "", ""),
    ("ptm", "ptm", "", ""),
    ("orientation", "ptm", property_name::ORIENTATION_ANGLE, ""),
    ("ptm-ordering", "ptm", property_name::PTM_ORDERING, ""),
    ("vorotop", "vorotop", "", ""),
    ("diamond", "diamond-structure", "", ""),
    ("centrosymmetry", "centrosymmetry", "", ""),
    ("coordination", "bond-order", property_name::COORDINATION, ""),
    ("q6", "bond-order", "", ""),
    ("bond-valence", "bond-orders", "", ""),
    ("voronoi-volume", "voronoi", "", ""),
    ("voronoi-coordination", "voronoi", property_name::VORONOI_COORDINATION, ""),
    ("surface", "surface-mesh", "", ""),
    ("planar-faults", "fcc-planar-faults", "", ""),
    ("chill-plus", "chill-plus", "", "O"),
    ("local-entropy", "local-entropy", "", ""),
    ("rings", "find-rings", "", ""),
];

fn property_choice_options() -> String {
    PROPERTY_CHOICES.iter().map(|c| c.0).collect::<Vec<_>>().join("|")
}

fn analyse(choice: &(&str, &str, &str, &str), structure: &Structure, cutoff: f64) -> Result<AtomProperty, String> {
    let (_, tool_name, property, element) = *choice;
    // Only the atoms of the element, when the structure has others too; the rest get no value.
    let members: Vec<usize> = if element.is_empty() { Vec::new() } else { (0..structure.atoms.len()).filter(|&i| structure.atoms[i].symbol == element).collect() };
    if members.is_empty() || members.len() == structure.atoms.len() {
        return analyse_per_atom(structure, tool_name, cutoff, property);
    }
    let mut subset = structure.clone();
    subset.atoms = members.iter().map(|&i| structure.atoms[i].clone()).collect();
    let on_subset = analyse_per_atom(&subset, tool_name, cutoff, property)?;
    let mut values = vec![f64::NAN; structure.atoms.len()];
    for (k, &i) in members.iter().enumerate() {
        values[i] = on_subset.values[k];
    }
    Ok(AtomProperty { name: on_subset.name, values })
}

fn compute_property(p: &Json, data: &mut PipelineData) -> R {
    let analysis = text(p, "compute-property", "analysis")?;
    let Some(choice) = PROPERTY_CHOICES.iter().find(|c| c.0 == analysis) else {
        return Err(format!("Unknown per-atom property {analysis}; choose from {}", property_choice_options()));
    };
    let property = analyse(choice, &data.structure, number(p, "compute-property", "cutoff")?)?;
    data.notes.push(format!("{} computed for {} atoms", property.name, property.values.len()));
    data.structure.atom_property = property.values;
    data.structure.atom_property_name = property.name;
    Ok(())
}

pub fn types() -> Vec<ModifierType> {
    let with_mode = |mut list: Vec<ModifierParameter>| {
        list.push(mode_parameter());
        list
    };
    let options = property_choice_options();
    vec![
        ModifierType::new(
            "select-overlapping",
            "Select overlapping atoms",
            "Selection",
            "Selects atoms closer than a distance to another atom (minimum image in periodic cells): the later atom of each pair, or both.",
            with_mode(params(&[("distance", "Distance (Angstrom)", "float", "1", ""), ("atoms", "Atoms to select", "choice", "one", "one|both")])),
            select_overlapping,
        ),
        ModifierType::new(
            "select-surface",
            "Select surface atoms",
            "Selection",
            "Selects the atoms on the surface of the solid: those a probe sphere of the given radius touches (as in Analysis > Surface Mesh); reports the surface area and solid volume.",
            with_mode(params(&[("probe_radius", "Probe radius (Angstrom)", "float", "4", "")])),
            select_surface,
        ),
        ModifierType::new(
            "select-region",
            "Select a LAMMPS region",
            "Selection",
            "Selects the atoms inside LAMMPS regions given by their region commands, separated by ';' (e.g. region top block INF INF INF INF 10 INF units box);\nunits lattice multiplies coordinates by the lattice spacing.",
            with_mode(params(&[
                ("regions", "Region commands", "string", "region r1 sphere 0 0 0 5 units box", ""),
                ("lattice_scale", "Lattice spacing (units lattice)", "float", "1", ""),
            ])),
            select_region,
        ),
        ModifierType::new(
            "assign-color",
            "Assign colour",
            "Modification",
            "Colours the selected (or all) atoms with a red green blue colour (components 0 to 1).",
            params(&[("color", "Colour (red green blue, 0-1)", "vector", "1 0 0", ""), ("target", "Atoms", "choice", "selected", "selected|all")]),
            assign_color,
        ),
        ModifierType::new(
            "ambient-occlusion",
            "Ambient occlusion",
            "Modification",
            "Darkens the atoms by how hidden they are from all directions (covalent radii), which brings out pores and grooves; the brightness becomes the per-atom property.",
            params(&[("intensity", "Intensity (0-1)", "float", "0.7", ""), ("directions", "Directions", "int", "40", "")]),
            ambient_occlusion_step,
        ),
        ModifierType::new("clear-property", "Remove per-atom property", "Analysis", "Removes the per-atom property (and its colouring).", Vec::new(), clear_property),
        ModifierType::new(
            "freeze-property",
            "Freeze per-atom property",
            "Analysis",
            "Keeps a copy of the per-atom property; later expressions compare with it as the variable frozen, e.g. compute-property cna | freeze-property |\nrelax ... | compute-property cna | select-expression \"property != frozen\" selects the atoms whose structure changed.",
            Vec::new(),
            freeze_property,
        ),
        ModifierType::new(
            "compute-property",
            "Compute per-atom property",
            "Analysis",
            "Computes a per-atom property for colouring and select-property: structure types (Ackland-Jones, adaptive CNA, template matching with its\norientation angle and chemical ordering, VoroTop, diamond), centrosymmetry, coordination, Steinhardt q6, bond valence, Voronoi volume or\ncoordination, surface atoms, fcc planar faults, Chill+ ice types (on the oxygens), local entropy or rings through each atom. The cutoff is\nused by analyses that need one; the others choose their own neighbour range.",
            params(&[("analysis", "Property", "choice", "structure-type", options.as_str()), ("cutoff", "Neighbour cutoff (Angstrom)", "float", "3", "")]),
            compute_property,
        ),
        ModifierType::new(
            "compute-expression",
            "Compute property from a formula",
            "Analysis",
            "Sets the per-atom property to a formula of the atom variables, e.g. sqrt(x^2 + y^2) or property * 2 (see select-expression for the variables);\nwith target selected, only selected atoms change.",
            params(&[
                ("expression", "Formula", "expression", "z", ""),
                ("name", "Property name (empty: the formula)", "string", "", ""),
                ("target", "Atoms", "choice", "all", "all|selected"),
            ]),
            compute_expression,
        ),
        ModifierType::new(
            "elastic-strain",
            "Elastic strain",
            "Analysis",
            "Sets the per-atom property to the elastic strain relative to the ideal lattice (polyhedral template matching): volumetric (trace/3), von Mises\nshear or one Green-Lagrange component. Atoms of other structures get no value.",
            params(&[
                ("lattice", "Lattice", "choice", "fcc", "fcc|hcp|bcc|sc|cubic-diamond|hex-diamond"),
                ("a", "Unstrained lattice constant (Angstrom)", "float", "3.615", ""),
                ("c_over_a", "Unstrained c/a (hexagonal)", "float", "1.633", ""),
                ("component", "Property", "choice", "volumetric", "volumetric|shear|xx|yy|zz|yz|xz|xy"),
            ]),
            elastic_strain_step,
        ),
        ModifierType::new(
            "displacements",
            "Displacement vectors",
            "Analysis",
            "Sets the per-atom property to the displacement from a reference structure file (same atoms in the same order; minimum image in periodic\ncells; optionally after mapping the cell deformation).",
            params(&[
                ("reference", "Reference structure file", "string", "", ""),
                ("component", "Property", "choice", "magnitude", "magnitude|x|y|z"),
                ("affine_mapping", "Affine mapping", "choice", "off", "off|to_reference|to_current"),
            ]),
            displacement_step,
        ),
        ModifierType::new(
            "bin-property",
            "Spatial binning",
            "Analysis",
            "Averages (or sums, counts, ...) the per-atom property over bins along a, b and c (1 = not binned) and gives each atom its bin's value;\nwithout a property the atoms are counted. bins 1 1 20 gives a profile along c.",
            params(&[("bins", "Bins along a, b, c", "ivec3", "1 1 20", ""), ("operation", "Operation", "choice", "mean", "mean|sum|density|min|max|count")]),
            bin_property,
        ),
        ModifierType::new(
            "reduce-property",
            "Property statistics",
            "Analysis",
            "Reports the count, sum, mean, standard deviation, minimum, median and maximum and a histogram of the per-atom property (all or selected atoms).",
            params(&[("target", "Atoms", "choice", "all", "all|selected"), ("bins", "Histogram bins", "int", "10", "")]),
            reduce_property,
        ),
        ModifierType::new(
            "segment-grains",
            "Segment grains",
            "Analysis",
            "Finds the grains of a polycrystal from local lattice orientations (polyhedral template matching), colours atoms by grain and sets the\nper-atom property to the grain id (0 for none).",
            params(&[
                ("threshold", "Misorientation threshold (degrees)", "float", "4", ""),
                ("min_size", "Minimum grain size (atoms)", "int", "100", ""),
                ("adopt", "Assign boundary atoms to grains", "bool", "true", ""),
            ]),
            segment_grains,
        ),
    ]
}
