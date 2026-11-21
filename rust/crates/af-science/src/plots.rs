//! Plots of scientific tool results (the C++ science/ResultPlots): x/y
//! series, labelled vertical markers and horizontal reference lines for each
//! tool with a curve to show, run overlays, histograms and CSV export.

use af_core::cxx::sort_by;
use af_core::json::Json;

type R<T> = Result<T, String>;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PlotSeries {
    pub name: String,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Markers instead of a line.
    pub points: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PlotMarker {
    pub x: f64,
    pub label: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PlotSpec {
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub series: Vec<PlotSeries>,
    /// Labelled vertical lines (e.g. k-path points).
    pub markers: Vec<PlotMarker>,
    /// Reference horizontal lines.
    pub horizontal: Vec<f64>,
    pub log_y: bool,
}

const NAN: f64 = f64::NAN;

fn numbers(value: Option<&Json>) -> Vec<f64> {
    match value {
        Some(Json::Array(items)) => items
            .iter()
            .map(|item| match item {
                Json::Number(n) => *n,
                Json::Bool(b) => f64::from(u8::from(*b)),
                _ => NAN,
            })
            .collect(),
        _ => Vec::new(),
    }
}

fn column(rows: Option<&Json>, k: usize) -> Vec<f64> {
    match rows {
        Some(Json::Array(items)) => items
            .iter()
            .map(|row| match row {
                Json::Array(r) if r.len() > k => match r[k] {
                    Json::Number(n) => n,
                    _ => NAN,
                },
                _ => NAN,
            })
            .collect(),
        _ => Vec::new(),
    }
}

fn indices(count: usize) -> Vec<f64> {
    (0..count).map(|i| i as f64).collect()
}

fn line(name: &str, x: Vec<f64>, y: Vec<f64>, points: bool) -> PlotSeries {
    PlotSeries { name: name.to_string(), x, y, points }
}

fn plot(title: &str, x_label: &str, y_label: &str, series: Vec<PlotSeries>) -> PlotSpec {
    PlotSpec { title: title.to_string(), x_label: x_label.to_string(), y_label: y_label.to_string(), series, ..Default::default() }
}

fn number(result: &Json, key: &str) -> f64 {
    match result.find(key) {
        Some(Json::Number(n)) => *n,
        _ => NAN,
    }
}

fn items(j: &Json) -> R<&Vec<Json>> {
    j.items().map_err(|e| e.to_string())
}

fn members(j: &Json) -> R<&Vec<(String, Json)>> {
    j.members().map_err(|e| e.to_string())
}

fn num(j: &Json) -> R<f64> {
    j.number().map_err(|e| e.to_string())
}

fn at<'a>(j: &'a Json, key: &str) -> R<&'a Json> {
    j.at(key).map_err(|e| e.to_string())
}

fn item(j: &Json, i: usize) -> R<&Json> {
    items(j)?.get(i).ok_or_else(|| "index out of range".to_string())
}

fn add_histogram(plots: &mut Vec<PlotSpec>, title: &str, x_label: &str, values: &[f64]) {
    if values.iter().any(|v| v.is_finite()) {
        plots.push(plot(title, x_label, "Atoms", vec![histogram(x_label, values, 40)]));
    }
}

/// Bar-like plot of the "counts" object of a classification (index order).
fn count_plot(result: &Json, plots: &mut Vec<PlotSpec>, title: &str) -> R<()> {
    if let Some(counts @ Json::Object(_)) = result.find("counts") {
        let (mut x, mut y) = (Vec::new(), Vec::new());
        for (t, (_, v)) in members(counts)?.iter().enumerate() {
            x.push(t as f64);
            y.push(num(v)?);
        }
        plots.push(plot(title, "Type", "Atoms", vec![line("atoms", x, y, true)]));
    }
    Ok(())
}

fn is_object(j: Option<&Json>) -> Option<&Json> {
    j.filter(|v| v.is_object())
}

fn is_array(j: Option<&Json>) -> Option<&Json> {
    j.filter(|v| v.is_array())
}

fn msd(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let lag = numbers(r.find("lag_fs"));
    let c = r.find("components_A2");
    plots.push(plot(
        "Mean-square displacement",
        "Lag time (fs)",
        "MSD (A^2)",
        vec![line("total", lag.clone(), numbers(r.find("msd_A2")), false), line("x", lag.clone(), column(c, 0), false), line("y", lag.clone(), column(c, 1), false), line("z", lag, column(c, 2), false)],
    ));
    Ok(())
}

fn diffusion(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let lag = numbers(r.find("lag_fs"));
    let range = numbers(r.find("fit_range_fs"));
    if !lag.is_empty() && range.len() == 2 {
        let (slope, intercept) = (number(r, "slope_A2_per_fs"), number(r, "intercept_A2"));
        let fit = vec![intercept + slope * range[0], intercept + slope * range[1]];
        plots.push(plot("Einstein diffusion fit", "Lag time (fs)", "MSD (A^2)", vec![line("MSD", lag, numbers(r.find("msd_A2")), false), line("linear fit", range, fit, false)]));
    }
    Ok(())
}

fn vacf(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let mut spec = plot("Velocity autocorrelation", "Lag time (fs)", "VACF", vec![line("VACF", numbers(r.find("lag_fs")), numbers(r.find("vacf")), false)]);
    spec.horizontal = vec![0.0];
    plots.push(spec);
    Ok(())
}

fn vibrational_spectrum(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    plots.push(plot("Vibrational spectrum", "Frequency (THz)", "Spectral density (1/THz)", vec![line("density", numbers(r.find("frequency_THz")), numbers(r.find("density_per_THz")), false)]));
    Ok(())
}

fn structure_factor(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let q = r.find("q_vectors_rad_per_A");
    let s = numbers(r.find("S_q"));
    let (qx, qy, qz) = (column(q, 0), column(q, 1), column(q, 2));
    let get = |v: &Vec<f64>, i: usize| v.get(i).copied().unwrap_or(NAN);
    let mut pairs: Vec<(f64, f64)> = (0..s.len()).map(|i| ((get(&qx, i) * get(&qx, i) + get(&qy, i) * get(&qy, i) + get(&qz, i) * get(&qz, i)).sqrt(), s[i])).collect();
    // std::pair's operator<.
    sort_by(&mut pairs, |a, b| a.0 < b.0 || (!(b.0 < a.0) && a.1 < b.1));
    let (x, y) = pairs.into_iter().unzip();
    plots.push(plot("Static structure factor", "|q| (rad/A)", "S(q)", vec![line("S(q)", x, y, true)]));
    Ok(())
}

fn local_strain(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    add_histogram(plots, "Non-affine displacement", "D2min (A^2)", &numbers(r.find("d2min_A2")));
    let mut volumetric = Vec::new();
    if let Some(strain) = is_array(r.find("green_lagrange_strain")) {
        for tensor in items(strain)? {
            let mut trace = 0.0;
            for k in 0..3 {
                let row = match tensor {
                    Json::Array(t) if t.len() == 3 => t[k].clone(),
                    _ => Json::Null,
                };
                trace += match &row {
                    Json::Array(v) => match v.get(k) {
                        Some(Json::Number(n)) => *n,
                        Some(_) => NAN,
                        None => return Err("index out of range".into()),
                    },
                    _ => NAN,
                };
            }
            volumetric.push(trace);
        }
    }
    add_histogram(plots, "Volumetric strain", "trace(E)", &volumetric);
    Ok(())
}

fn vasp_electronic(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let distance = numbers(r.find("distance"));
    if let Some(energies) = is_array(r.find("energies_minus_fermi_eV")) {
        if distance.len() > 1 {
            let spins = items(energies)?;
            let mut series = Vec::new();
            for (s, spin) in spins.iter().enumerate() {
                let count = items(item(spin, 0)?)?.len();
                for b in 0..count {
                    let prefix = if spins.len() > 1 { if s == 0 { "up " } else { "down " } } else { "band " };
                    series.push(line(&format!("{prefix}{}", b + 1), distance.clone(), column(Some(spin), b), false));
                }
            }
            let mut spec = plot("Band structure", "Wavevector distance", "E - E_F (eV)", series);
            if let Some(labels) = is_array(r.find("labels")) {
                for marker in items(labels)? {
                    spec.markers.push(PlotMarker { x: num(at(marker, "distance")?)?, label: at(marker, "label")?.string().map_err(|e| e.to_string())?.to_string() });
                }
            }
            spec.horizontal = vec![0.0];
            plots.push(spec.clone());
            if let Some(fat) = is_object(r.find("fat_band_weights")) {
                if fat.size() > 0 {
                    let (element, weights) = &members(fat)?[0];
                    let (mut x, mut y) = (Vec::new(), Vec::new());
                    let spin = &spins[0];
                    for (k, &d) in distance.iter().enumerate() {
                        let bands = items(item(spin, k)?)?;
                        for (b, energy) in bands.iter().enumerate() {
                            if num(item(item(item(weights, 0)?, k)?, b)?)? > 0.3 {
                                x.push(d);
                                y.push(num(energy)?);
                            }
                        }
                    }
                    let mut fat_spec = spec;
                    fat_spec.title = format!("Fat bands: states with >30% {element} character");
                    fat_spec.series.push(line(&format!("{element} > 30%"), x, y, true));
                    plots.push(fat_spec);
                }
            }
        }
    }
    if let Some(dos) = is_object(r.find("dos")) {
        let energy = numbers(dos.find("energy_minus_fermi_eV"));
        let mut series = Vec::new();
        if let Some(total) = dos.find("total") {
            for s in 0..total.size() {
                series.push(line(if s > 0 { "total (down)" } else { "total" }, energy.clone(), numbers(Some(item(total, s)?)), false));
            }
        }
        if let Some(projected) = is_object(dos.find("projected_by_element")) {
            for (element, data) in members(projected)? {
                series.push(line(element, energy.clone(), numbers(Some(item(data, 0)?)), false));
            }
        }
        let mut spec = plot("Density of states", "E - E_F (eV)", "DOS (states/eV)", series);
        spec.markers.push(PlotMarker { x: 0.0, label: "E_F".into() });
        plots.push(spec);
    }
    Ok(())
}

fn powder_xrd(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let mut spec = plot("Powder diffraction pattern", "2-theta (degrees)", "Intensity", vec![line("profile", numbers(r.find("two_theta_deg")), numbers(r.find("profile")), false)]);
    if let Some(peaks) = is_array(r.find("peaks")) {
        for peak in items(peaks)? {
            if num(at(peak, "intensity")?)? >= 5.0 {
                // The conventional cell's indices when known; digits run together
                // only when each is a single non-negative digit ("111", "1 10 0").
                let h = items(match peak.find("hkl_conventional") {
                    Some(c) => c,
                    None => at(peak, "hkl")?,
                })?;
                let values: Vec<f64> = h.iter().map(num).collect::<R<_>>()?;
                let compact = values.iter().all(|&v| (0.0..=9.0).contains(&v));
                let mut label = String::new();
                for v in values {
                    if !label.is_empty() && !compact {
                        label.push(' ');
                    }
                    label += &(v as i32).to_string();
                }
                spec.markers.push(PlotMarker { x: num(at(peak, "two_theta_deg")?)?, label });
            }
        }
    }
    plots.push(spec);
    Ok(())
}

fn electron_diffraction(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    plots.push(plot("Diffraction spots", "x (1/A)", "y (1/A)", vec![line("spots", numbers(r.find("spot_x_inv_A")), numbers(r.find("spot_y_inv_A")), true)]));
    Ok(())
}

/// Rows [value, count] as x and y; a single row gets a zero beside it.
fn row_pairs(rows: &Json) -> R<(Vec<f64>, Vec<f64>)> {
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for row in items(rows)? {
        x.push(num(item(row, 0)?)?);
        y.push(num(item(row, 1)?)?);
    }
    if x.len() == 1 {
        x.push(x[0] + 1.0);
        y.push(0.0);
    }
    Ok((x, y))
}

fn cluster_analysis(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    if let Some(h) = is_array(r.find("size_histogram")) {
        if h.size() > 0 {
            let (x, y) = row_pairs(h)?;
            plots.push(plot("Cluster size distribution", "Cluster size (atoms)", "Clusters", vec![line("clusters", x, y, true)]));
        }
    }
    Ok(())
}

fn void_analysis(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let (mut rank, mut volume) = (Vec::new(), Vec::new());
    if let Some(voids) = is_array(r.find("voids")) {
        for v in items(voids)? {
            rank.push((rank.len() + 1) as f64);
            volume.push(num(at(v, "accessible_volume_A3")?)?);
        }
    }
    if rank.len() >= 2 {
        plots.push(plot("Void volumes", "Void (largest first)", "Accessible volume (A^3)", vec![line("voids", rank, volume, true)]));
    }
    Ok(())
}

fn bond_order(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let mut series = Vec::new();
    if let Some(order) = is_object(r.find("order")) {
        for (name, values) in members(order)? {
            let h = histogram(name, &numbers(Some(values)), 40);
            if !h.x.is_empty() {
                series.push(h);
            }
        }
    }
    if !series.is_empty() {
        plots.push(plot("Steinhardt order distribution", "q_l", "Atoms", series));
    }
    Ok(())
}

fn work_function(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let distance = numbers(r.find("distance_A"));
    if !distance.is_empty() {
        let mut spec = plot("Planar potential", "Distance (A)", "Potential energy (eV)", vec![line("potential", distance, numbers(r.find("potential_eV")), false)]);
        spec.horizontal = vec![number(r, "vacuum_level_eV"), number(r, "fermi_eV")];
        plots.push(spec);
    }
    Ok(())
}

fn equation_of_state(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    plots.push(plot(
        "Birch-Murnaghan equation of state",
        "Volume (A^3)",
        "Energy (eV)",
        vec![line("data", numbers(r.find("volumes_A3")), numbers(r.find("energies_eV")), true), line("fit", numbers(r.find("fit_volumes_A3")), numbers(r.find("fit_energies_eV")), false)],
    ));
    Ok(())
}

fn phonon_dos(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    plots.push(plot("Phonon density of states", "Energy (eV)", "DOS (1/eV)", vec![line("DOS", numbers(r.find("energy_eV")), numbers(r.find("dos_per_eV")), false)]));
    Ok(())
}

fn thermodynamics(t: Vec<f64>, src: &Json, plots: &mut Vec<PlotSpec>) {
    plots.push(plot(
        "Harmonic free and internal energy",
        "Temperature (K)",
        "Energy (eV/cell)",
        vec![line("F", t.clone(), numbers(src.find("free_energy_eV")), false), line("U", t, numbers(src.find("internal_energy_eV")), false)],
    ));
}

fn harmonic_thermodynamics(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let t = numbers(r.find("temperature_K"));
    thermodynamics(t.clone(), r, plots);
    plots.push(plot(
        "Entropy and heat capacity",
        "Temperature (K)",
        "eV/(cell K)",
        vec![line("S", t.clone(), numbers(r.find("entropy_eV_per_K")), false), line("Cv", t, numbers(r.find("heat_capacity_eV_per_K")), false)],
    ));
    Ok(())
}

fn phonons(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let distance = numbers(r.find("distance_inv_A"));
    if let Some(bands) = is_array(r.find("frequencies_THz")) {
        if bands.size() == distance.len() && bands.size() > 1 {
            let count = item(bands, 0)?.size();
            let branches = (0..count).map(|b| line(&format!("branch {}", b + 1), distance.clone(), column(Some(bands), b), false)).collect();
            let mut spec = plot("Phonon dispersion", "Wavevector distance (1/A)", "Frequency (THz)", branches);
            if let Some(labels) = is_array(r.find("labels")) {
                for (i, l) in items(labels)?.iter().enumerate().take(distance.len()) {
                    if let Json::String(text) = l {
                        if !text.is_empty() {
                            spec.markers.push(PlotMarker { x: distance[i], label: if text == "GAMMA" { "G".into() } else { text.clone() } });
                        }
                    }
                }
            }
            spec.horizontal = vec![0.0];
            plots.push(spec);
        }
    }
    if let Some(dos) = is_object(r.find("dos")) {
        let energy: Vec<f64> = numbers(dos.find("energy_eV")).into_iter().map(|e| e / 4.135667696e-3).collect();
        let density: Vec<f64> = numbers(dos.find("dos_per_eV")).into_iter().map(|d| d * 4.135667696e-3).collect();
        plots.push(plot("Phonon density of states", "Frequency (THz)", "DOS (states/THz)", vec![line("DOS", energy, density, false)]));
    }
    if let Some(thermo) = is_object(r.find("thermodynamics")) {
        let t = numbers(thermo.find("temperature_K"));
        thermodynamics(t.clone(), thermo, plots);
        plots.push(plot("Heat capacity", "Temperature (K)", "Cv (eV/(cell K))", vec![line("Cv", t, numbers(thermo.find("heat_capacity_eV_per_K")), false)]));
    }
    Ok(())
}

fn relax(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let enthalpy = numbers(r.find("enthalpy_history_eV"));
    let force = numbers(r.find("max_force_history"));
    plots.push(plot("Relaxation energy", "Step", "E (+PV) (eV)", vec![line("energy", indices(enthalpy.len()), enthalpy, false)]));
    let mut spec = plot("Convergence", "Step", "Max force (eV/A)", vec![line("max force", indices(force.len()), force, false)]);
    spec.log_y = true;
    plots.push(spec);
    Ok(())
}

fn neb(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let mut energies = numbers(r.find("energies_eV"));
    if !energies.is_empty() {
        let start = energies[0];
        for e in energies.iter_mut() {
            *e -= start;
        }
        let mut x = numbers(r.find("reaction_coordinate_A"));
        if x.len() != energies.len() {
            x = indices(energies.len());
        }
        let mut spec = plot("Minimum-energy path", "Reaction coordinate (A)", "E - E(initial) (eV)", vec![line("path", x.clone(), energies.clone(), false), line("images", x, energies, true)]);
        spec.horizontal = vec![0.0];
        plots.push(spec);
    }
    Ok(())
}

fn dynamics(r: &Json, plots: &mut Vec<PlotSpec>, npt: bool) -> R<()> {
    let time = numbers(r.find("time_fs"));
    let mut temperature = plot("Temperature", "Time (fs)", "T (K)", vec![line("T", time.clone(), numbers(r.find("temperature_K")), false)]);
    temperature.horizontal = vec![number(r, "target_temperature_K")];
    let mut energy = plot("Total energy", "Time (fs)", "E (eV)", vec![line("E", time.clone(), numbers(r.find("total_energy_eV")), false)]);
    if let Some(Json::Number(start)) = r.find("production_start_fs") {
        temperature.markers.push(PlotMarker { x: *start, label: "NVE".into() });
        energy.markers.push(PlotMarker { x: *start, label: "NVE".into() });
    }
    plots.push(temperature);
    plots.push(energy);
    if npt {
        plots.push(plot("Volume", "Time (fs)", "V (A^3)", vec![line("V", time.clone(), numbers(r.find("volume_A3")), false)]));
        plots.push(plot("Pressure", "Time (fs)", "P (GPa)", vec![line("P", time, numbers(r.find("pressure_GPa")), false)]));
    }
    Ok(())
}

fn trajectory_structure(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let radius = numbers(r.find("r_A"));
    let mut rdf = vec![line("total", radius.clone(), numbers(r.find("g_total")), false)];
    if let Some(partial) = is_object(r.find("g_partial")) {
        if partial.size() > 1 {
            for (pair, values) in members(partial)? {
                rdf.push(line(pair, radius.clone(), numbers(Some(values)), false));
            }
        }
    }
    let mut spec = plot("Radial distribution function", "r (A)", "g(r)", rdf);
    spec.horizontal = vec![1.0];
    plots.push(spec);
    if let Some(distribution) = is_object(r.find("coordination_distribution")) {
        let mut series = Vec::new();
        for (element, rows) in members(distribution)? {
            let (x, y) = row_pairs(rows)?;
            series.push(line(element, x, y, true));
        }
        plots.push(plot("Coordination distribution", "Neighbours within the cutoff", "Fraction of atoms", series));
    }
    plots.push(plot("Bond-angle distribution", "Angle (degrees)", "Density (1/degree)", vec![line("angles", numbers(r.find("angle_deg")), numbers(r.find("angle_density_per_deg")), false)]));
    Ok(())
}

fn time_series(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let frame = numbers(r.find("frame"));
    if let Some(counts) = is_object(r.find("counts")) {
        if counts.size() > 0 {
            let series = members(counts)?.iter().map(|(name, values)| line(name, frame.clone(), numbers(Some(values)), false)).collect();
            plots.push(plot("Structure types versus frame", "Frame", "Atoms", series));
        }
    }
    if let Some(mean) = r.find("mean_property") {
        let title = match r.find("property") {
            Some(Json::String(name)) => format!("Mean {name}"),
            _ => "Mean property".into(),
        };
        plots.push(plot(&title, "Frame", "Mean", vec![line("mean", frame.clone(), numbers(Some(mean)), false)]));
    }
    plots.push(plot("Cell volume", "Frame", "Volume (A^3)", vec![line("volume", frame, numbers(r.find("volume_A3")), false)]));
    Ok(())
}

fn surface_mesh(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let pores = numbers(r.find("empty_region_volumes_A3"));
    if pores.len() > 1 {
        plots.push(plot("Empty regions (largest first)", "Region", "Volume (A^3)", vec![line("volume", indices(pores.len()), pores, true)]));
    }
    let solids = numbers(r.find("solid_region_volumes_A3"));
    if !solids.is_empty() {
        plots.push(plot("Solid regions (largest first)", "Region", "Volume (A^3)", vec![line("volume", indices(solids.len()), solids, true)]));
    }
    Ok(())
}

fn scatter(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let x = numbers(r.find("x"));
    let mut series = vec![line("atoms", x.clone(), numbers(r.find("y")), true)];
    let (slope, intercept) = (number(r, "slope"), number(r, "intercept"));
    if slope.is_finite() && !x.is_empty() {
        // std::minmax_element: the first smallest, the last largest.
        let mut low = x[0];
        let mut high = x[0];
        for &v in &x[1..] {
            if v < low {
                low = v;
            }
            if !(v < high) {
                high = v;
            }
        }
        series.push(line("linear fit", vec![low, high], vec![intercept + slope * low, intercept + slope * high], false));
    }
    plots.push(plot("Scatter plot", "x", "y", series));
    Ok(())
}

fn bond_order_assignment(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let mut orders = Vec::new();
    if let Some(bonds) = is_array(r.find("bonds")) {
        for bond in items(bonds)? {
            orders.push(num(at(bond, "order")?)?);
        }
    }
    if !orders.is_empty() {
        plots.push(plot("Bond orders", "Bond order", "Bonds", vec![histogram("Bond order", &orders, 40)]));
    }
    add_histogram(plots, "Bond-valence sums", "Sum of bond orders", &numbers(r.find("valence_sum")));
    Ok(())
}

fn voronoi(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    add_histogram(plots, "Atomic volume", "Volume (A^3)", &numbers(r.find("atomic_volume_A3")));
    add_histogram(plots, "Voronoi coordination", "Faces", &numbers(r.find("coordination")));
    if let Some(common) = is_array(r.find("most_common_indices")) {
        if common.size() > 0 {
            let (mut rank, mut atoms) = (Vec::new(), Vec::new());
            for c in items(common)? {
                rank.push((rank.len() + 1) as f64);
                atoms.push(num(at(c, "atoms")?)?);
            }
            plots.push(plot("Most common Voronoi indices (rank; see the report)", "Rank", "Atoms", vec![line("atoms", rank, atoms, true)]));
        }
    }
    Ok(())
}

fn ptm(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    count_plot(r, plots, "PTM structures (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico, 5 sc, 6 cubic dia., 7 hex dia.)")?;
    add_histogram(plots, "PTM RMSD", "RMSD", &numbers(r.find("rmsd")));
    add_histogram(plots, "Orientation angle", "Rotation angle (degrees)", &numbers(r.find("orientation_angle_deg")));
    Ok(())
}

fn grain_segmentation(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    if let Some(grains) = is_array(r.find("grains")) {
        if grains.size() > 0 {
            let (mut id, mut size) = (Vec::new(), Vec::new());
            for g in items(grains)? {
                id.push(num(at(g, "id")?)?);
                size.push(num(at(g, "atoms")?)?);
            }
            plots.push(plot("Grain sizes", "Grain", "Atoms", vec![line("atoms", id, size, true)]));
        }
    }
    if let Some(boundaries) = is_array(r.find("boundaries")) {
        let mut angles = Vec::new();
        for b in items(boundaries)? {
            if let Json::Number(n) = at(b, "misorientation_deg")? {
                angles.push(*n);
            }
        }
        add_histogram(plots, "Grain-boundary misorientation", "Misorientation (degrees)", &angles);
    }
    Ok(())
}

fn bond_distributions(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let lengths = numbers(r.find("length_A"));
    let angles = numbers(r.find("angle_deg"));
    let mut length_series = vec![line("all", lengths.clone(), numbers(r.find("bond_length_counts")), false)];
    if let Some(pairs) = is_object(r.find("bond_length_counts_by_pair")) {
        if pairs.size() > 1 {
            for (name, values) in members(pairs)? {
                length_series.push(line(name, lengths.clone(), numbers(Some(values)), false));
            }
        }
    }
    plots.push(plot("Bond lengths", "Length (A)", "Bonds", length_series));
    let mut angle_series = vec![line("all", angles.clone(), numbers(r.find("bond_angle_counts")), false)];
    if let Some(triplets) = is_object(r.find("bond_angle_counts_by_triplet")) {
        if triplets.size() > 1 && triplets.size() <= 8 {
            for (name, values) in members(triplets)? {
                angle_series.push(line(name, angles.clone(), numbers(Some(values)), false));
            }
        }
    }
    plots.push(plot("Bond angles", "Angle (degrees)", "Angles", angle_series));
    Ok(())
}

fn spatial_binning(r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    let bins = numbers(r.find("bins"));
    let Some(centers) = is_object(r.find("centers")) else { return Ok(()) };
    if bins.len() != 3 {
        return Ok(());
    }
    let binned: Vec<usize> = (0..3).filter(|&d| bins[d] > 1.0).collect();
    if binned.len() != 1 {
        return Ok(());
    }
    let name = ["a", "b", "c"][binned[0]];
    let Some(along) = centers.find(name) else { return Ok(()) };
    plots.push(plot(&format!("Profile along {name}"), "Distance (A)", "Value", vec![line("value", numbers(along.find("distance_A")), numbers(r.find("values")), false)]));
    Ok(())
}

/// The plots of a tool's result, or an error where the C++ would throw.
fn build(tool: &str, r: &Json, plots: &mut Vec<PlotSpec>) -> R<()> {
    match tool {
        "msd" => msd(r, plots),
        "diffusion" => diffusion(r, plots),
        "vacf" => vacf(r, plots),
        "vibrational-spectrum" => vibrational_spectrum(r, plots),
        "structure-factor" => structure_factor(r, plots),
        "local-strain" => local_strain(r, plots),
        "vasp-electronic" => vasp_electronic(r, plots),
        "powder-xrd" => powder_xrd(r, plots),
        "electron-diffraction" => electron_diffraction(r, plots),
        "cluster-analysis" => cluster_analysis(r, plots),
        "void-analysis" => void_analysis(r, plots),
        "structure-type" => count_plot(r, plots, "Structure types (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)"),
        "diamond-structure" => count_plot(r, plots, "Diamond types (0 other, 1 cubic, 2-3 its neighbours, 4 hexagonal, 5-6 its neighbours)"),
        "centrosymmetry" => {
            add_histogram(plots, "Centrosymmetry distribution", "CSP (A^2)", &numbers(r.find("centrosymmetry_A2")));
            Ok(())
        }
        "bond-order" => bond_order(r, plots),
        "wigner-seitz" => {
            add_histogram(plots, "Distance to assigned site", "Distance (A)", &numbers(r.find("distance_A")));
            Ok(())
        }
        "ptm" => ptm(r, plots),
        "voronoi" => voronoi(r, plots),
        "bond-orders" => bond_order_assignment(r, plots),
        "scatter-plot" => scatter(r, plots),
        "ambient-occlusion" => {
            add_histogram(plots, "Ambient brightness", "Brightness", &numbers(r.find("brightness")));
            Ok(())
        }
        "surface-mesh" => surface_mesh(r, plots),
        "trajectory-average" => {
            add_histogram(plots, "RMS fluctuation", "RMS fluctuation (A)", &numbers(r.find("rms_fluctuation_A")));
            Ok(())
        }
        "unwrap-trajectory" => {
            add_histogram(plots, "Net displacement", "|r(end) - r(start)| (A)", &numbers(r.find("net_displacement_A")));
            Ok(())
        }
        "time-series" => time_series(r, plots),
        "frame-difference" => {
            let label = match r.find("quantity") {
                Some(Json::String(q)) => q.clone(),
                _ => "Change".into(),
            };
            add_histogram(plots, "Change between frames", &label, &numbers(r.find("difference")));
            Ok(())
        }
        "cna" => count_plot(r, plots, "CNA structures (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)"),
        "vorotop" => count_plot(r, plots, "VoroTop structures (type index as in the report)"),
        "grain-segmentation" => grain_segmentation(r, plots),
        "elastic-strain" => {
            add_histogram(plots, "Volumetric elastic strain", "trace(E)/3", &numbers(r.find("volumetric_strain")));
            add_histogram(plots, "Elastic shear strain", "Von Mises shear strain", &numbers(r.find("shear_strain")));
            Ok(())
        }
        "fcc-planar-faults" => count_plot(r, plots, "Planar faults (0 non-hcp, 1 stacking fault, 2 twin, 3 multilayer, 4 other)"),
        "chill-plus" => count_plot(r, plots, "Chill+ (0 other, 1 hex ice, 2 cubic ice, 3 interf. ice, 4 hydrate, 5 interf. hydrate)"),
        "local-entropy" => {
            add_histogram(plots, "Local entropy", "s (kB)", &numbers(r.find("entropy_kB")));
            Ok(())
        }
        "find-rings" => {
            plots.push(plot("Ring-size distribution", "Ring size (atoms)", "Rings", vec![line("rings", numbers(r.find("ring_sizes")), numbers(r.find("ring_counts")), true)]));
            Ok(())
        }
        "bond-distributions" => bond_distributions(r, plots),
        "displacements" => {
            add_histogram(plots, "Displacement magnitudes", "|u| (A)", &numbers(r.find("magnitude_A")));
            Ok(())
        }
        "spatial-binning" => spatial_binning(r, plots),
        "property-histogram" => {
            plots.push(plot("Histogram", "Value", "Count", vec![line("count", numbers(r.find("centers")), numbers(r.find("counts")), false)]));
            Ok(())
        }
        "spatial-correlation" => {
            let mut spec = plot("Spatial correlation", "Distance (A)", "Normalized correlation", vec![line("C(r)", numbers(r.find("r_A")), numbers(r.find("normalized_correlation")), false)]);
            spec.horizontal = vec![0.0];
            plots.push(spec);
            Ok(())
        }
        "work-function" => work_function(r, plots),
        "equation-of-state" | "eos-scan" => equation_of_state(r, plots),
        "trajectory-structure" => trajectory_structure(r, plots),
        "phonon-dos" => phonon_dos(r, plots),
        "harmonic-thermodynamics" => harmonic_thermodynamics(r, plots),
        "phonons" => phonons(r, plots),
        "relax" => relax(r, plots),
        "neb" => neb(r, plots),
        "nvt" => dynamics(r, plots, false),
        "npt" => dynamics(r, plots, true),
        _ => Ok(()),
    }
}

/// Tool ids that have result plots, sorted (the C++ std::map order).
pub fn plotted_tools() -> Vec<&'static str> {
    let mut tools = vec![
        "msd", "diffusion", "vacf", "vibrational-spectrum", "structure-factor", "local-strain", "vasp-electronic", "powder-xrd", "electron-diffraction", "cluster-analysis", "void-analysis", "structure-type",
        "diamond-structure", "centrosymmetry", "bond-order", "wigner-seitz", "ptm", "voronoi", "bond-orders", "scatter-plot", "ambient-occlusion", "surface-mesh", "trajectory-average", "unwrap-trajectory",
        "time-series", "frame-difference", "cna", "vorotop", "grain-segmentation", "elastic-strain", "fcc-planar-faults", "chill-plus", "local-entropy", "find-rings", "bond-distributions", "displacements",
        "spatial-binning", "property-histogram", "spatial-correlation", "work-function", "equation-of-state", "eos-scan", "trajectory-structure", "phonon-dos", "harmonic-thermodynamics", "phonons", "relax",
        "neb", "nvt", "npt",
    ];
    tools.sort_unstable();
    tools
}

/// Fixed-width histogram of the finite values (bin centres and counts).
pub fn histogram(name: &str, values: &[f64], bins: i32) -> PlotSeries {
    let mut series = PlotSeries { name: name.to_string(), ..Default::default() };
    let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
    for &v in values {
        if v.is_finite() {
            low = low.min(v);
            high = high.max(v);
        }
    }
    if !(low <= high) || bins < 1 {
        return series;
    }
    if high == low {
        high = low + 1e-12 + low.abs() * 1e-9;
    }
    let width = (high - low) / bins as f64;
    let mut counts = vec![0.0; bins as usize];
    for &v in values {
        if v.is_finite() {
            let b = (bins - 1).min(((v - low) / width) as i32);
            counts[b as usize] += 1.0;
        }
    }
    for (b, &c) in counts.iter().enumerate() {
        series.x.push(low + (b as f64 + 0.5) * width);
        series.y.push(c);
    }
    series
}

/// Plots describing a tool result; empty when the result has no curve to
/// show. Plots without two finite points are left out.
pub fn result_plots(tool: &str, result: &Json) -> Result<Vec<PlotSpec>, String> {
    let mut plots = Vec::new();
    if !result.is_object() {
        return Ok(plots);
    }
    build(tool, result, &mut plots)?;
    plots.retain(|spec| spec.series.iter().map(|s| s.x.iter().zip(&s.y).filter(|(x, y)| x.is_finite() && y.is_finite()).count()).sum::<usize>() >= 2);
    Ok(plots)
}

/// One plot with the matching plots (same title) of earlier runs drawn on
/// the same axes; their series are named "label: series".
pub fn overlay_plots(current: &PlotSpec, current_label: &str, others: &[(String, Vec<PlotSpec>)]) -> PlotSpec {
    let mut combined = current.clone();
    combined.series.clear();
    let mut add = |label: &str, spec: &PlotSpec| {
        for s in &spec.series {
            let mut s = s.clone();
            s.name = format!("{label}: {}", s.name);
            combined.series.push(s);
        }
    };
    for (label, plots) in others {
        if let Some(spec) = plots.iter().find(|p| p.title == current.title) {
            add(label, spec);
        }
    }
    add(current_label, current);
    combined
}

/// CSV text of a plot's series (one x/y column pair per series).
pub fn plot_csv(plot: &PlotSpec) -> String {
    let mut text = String::new();
    let mut rows = 0;
    for (s, series) in plot.series.iter().enumerate() {
        if s > 0 {
            text.push(',');
        }
        text += &format!("{} [{}],{} [{}]", plot.x_label, series.name, plot.y_label, series.name);
        rows = rows.max(series.x.len());
    }
    text.push('\n');
    for r in 0..rows {
        for (s, series) in plot.series.iter().enumerate() {
            if s > 0 {
                text.push(',');
            }
            if r < series.x.len() && r < series.y.len() {
                text += &format!("{},{}", af_core::format::format_g(series.x[r], 10), af_core::format::format_g(series.y[r], 10));
            } else {
                text.push(',');
            }
        }
        text.push('\n');
    }
    text
}
