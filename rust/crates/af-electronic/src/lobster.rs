//! LOBSTER bonding output (the C++ electronic/Lobster): COHPCAR/COOPCAR/COBICAR
//! curves and ICOHPLIST/ICOOPLIST/ICOBILIST integrals, per-bond totals only.

use crate::text::{lines_of, stod_prefix, stoi_prefix, words};

/// Zero-based site indices.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CohpBond {
    pub atom1: i32,
    pub atom2: i32,
    pub length: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CohpData {
    pub spin_polarized: bool,
    pub fermi_energy: f64,
    /// Fermi level at 0 eV.
    pub energies: Vec<f64>,
    /// Excludes the "average" entry.
    pub bonds: Vec<CohpBond>,
    /// [spin][energy].
    pub average_cohp: Vec<Vec<f64>>,
    pub average_icohp: Vec<Vec<f64>>,
    /// [spin][bond][energy].
    pub cohp: Vec<Vec<Vec<f64>>>,
    pub icohp: Vec<Vec<Vec<f64>>>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IcohpEntry {
    pub atom1: i32,
    pub atom2: i32,
    pub length: f64,
    pub num_bonds: i32,
    pub icohp: [f64; 2],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct IcohpList {
    pub spin_polarized: bool,
    pub entries: Vec<IcohpEntry>,
}

fn read_all_lines(path: &str) -> Result<Vec<String>, String> {
    lines_of(path).ok_or_else(|| "Cannot open LOBSTER output file".to_string())
}

fn blank(s: &str) -> bool {
    words(s).is_empty()
}

/// "Fe12" -> 11 (zero-based site index).
fn site_index_from_label(label: &str) -> Result<i32, String> {
    let b = label.as_bytes();
    let mut start = b.len();
    while start > 0 && b[start - 1].is_ascii_digit() {
        start -= 1;
    }
    if start == b.len() {
        return Err(format!("LOBSTER site label has no trailing index: {label}"));
    }
    Ok(stoi_prefix(&label[start..])? - 1)
}

/// std::getline(in, part, ':') pieces.
fn colon_parts(s: &str) -> Vec<&str> {
    let mut parts: Vec<&str> = s.split(':').collect();
    if s.is_empty() || s.ends_with(':') {
        parts.pop();
    }
    parts
}

pub fn read_cohpcar(path: &str) -> Result<CohpData, String> {
    let lines = read_all_lines(path)?;
    if lines.len() < 4 {
        return Err("Truncated COHPCAR file".into());
    }
    let parameters = words(&lines[1]);
    if parameters.len() < 3 {
        return Err("Invalid COHPCAR parameter line".into());
    }
    let total_bonds = stoi_prefix(parameters[0])?; // includes the "average" entry
    let spin_flag = stoi_prefix(parameters[1])?;
    if total_bonds < 1 || (spin_flag != 1 && spin_flag != 2) {
        return Err("Invalid COHPCAR bond count or spin flag".into());
    }
    let n = (total_bonds - 1) as usize;
    let num_spins = spin_flag as usize;
    let mut result = CohpData { spin_polarized: num_spins == 2, fermi_energy: stod_prefix(parameters[parameters.len() - 1])?, ..Default::default() };
    if lines.len() < 3 + n {
        return Err("Truncated COHPCAR bond header block".into());
    }
    // Only bond totals are kept (headers without '['); raw_index is each kept
    // bond's position among all headers, which the data columns are keyed on.
    let mut raw_index = Vec::new();
    for bond in 0..n {
        let header = &lines[3 + bond];
        if header.contains('[') {
            continue;
        }
        let open = match header.rfind('(') {
            Some(open) if header.ends_with(')') => open,
            _ => return Err(format!("Invalid COHPCAR bond header: {header}")),
        };
        let length = stod_prefix(&header[open + 1..header.len() - 1])?;
        let prefix = header[..open].replace("->", ":");
        let parts = colon_parts(&prefix);
        if parts.len() < 3 {
            return Err(format!("Invalid COHPCAR bond header: {header}"));
        }
        result.bonds.push(CohpBond { atom1: site_index_from_label(parts[1])?, atom2: site_index_from_label(parts[2])?, length });
        raw_index.push(bond);
    }
    if result.bonds.is_empty() {
        return Err("COHPCAR file contains no bond-total headers".into());
    }
    let stored = result.bonds.len();
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for line in &lines[3 + n..] {
        if blank(line) {
            continue;
        }
        rows.push(words(line).iter().map(|t| stod_prefix(t)).collect::<Result<_, _>>()?);
    }
    if rows.is_empty() {
        return Err("COHPCAR file contains no data rows".into());
    }
    let expected = 1 + 2 * (n + 1) * num_spins;
    let count = rows.len();
    result.energies = vec![0.0; count];
    result.average_cohp = vec![vec![0.0; count]; num_spins];
    result.average_icohp = vec![vec![0.0; count]; num_spins];
    result.cohp = vec![vec![vec![0.0; count]; stored]; num_spins];
    result.icohp = vec![vec![vec![0.0; count]; stored]; num_spins];
    for (r, row) in rows.iter().enumerate() {
        if row.len() != expected {
            return Err("Inconsistent COHPCAR data row width".into());
        }
        result.energies[r] = row[0];
        for s in 0..num_spins {
            result.average_cohp[s][r] = row[1 + 2 * s * (n + 1)];
            result.average_icohp[s][r] = row[2 + 2 * s * (n + 1)];
            for b in 0..stored {
                result.cohp[s][b][r] = row[2 * (raw_index[b] + s * (n + 1)) + 3];
                result.icohp[s][b][r] = row[2 * (raw_index[b] + s * (n + 1)) + 4];
            }
        }
    }
    Ok(result)
}

pub fn read_icohplist(path: &str) -> Result<IcohpList, String> {
    let mut lines = read_all_lines(path)?;
    if lines.len() < 2 {
        return Err("Truncated ICOHPLIST file".into());
    }
    lines.remove(0); // the column-label header line
    while lines.last().is_some_and(|l| blank(l)) {
        lines.pop();
    }
    if lines.is_empty() {
        return Err("ICOHPLIST file contains no data".into());
    }
    let spin_polarized = lines[lines.len() / 2].contains("distance");
    // Only the bond totals are kept (second column without '_').
    let totals: Vec<&String> = lines.iter().filter(|l| !words(l).get(1).is_some_and(|t| t.contains('_'))).collect();
    if totals.is_empty() {
        return Err("ICOHPLIST file contains no bond-total rows".into());
    }
    let width = match words(totals[0]).len() {
        8 => 8, // LOBSTER >= 3.1.1 (with a translation vector)
        6 => 6, // LOBSTER <= 2.2.1
        _ => return Err("Unrecognized ICOHPLIST column count".into()),
    };
    let totals_spin_polarized = totals[totals.len() / 2].contains("distance");
    let num_bonds = if totals_spin_polarized { totals.len() / 2 } else { totals.len() };
    if num_bonds == 0 {
        return Err("ICOHPLIST file contains no data".into());
    }
    let mut result = IcohpList { spin_polarized, entries: Vec::new() };
    for bond in 0..num_bonds {
        let tokens = words(totals[bond]);
        if tokens.len() != width {
            return Err("Inconsistent ICOHPLIST row width".into());
        }
        let mut entry = IcohpEntry {
            atom1: site_index_from_label(tokens[1])?,
            atom2: site_index_from_label(tokens[2])?,
            length: stod_prefix(tokens[3])?,
            num_bonds: 1,
            icohp: [0.0; 2],
        };
        if width == 8 {
            entry.icohp[0] = stod_prefix(tokens[7])?;
        } else {
            entry.icohp[0] = stod_prefix(tokens[4])?;
            entry.num_bonds = stoi_prefix(tokens[5])?;
        }
        if spin_polarized {
            let down = totals.get(bond + num_bonds + 1).map(|l| words(l)).unwrap_or_default();
            if down.len() != width {
                return Err("Inconsistent ICOHPLIST spin-down row width".into());
            }
            entry.icohp[1] = stod_prefix(down[if width == 8 { 7 } else { 4 }])?;
        }
        result.entries.push(entry);
    }
    Ok(result)
}
