//! CIF 1.1 / CIF 2 / magnetic CIF (the C++ readCif): data blocks, items and
//! loops; the first block with atom sites, its cell, its symmetry operations
//! (listed, or from the Hall symbol, H-M name or number via spglib's table)
//! and the symmetry images of each site.

use std::collections::BTreeMap;

use af_core::{elements, AtomSite, Structure};
use glam::{DMat3, DVec3};

use crate::occupancy::resolve_partial_occupancy;
use crate::spacegroups;
use crate::text::lines;

#[derive(Default)]
struct Loop {
    tags: Vec<String>,
    values: Vec<String>,
}

impl Loop {
    fn column(&self, tag: &str) -> Option<usize> {
        self.tags.iter().position(|t| t == tag)
    }
}

#[derive(Default)]
struct Block {
    name: String,
    items: BTreeMap<String, String>,
    loops: Vec<Loop>,
}

impl Block {
    fn loop_with(&self, tag: &str) -> Option<&Loop> {
        self.loops.iter().find(|l| l.column(tag).is_some())
    }
    fn item(&self, tag: &str) -> String {
        if let Some(v) = self.items.get(tag) {
            return v.clone();
        }
        if let Some(l) = self.loop_with(tag) {
            if l.values.len() == l.tags.len() {
                return l.values[l.column(tag).unwrap()].clone();
            }
        }
        String::new()
    }
}

struct Token {
    text: String,
    quoted: bool,
}

fn is_space(c: u8) -> bool {
    crate::text::is_space(c as char)
}

fn parse(text: &str) -> Vec<Block> {
    let all = lines(text, true);
    let mut tokens: Vec<Token> = Vec::new();
    let mut n = 0;
    while n < all.len() {
        let line = all[n];
        n += 1;
        if line.starts_with(';') {
            // A text field, up to the next line starting with ';'.
            let mut field = line[1..].to_string();
            while n < all.len() {
                let l = all[n];
                n += 1;
                if l.starts_with(';') {
                    break;
                }
                field.push('\n');
                field.push_str(l);
            }
            tokens.push(Token { text: field, quoted: true });
            continue;
        }
        let mut line = line.to_string();
        let mut i = 0;
        loop {
            let b = line.as_bytes();
            while i < b.len() && is_space(b[i]) {
                i += 1;
            }
            if i >= b.len() || b[i] == b'#' {
                break;
            }
            if line[i..].starts_with("\"\"\"") || line[i..].starts_with("'''") {
                // CIF 2 triple-quoted string, across lines if needed.
                let delimiter = line[i..i + 3].to_string();
                let mut field = String::new();
                let mut from = i + 3;
                let mut close = line[from..].find(&delimiter).map(|c| from + c);
                let mut ended = false;
                while close.is_none() {
                    field.push_str(&line[from..]);
                    field.push('\n');
                    if n >= all.len() {
                        ended = true;
                        break;
                    }
                    line = all[n].to_string();
                    n += 1;
                    from = 0;
                    close = line.find(&delimiter);
                }
                if ended || close.is_none() {
                    tokens.push(Token { text: field, quoted: true });
                    break;
                }
                let c = close.unwrap();
                field.push_str(&line[from..c]);
                tokens.push(Token { text: field, quoted: true });
                i = c + 3;
                continue;
            }
            if b[i] == b'\'' || b[i] == b'"' {
                let quote = b[i];
                let mut close = i + 1;
                while close < b.len() && !(b[close] == quote && (close + 1 == b.len() || is_space(b[close + 1]))) {
                    close += 1;
                }
                let end = close.min(b.len());
                tokens.push(Token { text: line[i + 1..end].to_string(), quoted: true });
                i = close + 1;
                continue;
            }
            let start = i;
            while i < b.len() && !is_space(b[i]) {
                i += 1;
            }
            tokens.push(Token { text: line[start..i].to_string(), quoted: false });
        }
    }
    let tag_of = |v: &str| v.to_ascii_lowercase().replace('.', "_");
    let keyword = |t: &Token| {
        if t.quoted {
            return false;
        }
        let w = t.text.to_ascii_lowercase();
        w.starts_with("data_") || w == "loop_" || w.starts_with("save_") || w == "global_" || w == "stop_"
    };
    let is_tag = |t: &Token| !t.quoted && t.text.starts_with('_');
    let mut blocks: Vec<Block> = Vec::new();
    let mut i = 0;
    while i < tokens.len() {
        let token = &tokens[i];
        let word = if token.quoted { String::new() } else { token.text.to_ascii_lowercase() };
        if word.starts_with("data_") {
            blocks.push(Block { name: token.text[5..].to_string(), ..Default::default() });
            i += 1;
        } else if word == "loop_" {
            let mut l = Loop::default();
            i += 1;
            while i < tokens.len() && is_tag(&tokens[i]) {
                l.tags.push(tag_of(&tokens[i].text));
                i += 1;
            }
            while i < tokens.len() && !is_tag(&tokens[i]) && !keyword(&tokens[i]) {
                l.values.push(tokens[i].text.clone());
                i += 1;
            }
            if !l.tags.is_empty() {
                if blocks.is_empty() {
                    blocks.push(Block::default());
                }
                blocks.last_mut().unwrap().loops.push(l);
            }
        } else if is_tag(token) {
            let tag = tag_of(&token.text);
            i += 1;
            if i < tokens.len() && !is_tag(&tokens[i]) && !keyword(&tokens[i]) {
                if blocks.is_empty() {
                    blocks.push(Block::default());
                }
                blocks.last_mut().unwrap().items.insert(tag, tokens[i].text.clone());
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    blocks
}

/// A CIF number without its uncertainty; None for '?' and '.'.
fn cif_number(value: &str) -> Option<f64> {
    let text = &value[..value.find('(').unwrap_or(value.len())];
    if text.is_empty() || text == "?" || text == "." {
        return None;
    }
    crate::text::number_full(text)
}

/// A symmetry operation such as "-x+1/2, y, z" or "1/2+x,-y,z+0.5".
pub fn parse_symop(text: &str) -> Option<(DMat3, DVec3)> {
    let s: Vec<u8> = text.bytes().filter(|&c| !is_space(c) && c != b'\'' && c != b'"').map(|c| c.to_ascii_lowercase()).collect();
    let st = String::from_utf8_lossy(&s).to_string();
    let mut m = [[0.0f64; 3]; 3]; // m[row][column]
    let mut t = [0.0f64; 3];
    let mut row = 0usize;
    let mut i = 0usize;
    while i <= s.len() {
        if i == s.len() || s[i] == b',' {
            row += 1;
            i += 1;
            if row > 3 {
                return None;
            }
            continue;
        }
        if row > 2 {
            return None;
        }
        let mut sign = 1.0;
        while i < s.len() && (s[i] == b'+' || s[i] == b'-') {
            if s[i] == b'-' {
                sign = -sign;
            }
            i += 1;
        }
        if i >= s.len() {
            return None;
        }
        if (b'x'..=b'z').contains(&s[i]) {
            m[row][(s[i] - b'x') as usize] += sign;
            i += 1;
            continue;
        }
        if !s[i].is_ascii_digit() && s[i] != b'.' {
            return None;
        }
        let (mut value, used) = af_core::format::parse_prefix_f64(&st[i..]).unwrap_or((0.0, 0));
        i += used;
        if i < s.len() && s[i] == b'/' {
            i += 1;
            let (denominator, used) = af_core::format::parse_prefix_f64(&st[i..]).unwrap_or((0.0, 0));
            if used == 0 || denominator == 0.0 {
                return None;
            }
            i += used;
            value /= denominator;
        }
        if i < s.len() && s[i] == b'*' {
            i += 1;
        }
        if i < s.len() && (b'x'..=b'z').contains(&s[i]) {
            m[row][(s[i] - b'x') as usize] += sign * value;
            i += 1;
        } else {
            t[row] += sign * value;
        }
    }
    if row != 3 {
        return None;
    }
    let rotation = DMat3::from_cols(DVec3::new(m[0][0], m[1][0], m[2][0]), DVec3::new(m[0][1], m[1][1], m[2][1]), DVec3::new(m[0][2], m[1][2], m[2][2]));
    Some((rotation, DVec3::from(t)))
}

/// The element of a type symbol or label: Fe3+ and fe are Fe, CA2 is Ca,
/// OW1 is O, D and T are hydrogen.
fn cif_element(text: &str) -> i32 {
    let b = text.as_bytes();
    let letters = b.iter().take_while(|c| c.is_ascii_alphabetic()).count();
    if letters == 0 {
        return 0;
    }
    let first = (b[0] as char).to_ascii_uppercase();
    if letters >= 2 {
        let z = elements::atomic_number(&format!("{first}{}", (b[1] as char).to_ascii_lowercase()));
        if z > 0 {
            return z;
        }
    }
    if first == 'D' || first == 'T' {
        return 1;
    }
    elements::atomic_number(&first.to_string())
}

fn normalized_group_name(name: &str) -> String {
    name.chars().filter(|c| !crate::text::is_space(*c) && *c != '_' && *c != '\'' && *c != '"').map(|c| c.to_ascii_lowercase()).collect()
}

/// Hall numbers for a Hall symbol, H-M name or number, and the requested setting.
fn hall_numbers_for(hall_symbol: &str, hm_in: &str, number: i32, choice: &mut String) -> Vec<i32> {
    let types = spacegroups::types();
    if !hall_symbol.is_empty() {
        let wanted = normalized_group_name(hall_symbol);
        if let Some(t) = types.iter().find(|t| normalized_group_name(&t.hall_symbol) == wanted) {
            return vec![t.hall_number];
        }
    }
    let mut found = Vec::new();
    if !hm_in.is_empty() {
        let mut hm = hm_in.to_string();
        if let Some(colon) = hm.find(':') {
            *choice = normalized_group_name(&hm[colon + 1..]);
            hm.truncate(colon);
        } else {
            let trimmed = hm.trim_end_matches(|c: char| crate::text::is_space(c)).to_string();
            if let Some(space) = trimmed.rfind([' ', '\t']) {
                let last = trimmed.as_bytes()[trimmed.len() - 1];
                if trimmed.len() == space + 2 && b"ZSHR".contains(&last) {
                    *choice = (last as char).to_ascii_lowercase().to_string();
                    hm = trimmed[..space].to_string();
                }
            }
        }
        if choice == "z" {
            *choice = "2".into();
        }
        if choice == "s" {
            *choice = "1".into();
        }
        let wanted = normalized_group_name(&hm);
        for t in types {
            let mut names = vec![t.international_short.clone(), t.international_full.clone()];
            {
                let axes: Vec<&str> = t.international_full.split_whitespace().collect();
                let kept: Vec<&str> = axes.into_iter().filter(|a| *a != "1").collect();
                if kept.len() > 1 {
                    names.push(kept.iter().map(|a| format!("{a} ")).collect());
                }
            }
            let mut full = t.international.clone();
            while let Some(at) = full.find('=') {
                names.push(full[..at].to_string());
                full = full[at + 1..].to_string();
            }
            names.push(full);
            if names.iter().any(|n| normalized_group_name(n) == wanted) {
                found.push(t.hall_number);
            }
        }
    }
    if found.is_empty() && (1..=230).contains(&number) {
        found = types.iter().filter(|t| t.number == number).map(|t| t.hall_number).collect();
    }
    if !choice.is_empty() {
        let chosen: Vec<i32> = found.iter().copied().filter(|&h| normalized_group_name(&spacegroups::of_hall(h).choice) == *choice).collect();
        if !chosen.is_empty() {
            found = chosen;
        }
    }
    found
}

struct Site {
    atomic_number: i32,
    fractional: DVec3,
    occupancy: f64,
}

fn floor3(v: DVec3) -> DVec3 {
    DVec3::new(v.x.floor(), v.y.floor(), v.z.floor())
}

fn round3(v: DVec3) -> DVec3 {
    // glm::round rounds halves away from zero, as f64::round.
    DVec3::new(v.x.round(), v.y.round(), v.z.round())
}

/// Every symmetry image of each site in the cell; images of one site within
/// 0.1 A of each other are one atom at their mean position.
fn expand_sites(sites: &[Site], operations: &[(DMat3, DVec3)], cell: &DMat3) -> Vec<(DVec3, usize)> {
    let mut atoms = Vec::new();
    for (index, site) in sites.iter().enumerate() {
        struct Image {
            base: DVec3,
            sum: DVec3,
            count: i32,
        }
        let mut images: Vec<Image> = Vec::new();
        for (rotation, translation) in operations {
            let mut f = *rotation * site.fractional + *translation;
            f -= floor3(f);
            let mut merged = false;
            for image in images.iter_mut() {
                let mut d = f - image.base;
                d -= round3(d);
                if (*cell * d).length() < 0.1 {
                    image.sum += d;
                    image.count += 1;
                    merged = true;
                    break;
                }
            }
            if !merged {
                images.push(Image { base: f, sum: DVec3::ZERO, count: 1 });
            }
        }
        for image in images {
            let mut f = image.base + image.sum / image.count as f64;
            f -= floor3(f);
            atoms.push((f, index));
        }
    }
    atoms
}

/// The cell of lengths and angles (degrees) with a along x and b in the xy plane.
pub fn cell_from_parameters(length: [f64; 3], angle: [f64; 3]) -> Option<DMat3> {
    let cosine = |degrees: f64| {
        if (degrees - 90.0).abs() < 1e-10 {
            0.0
        } else if (degrees - 60.0).abs() < 1e-10 {
            0.5
        } else if (degrees - 120.0).abs() < 1e-10 {
            -0.5
        } else {
            (degrees * std::f64::consts::PI / 180.0).cos()
        }
    };
    let (ca, cb, cg) = (cosine(angle[0]), cosine(angle[1]), cosine(angle[2]));
    let sg = (1.0 - cg * cg).sqrt();
    if sg < 1e-9 {
        return None;
    }
    let cy = (ca - cb * cg) / sg;
    let cz2 = 1.0 - cb * cb - cy * cy;
    if cz2 <= 1e-12 {
        return None;
    }
    Some(DMat3::from_cols(
        DVec3::new(length[0], 0.0, 0.0),
        DVec3::new(length[1] * cg, length[1] * sg, 0.0),
        DVec3::new(length[2] * cb, length[2] * cy, length[2] * cz2.sqrt()),
    ))
}

/// Ok(true) when read, Ok(false) when the file has no atom sites, Err on an
/// error. Sites with only Cartesian coordinates (_atom_site_Cartn_x, which
/// the C++ handed to Open Babel) are converted to fractional ones in the
/// CIF's standard orthogonal frame (a along x, b in the xy plane).
pub fn read(text: &str, structure: &mut Structure) -> Result<bool, String> {
    let blocks = parse(text);
    let sites_of = |b: &Block| b.loop_with("_atom_site_fract_x").is_some() || b.loop_with("_atom_site_cartn_x").is_some();
    let with_atoms: Vec<&Block> = blocks.iter().filter(|b| sites_of(b)).collect();
    let Some(&block) = with_atoms.first() else { return Ok(false) };
    let fail = |message: String| Err(format!("Invalid CIF: {message}"));
    let mut length = [0.0; 3];
    let mut angle = [90.0; 3];
    let length_tags = ["_cell_length_a", "_cell_length_b", "_cell_length_c"];
    let angle_tags = ["_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"];
    for k in 0..3 {
        match cif_number(&block.item(length_tags[k])) {
            Some(v) if v > 0.0 => length[k] = v,
            _ => return fail(format!("the cell length {} is missing.", length_tags[k])),
        }
        let value = block.item(angle_tags[k]);
        if !value.is_empty() {
            match cif_number(&value) {
                Some(v) if v > 0.0 && v < 180.0 => angle[k] = v,
                _ => return fail(format!("{} is not an angle.", angle_tags[k])),
            }
        }
    }
    let Some(cell) = cell_from_parameters(length, angle) else { return fail("the cell angles do not make a cell.".into()) };

    let cartesian = block.loop_with("_atom_site_fract_x").is_none();
    let l = block.loop_with(if cartesian { "_atom_site_cartn_x" } else { "_atom_site_fract_x" }).unwrap();
    let to_fractional = cell.inverse();
    let columns = l.tags.len();
    if l.values.len() % columns != 0 {
        return fail("the atom-site loop has an incomplete row.".into());
    }
    let label = l.column("_atom_site_label");
    let type_column = l.column("_atom_site_type_symbol");
    let axis = |k: &str| l.column(&format!("_atom_site_{}_{k}", if cartesian { "cartn" } else { "fract" }));
    let (fx, fy, fz) = (axis("x"), axis("y"), axis("z"));
    let occupancy_column = l.column("_atom_site_occupancy");
    let flag = l.column("_atom_site_calc_flag");
    if fy.is_none() || fz.is_none() {
        return fail(if cartesian { "the atom sites need Cartn_x, Cartn_y and Cartn_z." } else { "the atom sites need fract_x, fract_y and fract_z." }.into());
    }
    if label.is_none() && type_column.is_none() {
        return fail("the atom sites need a label or a type symbol.".into());
    }
    let mut sites = Vec::new();
    let (mut unplaced, mut empty) = (0usize, 0usize);
    let mut row = 0;
    while row < l.values.len() {
        let value = |c: Option<usize>| c.map_or(String::new(), |c| l.values[row + c].clone());
        let name = if value(label).is_empty() { value(type_column) } else { value(label) };
        if value(flag) == "dum" {
            row += columns;
            continue;
        }
        let (Some(x), Some(y), Some(z)) = (cif_number(&value(fx)), cif_number(&value(fy)), cif_number(&value(fz))) else {
            unplaced += 1;
            row += columns;
            continue;
        };
        let mut z_number = cif_element(&value(type_column));
        if z_number == 0 {
            z_number = cif_element(&value(label));
        }
        if z_number == 0 {
            return fail(format!("the element of atom site '{name}' is unknown."));
        }
        let mut occupancy = 1.0;
        if let Some(o) = cif_number(&value(occupancy_column)) {
            if o <= 0.0 {
                empty += 1;
                row += columns;
                continue;
            }
            occupancy = o;
        }
        let fractional = if cartesian { to_fractional * DVec3::new(x, y, z) } else { DVec3::new(x, y, z) };
        sites.push(Site { atomic_number: z_number, fractional, occupancy });
        row += columns;
    }

    let mut operations: Vec<(DMat3, DVec3)> = Vec::new();
    let mut symops = block.loop_with("_space_group_symop_operation_xyz");
    let mut symop_tag = "_space_group_symop_operation_xyz";
    if symops.is_none() {
        symops = block.loop_with("_symmetry_equiv_pos_as_xyz");
        symop_tag = "_symmetry_equiv_pos_as_xyz";
    }
    let column_values = |source: Option<&Loop>, tag: &str| -> Vec<String> {
        let mut v = Vec::new();
        if let Some(s) = source {
            let mut r = s.column(tag).unwrap();
            while r < s.values.len() {
                v.push(s.values[r].clone());
                r += s.tags.len();
            }
        }
        v
    };
    let spatial = |op: &str| -> String {
        if op.matches(',').count() == 3 {
            op[..op.rfind(',').unwrap()].to_string()
        } else {
            op.to_string()
        }
    };
    let magnetic = column_values(block.loop_with("_space_group_symop_magn_operation_xyz"), "_space_group_symop_magn_operation_xyz");
    let mut centerings = column_values(block.loop_with("_space_group_symop_magn_centering_xyz"), "_space_group_symop_magn_centering_xyz");
    let mut symop_texts: Vec<String> = Vec::new();
    if !magnetic.is_empty() {
        if centerings.is_empty() {
            centerings.push("x,y,z".into());
        }
        for centering in &centerings {
            let Some((sr, st)) = parse_symop(&spatial(centering)) else { return fail(format!("cannot read the centering '{centering}'.")) };
            for op in &magnetic {
                let Some((br, bt)) = parse_symop(&spatial(op)) else { return fail(format!("cannot read the symmetry operation '{op}'.")) };
                operations.push((sr * br, sr * bt + st));
            }
        }
    } else if let Some(s) = symops {
        symop_texts = column_values(Some(s), symop_tag);
    } else {
        let single = if let Some(v) = block.items.get("_space_group_symop_operation_xyz") {
            v.clone()
        } else {
            block.item("_symmetry_equiv_pos_as_xyz")
        };
        if !single.is_empty() {
            symop_texts.push(single);
        }
    }
    if magnetic.is_empty() {
        for symop in &symop_texts {
            let Some(op) = parse_symop(symop) else { return fail(format!("cannot read the symmetry operation '{symop}'.")) };
            operations.push(op);
        }
    }
    let mut note = String::new();
    if operations.is_empty() {
        let mut hall = block.item("_space_group_name_hall");
        if hall.is_empty() {
            hall = block.item("_symmetry_space_group_name_hall");
        }
        let mut hm = block.item("_space_group_name_h-m_alt");
        if hm.is_empty() {
            hm = block.item("_symmetry_space_group_name_h-m");
        }
        let mut number_text = block.item("_space_group_it_number");
        if number_text.is_empty() {
            number_text = block.item("_symmetry_int_tables_number");
        }
        let number = cif_number(&number_text).map_or(0, |v| v as i32);
        let named = !hall.is_empty() || (!hm.is_empty() && hm != "?") || number > 1;
        if named {
            let mut choice = String::new();
            let halls = hall_numbers_for(if hall == "?" { "" } else { &hall }, if hm == "?" { "" } else { &hm }, number, &mut choice);
            if halls.is_empty() {
                let shown = if !hm.is_empty() {
                    hm.clone()
                } else if !hall.is_empty() {
                    hall.clone()
                } else {
                    number_text.clone()
                };
                return fail(format!(
                    "the space group '{shown}' is not recognised; list its symmetry operations (_space_group_symop_operation_xyz)."
                ));
            }
            let choice_of = |h: i32| normalized_group_name(&spacegroups::of_hall(h).choice);
            let with_choice = |wanted: &str| halls.iter().copied().find(|&h| choice_of(h) == wanted).unwrap_or(0);
            let mut hall_number = halls[0];
            if halls.len() > 1 && with_choice("h") != 0 && with_choice("r") != 0 {
                hall_number = if (angle[2] - 120.0).abs() < 0.5 && (length[0] - length[1]).abs() < 1e-3 * length[0] {
                    with_choice("h")
                } else {
                    with_choice("r")
                };
            } else if halls.len() > 1 && with_choice("1") != 0 && with_choice("2") != 0 {
                let one = expand_sites(&sites, &spacegroups::of_hall(with_choice("1")).operations, &cell).len();
                let two = expand_sites(&sites, &spacegroups::of_hall(with_choice("2")).operations, &cell).len();
                hall_number = if one < two { with_choice("1") } else { with_choice("2") };
                if one == two {
                    note = "The file does not say which origin choice its space group uses; origin choice 2 is assumed.".into();
                }
            }
            operations = spacegroups::of_hall(hall_number).operations.clone();
        }
    }
    if operations.is_empty() {
        operations.push((DMat3::IDENTITY, DVec3::ZERO));
    }
    let atoms = expand_sites(&sites, &operations, &cell);
    if atoms.len() > 5_000_000 {
        return fail("it has more than 5000000 atoms.".into());
    }
    let mut occupancies = Vec::with_capacity(atoms.len());
    for (fractional, index) in &atoms {
        let z = sites[*index].atomic_number;
        let mut site = AtomSite::new(elements::symbol(z), (cell * *fractional).to_array());
        site.atomic_number = z;
        structure.atoms.push(site);
        occupancies.push(sites[*index].occupancy);
    }
    structure.has_unit_cell = true;
    for r in 0..3 {
        structure.cell_vectors[r] = cell.col(r).to_array();
    }
    resolve_partial_occupancy(structure, &occupancies);
    let mut add = |more: String| {
        if !note.is_empty() {
            note.push(' ');
        }
        note.push_str(&more);
    };
    if !structure.load_note.is_empty() {
        add(structure.load_note.clone());
    }
    if unplaced > 0 {
        add(format!(
            "{unplaced}{} no coordinates and {} left out.",
            if unplaced == 1 { " atom site has" } else { " atom sites have" },
            if unplaced == 1 { "is" } else { "are" }
        ));
    }
    if empty > 0 {
        add(format!(
            "{empty}{} occupancy 0 and {} left out.",
            if empty == 1 { " atom site has" } else { " atom sites have" },
            if empty == 1 { "is" } else { "are" }
        ));
    }
    if with_atoms.len() > 1 {
        add(format!("The file holds {} structures; the first, '{}', is opened.", with_atoms.len(), block.name));
    }
    structure.load_note = note;
    Ok(true)
}

#[cfg(test)]
mod cartesian_tests {
    use super::*;

    #[test]
    fn cartesian_sites_match_fractional_ones() {
        let head = "data_x\n_cell_length_a 4\n_cell_length_b 5\n_cell_length_c 6\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 120\n";
        let fractional = format!("{head}loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nCu 0.25 0.5 0.1\nO 0 0 0.5\n");
        let (mut a, mut b) = (Structure::default(), Structure::default());
        assert!(read(&fractional, &mut a).unwrap());
        // The same sites in Cartesian coordinates of the standard frame.
        let cell = cell_from_parameters([4.0, 5.0, 6.0], [90.0, 90.0, 120.0]).unwrap();
        let p = cell * DVec3::new(0.25, 0.5, 0.1);
        let q = cell * DVec3::new(0.0, 0.0, 0.5);
        let cartesian = format!("{head}loop_\n_atom_site_label\n_atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\nCu {} {} {}\nO {} {} {}\n", p.x, p.y, p.z, q.x, q.y, q.z);
        assert!(read(&cartesian, &mut b).unwrap());
        assert_eq!(a.atoms.len(), b.atoms.len());
        for (x, y) in a.atoms.iter().zip(&b.atoms) {
            assert_eq!(x.atomic_number, y.atomic_number);
            assert!((x.x - y.x).abs() < 1e-9 && (x.y - y.y).abs() < 1e-9 && (x.z - y.z).abs() < 1e-9);
        }
    }
}
