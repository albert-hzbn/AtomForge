//! Element symbols, names, covalent radii, masses and colours (the C++
//! util/ElementData and util/ElementColourSchemes).

use crate::element_tables as t;

/// The chemical symbol for atomic number z (1-118), or "?".
pub fn symbol(z: i32) -> &'static str {
    if (1..=118).contains(&z) {
        t::SYMBOLS[z as usize]
    } else {
        "?"
    }
}

/// The element name for z (1-118), or "Unknown".
pub fn name(z: i32) -> &'static str {
    if (1..=118).contains(&z) {
        t::NAMES[z as usize]
    } else {
        "Unknown"
    }
}

/// Atomic number of a symbol (exact case, as the C++ lookups), or 0.
pub fn atomic_number(symbol: &str) -> i32 {
    (1..=118).find(|&z| t::SYMBOLS[z as usize] == symbol).map_or(0, |z| z as i32)
}

/// Atomic number of a symbol in any case ("cu", "CU" -> 29), or 0.
pub fn atomic_number_any_case(symbol: &str) -> i32 {
    let mut chars = symbol.trim().chars();
    let Some(first) = chars.next() else { return 0 };
    let normal: String = std::iter::once(first.to_ascii_uppercase()).chain(chars.map(|c| c.to_ascii_lowercase())).collect();
    atomic_number(&normal)
}

/// Covalent radii (Angstrom; Cordero et al. 2008) indexed by z (0-118).
pub fn literature_covalent_radii() -> Vec<f32> {
    t::COVALENT_RADII.to_vec()
}

/// Covalent radius of z (1.6 outside the table, as the C++ default).
pub fn covalent_radius(z: i32) -> f32 {
    if (0..=118).contains(&z) {
        t::COVALENT_RADII[z as usize]
    } else {
        1.6
    }
}

/// Standard atomic weight of z; 1.0 out of range.
pub fn atomic_mass(z: i32) -> f64 {
    if (1..=118).contains(&z) {
        t::MASSES[z as usize]
    } else {
        1.0
    }
}

/// Default (CPK-style) colour; z is clamped to 1-118.
pub fn default_colour(z: i32) -> [f32; 3] {
    t::DEFAULT_COLOURS[z.clamp(1, 118) as usize]
}

/// Default colours indexed by z (0-118; index 0 is grey, as the C++ table
/// makeDefaultElementColors, unlike default_colour(0) which clamps to H).
pub fn default_colours() -> Vec<[f32; 3]> {
    std::iter::once([0.7; 3]).chain((1..=118).map(default_colour)).collect()
}

/// A named element colour scheme.
#[derive(Clone, Copy, Debug)]
pub struct ColourScheme {
    pub name: &'static str,
    pub description: &'static str,
}

/// AtomForge, Jmol, RasMol CPK and Greyscale.
pub fn colour_schemes() -> Vec<ColourScheme> {
    t::SCHEME_NAMES.iter().zip(t::SCHEME_DESCRIPTIONS.iter()).map(|(&name, &description)| ColourScheme { name, description }).collect()
}

/// Colours of every element (0-118) in a scheme; the defaults for an unknown scheme.
pub fn scheme_colours(scheme: i32) -> Vec<[f32; 3]> {
    if scheme >= 0 && (scheme as usize) < t::SCHEME_COLOURS.len() {
        t::SCHEME_COLOURS[scheme as usize].to_vec()
    } else {
        default_colours()
    }
}

/// The scheme a colour table equals, or -1.
pub fn matching_scheme(colours: &[[f32; 3]]) -> i32 {
    (0..t::SCHEME_COLOURS.len()).find(|&i| t::SCHEME_COLOURS[i][..] == *colours).map_or(-1, |i| i as i32)
}
