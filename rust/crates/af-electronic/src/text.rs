//! Text parsing as the C++ library does it (std::stod / std::stoi on a
//! prefix, istream words, text-mode lines).

use af_core::format::stod;

/// std::stod: the leading number; Err("stod") without one or out of range.
pub fn stod_prefix(s: &str) -> Result<f64, String> {
    stod(s).map(|(v, _)| v).map_err(|_| "stod".to_string())
}

/// std::stoi: the leading base-10 integer; Err("stoi") without one or out of range.
pub fn stoi_prefix(s: &str) -> Result<i32, String> {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && is_space(b[i]) {
        i += 1;
    }
    let start = i;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let digits = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == digits {
        return Err("stoi".into());
    }
    match s[start..i].parse::<i64>() {
        Ok(v) if v >= i64::from(i32::MIN) && v <= i64::from(i32::MAX) => Ok(v as i32),
        _ => Err("stoi".into()),
    }
}

/// C isspace.
pub fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\x0B' | b'\x0C' | b'\r')
}

/// istream >> std::string words.
pub fn words(s: &str) -> Vec<&str> {
    s.split(|c: char| c.is_ascii() && is_space(c as u8)).filter(|w| !w.is_empty()).collect()
}

/// The lines of a file read with std::getline (text mode), a final "\r" stripped.
pub fn lines_of(path: &str) -> Option<Vec<String>> {
    let bytes = std::fs::read(path).ok()?;
    let text = String::from_utf8_lossy(&bytes).replace("\r\n", "\n");
    let mut lines: Vec<String> = text.split('\n').map(|l| l.strip_suffix('\r').unwrap_or(l).to_string()).collect();
    if text.ends_with('\n') || text.is_empty() {
        lines.pop();
    }
    Some(lines)
}
