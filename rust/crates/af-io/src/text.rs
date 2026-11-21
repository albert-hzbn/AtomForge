//! Reading text as the C++ readers do: std::getline lines, istream words, and
//! the strtod / std::stod / std::stoll conversions with their prefix rules.

use af_core::format::parse_prefix_f64;

/// Lines as std::getline returns them (a trailing newline ends the last line;
/// no empty line follows it), with a trailing '\r' removed when `strip_cr`.
pub fn lines(text: &str, strip_cr: bool) -> Vec<&str> {
    let mut out = Vec::new();
    let mut rest = text;
    while !rest.is_empty() {
        let (line, next) = match rest.find('\n') {
            Some(i) => (&rest[..i], &rest[i + 1..]),
            None => (rest, ""),
        };
        out.push(if strip_cr { line.strip_suffix('\r').unwrap_or(line) } else { line });
        rest = next;
    }
    out
}

/// isspace in the C locale.
pub fn is_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')
}

/// Words as `istream >> std::string` reads them.
pub fn words(line: &str) -> Vec<&str> {
    line.split(is_space).filter(|w| !w.is_empty()).collect()
}

/// strtod over the whole word (nothing may follow), finite: the readers'
/// `number` helpers.
pub fn number_full(word: &str) -> Option<f64> {
    let (v, used) = parse_prefix_f64(word)?;
    (used == word.len() && v.is_finite()).then_some(v)
}

/// std::stod: leading spaces skipped, a number prefix converted; an error
/// when none (or out of range).
pub fn stod(text: &str) -> Result<f64, String> {
    let (v, _) = parse_prefix_f64(text).ok_or_else(|| format!("stod: no conversion of '{text}'"))?;
    // An overflow (ERANGE) throws out_of_range; a written "inf" or "nan" is a value.
    let written = text.trim_start().trim_start_matches(['+', '-']).as_bytes();
    if v.is_infinite() && !written.get(..3).is_some_and(|w| w.eq_ignore_ascii_case(b"inf")) {
        return Err(format!("stod: '{text}' is out of range"));
    }
    Ok(v)
}

/// std::stod and the number of characters used.
pub fn stod_used(text: &str) -> Result<(f64, usize), String> {
    parse_prefix_f64(text).ok_or_else(|| format!("stod: no conversion of '{text}'"))
}

/// strtol / strtoll in base 10: the value and the byte offset after it (0
/// when no digits; the value is then 0).
pub fn strtol(text: &str) -> (i64, usize) {
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() && is_space(bytes[i] as char) {
        i += 1;
    }
    let mut negative = false;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        negative = bytes[i] == b'-';
        i += 1;
    }
    let start = i;
    let mut value: i64 = 0;
    let mut overflow = false;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        match value.checked_mul(10).and_then(|v| v.checked_add((bytes[i] - b'0') as i64)) {
            Some(v) => value = v,
            None => overflow = true,
        }
        i += 1;
    }
    if i == start {
        return (0, 0);
    }
    if overflow {
        value = i64::MAX;
    }
    (if negative { -value } else { value }, i)
}

/// std::stoll: a prefix integer, an error when none.
pub fn stoll(text: &str) -> Result<i64, String> {
    let (v, used) = strtol(text);
    if used == 0 {
        Err(format!("stoll: no conversion of '{text}'"))
    } else {
        Ok(v)
    }
}

/// std::stoi: as stoll, and an error outside int.
pub fn stoi(text: &str) -> Result<i32, String> {
    let v = stoll(text)?;
    i32::try_from(v).map_err(|_| format!("stoi: '{text}' is out of range"))
}

/// std::atoi.
pub fn atoi(text: &str) -> i32 {
    strtol(text).0 as i32
}

pub fn lower(s: &str) -> String {
    s.to_ascii_lowercase()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn getline_lines() {
        assert_eq!(lines("a\nb\n", true), vec!["a", "b"]);
        assert_eq!(lines("a\r\n\nb", true), vec!["a", "", "b"]);
        assert!(lines("", true).is_empty());
    }

    #[test]
    fn conversions() {
        assert_eq!(stod(" 2.5x").unwrap(), 2.5);
        assert!(stod("x").is_err());
        assert_eq!(stoll("12abc").unwrap(), 12);
        assert_eq!(strtol("-7 rest"), (-7, 2));
        assert_eq!(number_full("1e3"), Some(1000.0));
        assert_eq!(number_full("1e3x"), None);
    }
}
