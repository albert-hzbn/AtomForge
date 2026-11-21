//! Number formatting and parsing that match the C++ streams and C library,
//! so text files written by the port are identical to the original's.

/// `printf("%.{precision}g", value)`: the shortest of fixed and scientific
/// notation with `precision` significant digits, trailing zeros removed; as a
/// C++ `ostream` with `precision(p)` prints doubles.
pub fn format_g(value: f64, precision: usize) -> String {
    let p = precision.max(1);
    if value.is_nan() {
        return if value.is_sign_negative() { "-nan".into() } else { "nan".into() };
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf".into() } else { "inf".into() };
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0".into() } else { "0".into() };
    }
    // The decimal exponent after rounding to p significant digits.
    let sci = format!("{:.*e}", p - 1, value);
    let (mantissa, exponent) = sci.split_once('e').expect("scientific form");
    let exponent: i32 = exponent.parse().expect("exponent");
    if exponent < -4 || exponent >= p as i32 {
        let mantissa = strip_zeros(mantissa);
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.abs())
    } else {
        let decimals = (p as i32 - 1 - exponent).max(0) as usize;
        strip_zeros(&format!("{:.*}", decimals, value)).to_string()
    }
}

fn strip_zeros(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
    }
}

/// `strtod`: the number at the start of `text` (after spaces), and the
/// number of bytes it used; `None` when no number starts there.
pub fn parse_prefix_f64(text: &str) -> Option<(f64, usize)> {
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() && (bytes[i] as char).is_ascii_whitespace() {
        i += 1;
    }
    let start = i;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    // Only the next few bytes can spell infinity/inf/nan (copying the whole
    // rest made parsing long texts quadratic).
    let lower = bytes[i..bytes.len().min(i + 8)].to_ascii_lowercase();
    for word in ["infinity", "inf", "nan"] {
        if lower.starts_with(word.as_bytes()) {
            let end = i + word.len();
            return text[start..end].parse::<f64>().ok().map(|v| (v, end));
        }
    }
    let digits_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    let mut digits = i - digits_start;
    if i < bytes.len() && bytes[i] == b'.' {
        i += 1;
        let fraction_start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        digits += i - fraction_start;
    }
    if digits == 0 {
        return None;
    }
    let mantissa_end = i;
    if i < bytes.len() && (bytes[i] == b'e' || bytes[i] == b'E') {
        let mut j = i + 1;
        if j < bytes.len() && (bytes[j] == b'+' || bytes[j] == b'-') {
            j += 1;
        }
        let exponent_digits = j;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
        if j > exponent_digits {
            i = j;
        }
    }
    let end = if i > mantissa_end { i } else { mantissa_end };
    text[start..end].parse::<f64>().ok().map(|v| (v, end))
}

/// A hexadecimal float ("0x1.8p3") at `i` (after the sign): value and end.
fn hex_float(bytes: &[u8], i: usize) -> Option<(f64, usize)> {
    if i + 2 > bytes.len() || bytes[i] != b'0' || (bytes[i + 1] | 0x20) != b'x' {
        return None;
    }
    let mut j = i + 2;
    let mut mantissa = 0u64;
    let mut exponent = 0i64;
    let mut digits = 0;
    let mut sticky = false;
    let mut add = |d: u64, fraction: bool, mantissa: &mut u64, exponent: &mut i64| {
        if *mantissa >> 56 == 0 {
            *mantissa = *mantissa << 4 | d;
            if fraction {
                *exponent -= 4;
            }
        } else {
            sticky |= d != 0;
            if !fraction {
                *exponent += 4;
            }
        }
    };
    while j < bytes.len() && bytes[j].is_ascii_hexdigit() {
        add((bytes[j] as char).to_digit(16).unwrap() as u64, false, &mut mantissa, &mut exponent);
        digits += 1;
        j += 1;
    }
    if j < bytes.len() && bytes[j] == b'.' {
        j += 1;
        while j < bytes.len() && bytes[j].is_ascii_hexdigit() {
            add((bytes[j] as char).to_digit(16).unwrap() as u64, true, &mut mantissa, &mut exponent);
            digits += 1;
            j += 1;
        }
    }
    if digits == 0 {
        return None;
    }
    if j < bytes.len() && (bytes[j] | 0x20) == b'p' {
        let mut k = j + 1;
        let negative = k < bytes.len() && bytes[k] == b'-';
        if k < bytes.len() && (bytes[k] == b'+' || bytes[k] == b'-') {
            k += 1;
        }
        let start = k;
        let mut e = 0i64;
        while k < bytes.len() && bytes[k].is_ascii_digit() {
            e = (e * 10 + (bytes[k] - b'0') as i64).min(1 << 40);
            k += 1;
        }
        if k > start {
            exponent += if negative { -e } else { e };
            j = k;
        }
    }
    // The low bit records discarded digits so the one rounding is correct.
    let m = mantissa | sticky as u64;
    let value = (m as f64) * 2f64.powi(exponent.clamp(-1100, 1100) as i32);
    Some((value, j))
}

/// C `strtod`: the number at the start of `text` (after spaces; decimal,
/// hexadecimal, inf, nan) and the bytes used, (0, 0) when none.
pub fn strtod(text: &str) -> (f64, usize) {
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() && (bytes[i] as char).is_ascii_whitespace() {
        i += 1;
    }
    let negative = i < bytes.len() && bytes[i] == b'-';
    let after_sign = if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') { i + 1 } else { i };
    if let Some((value, end)) = hex_float(bytes, after_sign) {
        return (if negative { -value } else { value }, end);
    }
    parse_prefix_f64(text).unwrap_or((0.0, 0))
}

/// Why `std::stod` throws.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StodError {
    /// No number at the start (std::invalid_argument).
    Invalid,
    /// Overflow or underflow to zero (std::out_of_range).
    Range,
}

/// `std::stod`: leading spaces skipped, the longest number prefix converted
/// (decimal, hexadecimal, inf, nan) with the bytes used.
pub fn stod(text: &str) -> Result<(f64, usize), StodError> {
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() && (bytes[i] as char).is_ascii_whitespace() {
        i += 1;
    }
    let negative = i < bytes.len() && bytes[i] == b'-';
    let after_sign = if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') { i + 1 } else { i };
    if let Some((value, end)) = hex_float(bytes, after_sign) {
        let value = if negative { -value } else { value };
        if value.is_infinite() {
            return Err(StodError::Range);
        }
        return Ok((value, end));
    }
    let (value, used) = parse_prefix_f64(text).ok_or(StodError::Invalid)?;
    let written = &bytes[after_sign..used];
    let word = written.first().is_some_and(|c| c.is_ascii_alphabetic());
    if !word {
        let nonzero = written.iter().take_while(|&&c| c != b'e' && c != b'E').any(|c| (b'1'..=b'9').contains(c));
        if value.is_infinite() || (value == 0.0 && nonzero) {
            return Err(StodError::Range);
        }
    }
    Ok((value, used))
}

/// `strtof`: as [`parse_prefix_f64`], rounded once to f32.
pub fn parse_prefix_f32(text: &str) -> Option<(f32, usize)> {
    let (_, end) = parse_prefix_f64(text)?;
    let start = text.len() - text.trim_start().len();
    text[start..end].parse::<f32>().ok().map(|v| (v, end))
}

/// Whitespace-separated numbers read as `std::istream >> float` does: each
/// word must start with a number.
pub fn read_f32_words(text: &str, count: usize) -> Option<Vec<f32>> {
    let mut out = Vec::with_capacity(count);
    let mut rest = text;
    for _ in 0..count {
        let (value, used) = parse_prefix_f32(rest)?;
        out.push(value);
        rest = &rest[used..];
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn printf_g_examples() {
        assert_eq!(format_g(0.4f32 as f64, 9), "0.400000006");
        assert_eq!(format_g(1.25, 9), "1.25");
        assert_eq!(format_g(1e-5, 6), "1e-05");
        assert_eq!(format_g(123456789.0, 6), "1.23457e+08");
        assert_eq!(format_g(100000.0, 6), "100000");
        assert_eq!(format_g(1000000.0, 6), "1e+06");
    }

    #[test]
    fn strtod_prefixes() {
        assert_eq!(parse_prefix_f64("  1.5abc"), Some((1.5, 5)));
        assert_eq!(parse_prefix_f64("-2e3x"), Some((-2000.0, 4)));
        assert_eq!(parse_prefix_f64("3e"), Some((3.0, 1)));
        assert_eq!(parse_prefix_f64(".5"), Some((0.5, 2)));
        assert_eq!(parse_prefix_f64("abc"), None);
        assert_eq!(parse_prefix_f64("."), None);
    }
}
