//! The molecule library against the C++ one (tests/golden/molecules.txt,
//! from tools/golden/golden_molecules.cpp: name, element and the position as
//! C hex floats), bit for bit.

fn hex(x: f64) -> String {
    // C's %a: [-]0x1.<hex digits, trailing zeros trimmed>p<exponent>, zero as 0x0p+0.
    if x == 0.0 {
        return if x.is_sign_negative() { "-0x0p+0".into() } else { "0x0p+0".into() };
    }
    let b = x.to_bits();
    let sign = if b >> 63 == 1 { "-" } else { "" };
    let e = ((b >> 52) & 0x7ff) as i64 - 1023;
    let mut digits = format!("{:013x}", b & ((1u64 << 52) - 1));
    while digits.ends_with('0') {
        digits.pop();
    }
    let es = if e >= 0 { format!("+{e}") } else { format!("{e}") };
    if digits.is_empty() {
        format!("{sign}0x1p{es}")
    } else {
        format!("{sign}0x1.{digits}p{es}")
    }
}

#[cfg_attr(not(windows), ignore = "bit for bit against C++ references made on Windows; other maths libraries differ in the last digits")]
#[test]
fn molecule_library_equals_cpp() {
    let golden = std::fs::read_to_string(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/molecules.txt")).unwrap();
    let mut lines = Vec::new();
    for name in af_pipeline::molecule_library::molecule_names() {
        for a in af_pipeline::molecule_library::molecule(&name).unwrap() {
            lines.push(format!("{} {} {} {} {}", name, a.symbol, hex(a.position[0]), hex(a.position[1]), hex(a.position[2])));
        }
    }
    let expected: Vec<&str> = golden.lines().collect();
    assert_eq!(lines.len(), expected.len());
    for (a, b) in lines.iter().zip(&expected) {
        assert_eq!(a, b);
    }
}
