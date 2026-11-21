//! af_core::crmath against GCC's constant-folded (MPFR, correctly rounded)
//! sin and cos: tests/golden/folded_trig.txt holds argument, sin, cos as C
//! hex floats, written by a C++ program whose 3000 calls GCC folded.

fn parse_hex(text: &str) -> f64 {
    let (negative, body) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let body = body.strip_prefix("0x").expect("hex float");
    let (mantissa, exponent) = body.split_once('p').expect("exponent");
    let (int_part, frac_part) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let mut value = u64::from_str_radix(int_part, 16).unwrap() as f64;
    let mut scale = 1.0 / 16.0;
    for c in frac_part.chars() {
        value += c.to_digit(16).unwrap() as f64 * scale;
        scale /= 16.0;
    }
    let v = value * 2f64.powi(exponent.parse::<i32>().unwrap());
    if negative {
        -v
    } else {
        v
    }
}

#[test]
fn correctly_rounded_trig_equals_gcc_constant_folding() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/folded_trig.txt");
    let text = std::fs::read_to_string(root).unwrap();
    let mut differences = Vec::new();
    let mut runtime_differences = 0;
    for line in text.lines() {
        let f: Vec<f64> = line.split_whitespace().map(parse_hex).collect();
        let (x, s, c) = (f[0], f[1], f[2]);
        if af_core::crmath::sin(x).to_bits() != s.to_bits() || af_core::crmath::cos(x).to_bits() != c.to_bits() {
            differences.push(format!("{x:e}: sin {:e} vs {s:e}, cos {:e} vs {c:e}", af_core::crmath::sin(x), af_core::crmath::cos(x)));
        }
        if x.sin().to_bits() != s.to_bits() || x.cos().to_bits() != c.to_bits() {
            runtime_differences += 1;
        }
    }
    println!("{runtime_differences} of {} arguments differ with the runtime libm", text.lines().count());
    assert!(differences.is_empty(), "{} differences:\n{}", differences.len(), differences.join("\n"));
}
