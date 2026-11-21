//! af_core::eigen::self_adjoint_eigen3 against Eigen's
//! SelfAdjointEigenSolver<Matrix3d> (tests/golden/self_adjoint_eigen.txt from
//! tools/golden/golden_selfadjoint.cpp: matrix | eigenvalues, eigenvectors).

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
fn self_adjoint_eigen_equals_eigen() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/self_adjoint_eigen.txt");
    let text = std::fs::read_to_string(path).unwrap();
    let mut differences = 0;
    for (n, line) in text.lines().enumerate() {
        let (input, output) = line.split_once('|').unwrap();
        let m: Vec<f64> = input.split_whitespace().map(parse_hex).collect();
        let expected: Vec<f64> = output.split_whitespace().map(parse_hex).collect();
        let matrix = [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]];
        let (values, vectors) = af_core::eigen::self_adjoint_eigen3(&matrix).expect("converges");
        let got: Vec<f64> = values.iter().chain(vectors.iter().flatten()).cloned().collect();
        if got.iter().zip(&expected).any(|(a, b)| a.to_bits() != b.to_bits()) {
            differences += 1;
            if differences <= 3 {
                println!("case {n}: {got:?}\nvs {expected:?}");
            }
        }
    }
    assert_eq!(differences, 0, "{differences} of {} decompositions differ", text.lines().count());
}
