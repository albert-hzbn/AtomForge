//! Correctly rounded sin and cos, as GCC computes them (with MPFR) when it
//! folds a call on a compile-time constant: `std::sin(2 * pi * k / 3)` in an
//! unrolled loop gives the correctly rounded value, which the runtime libm
//! may miss by an ulp. Use these where the C++ argument is a constant.
//!
//! Double-double arithmetic: the argument reduced by a three-part pi/2, then
//! Taylor series to about 1e-32, rounded once at the end.

#[derive(Clone, Copy, Debug)]
struct Dd {
    hi: f64,
    lo: f64,
}

fn two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    Dd { hi: s, lo: err }
}

fn quick_two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    Dd { hi: s, lo: b - (s - a) }
}

fn two_prod(a: f64, b: f64) -> Dd {
    let p = a * b;
    Dd { hi: p, lo: a.mul_add(b, -p) }
}

fn add(a: Dd, b: Dd) -> Dd {
    let s = two_sum(a.hi, b.hi);
    let t = two_sum(a.lo, b.lo);
    let s = quick_two_sum(s.hi, s.lo + t.hi);
    quick_two_sum(s.hi, s.lo + t.lo)
}

fn neg(a: Dd) -> Dd {
    Dd { hi: -a.hi, lo: -a.lo }
}

fn mul(a: Dd, b: Dd) -> Dd {
    let p = two_prod(a.hi, b.hi);
    quick_two_sum(p.hi, p.lo + (a.hi * b.lo + a.lo * b.hi))
}

fn div_f(a: Dd, d: f64) -> Dd {
    let q1 = a.hi / d;
    let p = two_prod(q1, d);
    let r = two_sum(a.hi, -p.hi);
    let q2 = (r.hi + (r.lo - p.lo + a.lo)) / d;
    quick_two_sum(q1, q2)
}

const PIO2: [f64; 3] = [1.5707963267948966, 6.123233995736766e-17, -1.4973849048591698e-33];

/// x reduced to [-pi/4, pi/4] as a double-double, and the quadrant.
fn reduce(x: f64) -> (Dd, i64) {
    let k = (x / PIO2[0]).round();
    let mut r = Dd { hi: x, lo: 0.0 };
    for part in PIO2 {
        r = add(r, neg(two_prod(k, part)));
    }
    (r, k as i64)
}

/// sin and cos of a reduced argument (|r| <= pi/4).
fn kernel(r: Dd) -> (Dd, Dd) {
    let r2 = mul(r, r);
    // sin: r - r^3/3! + ...
    let mut term = r;
    let mut s = r;
    let mut n = 1.0;
    for _ in 0..16 {
        term = div_f(mul(term, neg(r2)), (n + 1.0) * (n + 2.0));
        n += 2.0;
        s = add(s, term);
    }
    // cos: 1 - r^2/2! + ...
    let mut term = Dd { hi: 1.0, lo: 0.0 };
    let mut c = term;
    let mut n = 0.0;
    for _ in 0..16 {
        term = div_f(mul(term, neg(r2)), (n + 1.0) * (n + 2.0));
        n += 2.0;
        c = add(c, term);
    }
    (s, c)
}

fn round(d: Dd) -> f64 {
    d.hi + d.lo
}

/// Correctly rounded sin (|x| up to about 1e6).
pub fn sin(x: f64) -> f64 {
    if !x.is_finite() || x == 0.0 {
        return x.sin();
    }
    let (r, k) = reduce(x);
    let (s, c) = kernel(r);
    round(match k.rem_euclid(4) {
        0 => s,
        1 => c,
        2 => neg(s),
        _ => neg(c),
    })
}

/// Correctly rounded cos (|x| up to about 1e6).
pub fn cos(x: f64) -> f64 {
    if !x.is_finite() {
        return x.cos();
    }
    let (r, k) = reduce(x);
    let (s, c) = kernel(r);
    round(match k.rem_euclid(4) {
        0 => c,
        1 => neg(s),
        2 => neg(c),
        _ => s,
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn known_values() {
        let pi = std::f64::consts::PI;
        assert_eq!(super::cos(0.0), 1.0);
        assert_eq!(super::sin(pi / 6.0), 0.49999999999999994);
        assert_eq!(super::cos(pi / 3.0), 0.5000000000000001);
        for k in 0..2000 {
            let x = -7.0 + k as f64 * 0.00731;
            assert!((super::sin(x) - x.sin()).abs() <= 2.3e-16, "{x}");
            assert!((super::cos(x) - x.cos()).abs() <= 2.3e-16, "{x}");
        }
    }
}
