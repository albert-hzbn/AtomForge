//! The C++ random-number machinery the C++ code uses, reproduced exactly so
//! seeded results match the C++ build (GCC's libstdc++): std::mt19937 and
//! std::mt19937_64, generate_canonical, uniform_real_distribution,
//! uniform_int_distribution (Lemire's nearly divisionless method),
//! normal_distribution (Marsaglia's polar method) and std::shuffle.

/// A uniform random bit generator with values in [0, max].
pub trait Engine {
    fn next_u64(&mut self) -> u64;
    fn max(&self) -> u64;
}

/// std::mt19937 (32-bit Mersenne Twister).
#[derive(Clone)]
pub struct Mt19937 {
    state: [u32; 624],
    index: usize,
}

impl Mt19937 {
    pub fn new(seed: u32) -> Self {
        let mut state = [0u32; 624];
        state[0] = seed;
        for i in 1..624 {
            state[i] = 1812433253u32.wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30)).wrapping_add(i as u32);
        }
        Self { state, index: 624 }
    }

    fn twist(&mut self) {
        for i in 0..624 {
            let y = (self.state[i] & 0x8000_0000) | (self.state[(i + 1) % 624] & 0x7fff_ffff);
            let mut v = self.state[(i + 397) % 624] ^ (y >> 1);
            if y & 1 != 0 {
                v ^= 0x9908_b0df;
            }
            self.state[i] = v;
        }
        self.index = 0;
    }

    pub fn next_u32(&mut self) -> u32 {
        if self.index >= 624 {
            self.twist();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// discard(n)
    pub fn discard(&mut self, n: u64) {
        for _ in 0..n {
            self.next_u32();
        }
    }
}

impl Engine for Mt19937 {
    fn next_u64(&mut self) -> u64 {
        self.next_u32() as u64
    }
    fn max(&self) -> u64 {
        u32::MAX as u64
    }
}

/// std::mt19937_64.
#[derive(Clone)]
pub struct Mt19937_64 {
    state: [u64; 312],
    index: usize,
}

impl Mt19937_64 {
    pub fn new(seed: u64) -> Self {
        let mut state = [0u64; 312];
        state[0] = seed;
        for i in 1..312 {
            state[i] = 6364136223846793005u64.wrapping_mul(state[i - 1] ^ (state[i - 1] >> 62)).wrapping_add(i as u64);
        }
        Self { state, index: 312 }
    }

    fn twist(&mut self) {
        const UPPER: u64 = 0xffff_ffff_8000_0000;
        const LOWER: u64 = 0x7fff_ffff;
        for i in 0..312 {
            let y = (self.state[i] & UPPER) | (self.state[(i + 1) % 312] & LOWER);
            let mut v = self.state[(i + 156) % 312] ^ (y >> 1);
            if y & 1 != 0 {
                v ^= 0xb502_6f5a_a966_19e9;
            }
            self.state[i] = v;
        }
        self.index = 0;
    }

    pub fn next(&mut self) -> u64 {
        if self.index >= 312 {
            self.twist();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= (y >> 29) & 0x5555_5555_5555_5555;
        y ^= (y << 17) & 0x71d6_7fff_eda6_0000;
        y ^= (y << 37) & 0xfff7_eee0_0000_0000;
        y ^ (y >> 43)
    }
}

impl Engine for Mt19937_64 {
    fn next_u64(&mut self) -> u64 {
        self.next()
    }
    fn max(&self) -> u64 {
        u64::MAX
    }
}

/// std::generate_canonical<double, 53>: one 64-bit draw, or two 32-bit draws.
pub fn generate_canonical<E: Engine>(engine: &mut E) -> f64 {
    let ret = if engine.max() == u64::MAX {
        let tmp = 18446744073709551616.0f64; // 2^64
        engine.next_u64() as f64 / tmp
    } else {
        // r = 2^32, m = 2 terms.
        let r = 4294967296.0f64;
        let mut sum = 0.0f64;
        let mut tmp = 1.0f64;
        for _ in 0..2 {
            sum += engine.next_u64() as f64 * tmp;
            tmp *= r;
        }
        sum / tmp
    };
    if ret >= 1.0 {
        1.0f64.next_down()
    } else {
        ret
    }
}

/// std::uniform_real_distribution<double>(a, b).
pub fn uniform_real<E: Engine>(engine: &mut E, a: f64, b: f64) -> f64 {
    generate_canonical(engine) * (b - a) + a
}

/// std::generate_canonical<float, 24>: one draw (of 32 or 64 bits) as a
/// float, divided by the engine's range.
pub fn generate_canonical_f32<E: Engine>(engine: &mut E) -> f32 {
    let range = if engine.max() == u64::MAX { 18446744073709551616.0f32 } else { 4294967296.0f32 };
    let ret = engine.next_u64() as f32 / range;
    if ret >= 1.0 {
        1.0f32.next_down()
    } else {
        ret
    }
}

/// std::uniform_real_distribution<float>(a, b).
pub fn uniform_real_f32<E: Engine>(engine: &mut E, a: f32, b: f32) -> f32 {
    generate_canonical_f32(engine) * (b - a) + a
}

/// Lemire's nearly divisionless reduction (libstdc++ _S_nd).
fn lemire64<E: Engine>(engine: &mut E, range: u64) -> u64 {
    let mut product = engine.next_u64() as u128 * range as u128;
    let mut low = product as u64;
    if low < range {
        let threshold = range.wrapping_neg() % range;
        while low < threshold {
            product = engine.next_u64() as u128 * range as u128;
            low = product as u64;
        }
    }
    (product >> 64) as u64
}

fn lemire32<E: Engine>(engine: &mut E, range: u32) -> u64 {
    let mut product = engine.next_u64() * range as u64;
    let mut low = product as u32;
    if low < range {
        let threshold = range.wrapping_neg() % range;
        while low < threshold {
            product = engine.next_u64() * range as u64;
            low = product as u32;
        }
    }
    product >> 32
}

/// std::uniform_int_distribution over [a, b] as unsigned offsets (the
/// libstdc++ algorithm; b - a must fit in u64).
pub fn uniform_u64<E: Engine>(engine: &mut E, a: u64, b: u64) -> u64 {
    let urngrange = engine.max();
    let urange = b.wrapping_sub(a);
    let ret = if urngrange > urange {
        let uerange = urange + 1;
        if urngrange == u64::MAX {
            lemire64(engine, uerange)
        } else if urngrange == u32::MAX as u64 {
            lemire32(engine, uerange as u32)
        } else {
            let scaling = urngrange / uerange;
            let past = uerange * scaling;
            loop {
                let r = engine.next_u64();
                if r < past {
                    break r / scaling;
                }
            }
        }
    } else if urngrange < urange {
        // Upscaling (a 32-bit engine for a 64-bit range).
        let uerngrange = urngrange + 1;
        loop {
            let high = uniform_u64(engine, 0, urange / uerngrange);
            let tmp = uerngrange.wrapping_mul(high);
            let r = tmp.wrapping_add(engine.next_u64());
            if r <= urange && r >= tmp {
                break r;
            }
        }
    } else {
        engine.next_u64()
    };
    ret.wrapping_add(a)
}

/// std::uniform_int_distribution<long long>(a, b).
pub fn uniform_i64<E: Engine>(engine: &mut E, a: i64, b: i64) -> i64 {
    // libstdc++ works on the unsigned type: (b - a) and the offset from a.
    let offset = uniform_u64(engine, 0, (b as u64).wrapping_sub(a as u64));
    (a as u64).wrapping_add(offset) as i64
}

/// std::normal_distribution<double>: pairs from Marsaglia's polar method,
/// the second value kept for the next call.
#[derive(Clone, Debug)]
pub struct Normal {
    pub mean: f64,
    pub stddev: f64,
    saved: Option<f64>,
}

impl Normal {
    pub fn new(mean: f64, stddev: f64) -> Self {
        Self { mean, stddev, saved: None }
    }

    pub fn sample<E: Engine>(&mut self, engine: &mut E) -> f64 {
        let value = if let Some(v) = self.saved.take() {
            v
        } else {
            let (x, y, r2) = loop {
                let x = 2.0 * generate_canonical(engine) - 1.0;
                let y = 2.0 * generate_canonical(engine) - 1.0;
                let r2 = x * x + y * y;
                if !(r2 > 1.0 || r2 == 0.0) {
                    break (x, y, r2);
                }
            };
            let mult = (-2.0 * r2.ln() / r2).sqrt();
            self.saved = Some(x * mult);
            y * mult
        };
        value * self.stddev + self.mean
    }
}

/// std::shuffle (libstdc++): two swap positions per draw when the engine's
/// range allows it.
pub fn shuffle<T, E: Engine>(v: &mut [T], engine: &mut E) {
    let n = v.len() as u64;
    if n == 0 {
        return;
    }
    let urngrange = engine.max();
    if urngrange / n >= n {
        let mut i = 1u64;
        if n % 2 == 0 {
            let j = uniform_u64(engine, 0, 1);
            v.swap(i as usize, j as usize);
            i += 1;
        }
        while i != n {
            let swap_range = i + 1;
            let x = uniform_u64(engine, 0, swap_range * (swap_range + 1) - 1);
            let (first, second) = (x / (swap_range + 1), x % (swap_range + 1));
            v.swap(i as usize, first as usize);
            i += 1;
            v.swap(i as usize, second as usize);
            i += 1;
        }
        return;
    }
    for i in 1..n {
        let j = uniform_u64(engine, 0, i);
        v.swap(i as usize, j as usize);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_mersenne_values() {
        // The C++ standard: the 10000th value of a default-seeded engine.
        let mut a = Mt19937::new(5489);
        a.discard(9999);
        assert_eq!(a.next_u32(), 4123659995);
        let mut b = Mt19937_64::new(5489);
        for _ in 0..9999 {
            b.next();
        }
        assert_eq!(b.next(), 9981545732273789042);
    }
}
