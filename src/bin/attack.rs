//! Demo attack that recovers s(1) using CRT projections and trivial solvers.
//!
//! The implementation mirrors the reference pseudo attack provided by the user
//! but reuses the crate's core types such as [`Polynomial`], [`Cipher`] and
//! [`ArithChannel`].

use aces_core::{ArithChannel, Cipher, Polynomial};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::HashMap;

// ---- Parameters ----
const N: usize = 16;
const Q: u128 = (1u128 << 60) + 1; // 16^15 + 1
const Q_FACTORS: [u128; 4] = [17, 241, 61_681, 4_562_284_561];
const M_SAMPLES: usize = 400;
const IDEAL_BOUND: i64 = 2; // e = (X-1)r の r の係数幅

// ---- modular helpers ----
#[inline]
fn mod_sub(a: u128, b: u128, m: u128) -> u128 {
    if a >= b {
        a - b
    } else {
        m - (b - a)
    }
}

#[inline]
fn mod_mul(a: u128, b: u128, m: u128) -> u128 {
    ((a as u128).wrapping_mul(b as u128)) % m
}

fn mod_inverse(a: u128, m: u128) -> Option<u128> {
    use aces_core::arith_channel::extended_gcd2;
    let (g, x, _) = extended_gcd2(a % m, m);
    if g != 1 {
        return None;
    }
    let mut x = x % m as i128;
    if x < 0 {
        x += m as i128;
    }
    Some(x as u128)
}

fn crt_pair(x1: u128, m1: u128, x2: u128, m2: u128) -> Option<(u128, u128)> {
    use aces_core::arith_channel::extended_gcd2;
    let (g, s, t) = extended_gcd2(m1, m2);
    if g != 1 {
        return None;
    }
    let modulus = m1.checked_mul(m2)?;
    let s_u = ((s % m2 as i128) + m2 as i128) as u128 % m2;
    let t_u = ((t % m1 as i128) + m1 as i128) as u128 % m1;
    let term1 = mod_mul(mod_mul(x1 % m1, t_u, m1), m2, modulus);
    let term2 = mod_mul(mod_mul(x2 % m2, s_u, m2), m1, modulus);
    Some(((term1 + term2) % modulus, modulus))
}

fn crt_all(residues: &[u128], moduli: &[u128]) -> Option<u128> {
    assert_eq!(residues.len(), moduli.len());
    let mut x = residues[0];
    let mut m = moduli[0];
    for (&r, &modulus) in residues.iter().zip(moduli).skip(1) {
        let (nx, nm) = crt_pair(x, m, r, modulus)?;
        x = nx;
        m = nm;
    }
    Some(x)
}

// ---- Polynomial helpers ----
fn random_poly(rng: &mut StdRng) -> Polynomial {
    let coeffs: Vec<u128> = (0..N).map(|_| rng.gen_range(0..Q)).collect();
    Polynomial::new(coeffs, Q)
}

fn reduce(poly: Polynomial, modulus_poly: &Polynomial) -> Polynomial {
    &poly % modulus_poly
}

fn poly_add_mod(a: &Polynomial, b: &Polynomial, modulus_poly: &Polynomial) -> Polynomial {
    let sum = a + b;
    reduce(sum, modulus_poly)
}

fn poly_mul_mod(a: &Polynomial, b: &Polynomial, modulus_poly: &Polynomial) -> Polynomial {
    let product = a * b;
    reduce(product, modulus_poly)
}

fn random_poly_reduced(rng: &mut StdRng, modulus_poly: &Polynomial) -> Polynomial {
    let poly = random_poly(rng);
    reduce(poly, modulus_poly)
}

fn e_ideal_xminus1(rng: &mut StdRng, bound: i64, modulus_poly: &Polynomial) -> Polynomial {
    // Generate r with coefficients in [-bound, bound]
    let mut r_coeffs = Vec::with_capacity(N);
    for _ in 0..N {
        let x: i64 = rng.gen_range(-bound..=bound);
        let xi = ((x % Q as i64) + Q as i64) as u128 % Q;
        r_coeffs.push(xi);
    }
    let r_poly = Polynomial::new(r_coeffs, Q);
    let x_poly = Polynomial::new(vec![0, 1], Q);
    let xr = poly_mul_mod(&x_poly, &r_poly, modulus_poly);
    let diff = &xr - &r_poly;
    reduce(diff, modulus_poly)
}

fn solve_trivial_majority(pairs: &[(u128, u128)], modulus: u128) -> Option<u128> {
    let mut counts: HashMap<u128, usize> = HashMap::new();
    for &(a, b) in pairs {
        let a = a % modulus;
        if a == 0 {
            continue;
        }
        if let Some(inv) = mod_inverse(a, modulus) {
            let g = mod_mul(b % modulus, inv, modulus);
            *counts.entry(g).or_default() += 1;
        }
    }
    counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(value, _)| value)
}

fn main() {
    let mut rng = StdRng::seed_from_u64(0xC0FFEE);

    // Arithmetic channel configured for our parameters.
    let mut chan = ArithChannel::new(2, Q, N, Q_FACTORS.len());
    // Override modulus polynomial to enforce X^N - 1 structure.
    let mut modulus_coeffs = vec![0u128; N + 1];
    modulus_coeffs[0] = mod_sub(0, 1, Q); // -1 mod Q => Q - 1
    modulus_coeffs[N] = 1;
    chan.u = Polynomial::new(modulus_coeffs, Q);

    // Secret & samples
    let s = random_poly_reduced(&mut rng, &chan.u);
    let s1 = chan.eval(&s);

    let mut samples: Vec<Cipher> = Vec::with_capacity(M_SAMPLES);
    for _ in 0..M_SAMPLES {
        let a = random_poly_reduced(&mut rng, &chan.u);
        let e = e_ideal_xminus1(&mut rng, IDEAL_BOUND, &chan.u); // ★ e(1)=0
        let prod = poly_mul_mod(&a, &s, &chan.u);
        let b = poly_add_mod(&prod, &e, &chan.u);
        samples.push(Cipher {
            dec: vec![a],
            enc: b,
            level: 0,
        });
    }

    println!(
        "Ring: Z_{}[X]/(X^{}) with X^{}-1; factors = {:?}",
        Q, N, N, Q_FACTORS
    );
    println!("Mode: ideal error e=(X-1)r → e(1)=0");
    println!("True s(1) mod q = {}", s1);

    // Projection & trivial solve at each q_i
    let mut residues = Vec::new();
    for &qi in &Q_FACTORS {
        let mut pairs = Vec::with_capacity(samples.len());
        for sample in &samples {
            let ap = chan.eval(&sample.dec[0]) % qi;
            let bp = chan.eval(&sample.enc) % qi;
            pairs.push((ap, bp));
        }
        match solve_trivial_majority(&pairs, qi) {
            Some(g) => {
                println!("qi = {:>11}: recovered s(1) ≡ {} (mod {})", qi, g, qi);
                residues.push(g);
            }
            None => {
                println!("qi = {:>11}: could NOT recover", qi);
                return;
            }
        }
    }

    let s1_rec = crt_all(&residues, &Q_FACTORS).expect("CRT failed");
    println!("CRT recombined s(1) mod q = {}", s1_rec);
    println!(
        "Match = {}",
        if s1_rec % Q == s1 % Q {
            "SUCCESS"
        } else {
            "FAIL"
        }
    );
}
