//! Attack attempting to recover s(1) using the actual ACES parameters.
//!
//! The code instantiates the scheme with the same parameter generation
//! pipeline used by the library (real `f₀`, `f₁` generation), collects a
//! handful of encryptions and solves the resulting linear systems modulo
//! each prime-power factor of `q`.
//!
//! Besides the core key-recovery logic, we dump a few useful diagnostics:
//!   * whether the vanishing noise `e` added to each `f₁` row satisfies
//!     `e(1) = 0` (spoiler: it equals `p`, not 0),
//!   * the gcd pattern of the public key evaluations — every column is a
//!     multiple of the same prime ideal factor, and
//!   * per-factor success when recombining the recovered secrets via CRT.

use aces_core::{Aces, ArithChannel};

const P: u128 = 32;
const Q: u128 = (1u128 << 40) + 1; // 2^40 + 1
const DIM: usize = 5;
const N: usize = 5;

fn mod_add(a: u128, b: u128, modulus: u128) -> u128 {
    let sum = (a % modulus) + (b % modulus);
    if sum >= modulus {
        sum - modulus
    } else {
        sum
    }
}

fn mod_sub(a: u128, b: u128, modulus: u128) -> u128 {
    if a >= b {
        (a - b) % modulus
    } else {
        (modulus - ((b - a) % modulus)) % modulus
    }
}

fn mod_mul(a: u128, b: u128, modulus: u128) -> u128 {
    ((a % modulus) * (b % modulus)) % modulus
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

fn mod_inverse(a: u128, modulus: u128) -> Option<u128> {
    use aces_core::arith_channel::extended_gcd2;
    let (g, x, _) = extended_gcd2(a % modulus, modulus);
    if g != 1 {
        return None;
    }
    let mut x = x % modulus as i128;
    if x < 0 {
        x += modulus as i128;
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
    if residues.is_empty() {
        return Some(0);
    }
    let mut x = residues[0] % moduli[0];
    let mut m = moduli[0];
    for (&residue, &modulus) in residues.iter().zip(moduli).skip(1) {
        let (nx, nm) = crt_pair(x, m, residue, modulus)?;
        x = nx % nm;
        m = nm;
    }
    Some(x % m)
}

fn factor_into_coprime(mut n: u128) -> Vec<u128> {
    let mut factors = Vec::new();
    let mut d: u128 = 2;
    while d * d <= n {
        let mut count = 0u32;
        while n % d == 0 {
            n /= d;
            count += 1;
        }
        if count > 0 {
            let mut power = 1u128;
            for _ in 0..count {
                power = power.saturating_mul(d);
            }
            factors.push(power);
        }
        d += if d == 2 { 1 } else { 2 };
    }
    if n > 1 {
        factors.push(n);
    }
    factors
}

fn solve_linear_system_mod(matrix: &[Vec<u128>], rhs: &[u128], modulus: u128) -> Option<Vec<u128>> {
    let n = matrix.len();
    if n == 0 || rhs.len() != n {
        return None;
    }
    if matrix.iter().any(|row| row.len() != n) {
        return None;
    }

    let mut aug = vec![vec![0u128; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = matrix[i][j] % modulus;
        }
        aug[i][n] = rhs[i] % modulus;
    }

    for col in 0..n {
        let mut pivot = col;
        while pivot < n {
            let entry = aug[pivot][col] % modulus;
            if entry != 0 && gcd_u128(entry, modulus) == 1 {
                break;
            }
            pivot += 1;
        }
        if pivot == n {
            return None;
        }
        if pivot != col {
            aug.swap(pivot, col);
        }

        let pivot_val = aug[col][col] % modulus;
        let inv = mod_inverse(pivot_val, modulus)?;
        for j in col..=n {
            aug[col][j] = mod_mul(aug[col][j], inv, modulus);
        }

        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col] % modulus;
            if factor == 0 {
                continue;
            }
            for j in col..=n {
                let term = mod_mul(factor, aug[col][j], modulus);
                aug[row][j] = mod_sub(aug[row][j], term, modulus);
            }
        }
    }

    Some(aug.into_iter().map(|row| row[n] % modulus).collect())
}

fn main() {
    let chan = ArithChannel::new(P, Q, DIM, N);
    let (aces, secret_key) = Aces::generate_keypair(&chan);

    println!("Parameters: p = {}, q = {}, dim = {}, n = {}", P, Q, DIM, N);
    let q_factors = factor_into_coprime(Q);
    println!("Prime-power factorisation of q: {:?}", q_factors);

    // ---- Evaluate public key pieces ----
    let secret_eval: Vec<u128> = secret_key.iter().map(|poly| chan.eval(poly) % Q).collect();
    println!("True secret evaluations s_i(1): {:?}", secret_eval);

    let f0_eval: Vec<Vec<u128>> = aces
        .f0
        .iter()
        .map(|row| row.iter().map(|poly| chan.eval(poly) % Q).collect())
        .collect();
    let f1_eval: Vec<u128> = aces.f1.iter().map(|poly| chan.eval(poly) % Q).collect();

    println!("f0 evaluations (by row): {:?}", f0_eval);
    println!("f1 evaluations: {:?}", f1_eval);

    // ---- Check e(1) = 0? ----
    let mut e_evals = Vec::with_capacity(N);
    for i in 0..N {
        let mut accum = 0u128;
        for j in 0..DIM {
            let term = mod_mul(f0_eval[i][j], secret_eval[j], Q);
            accum = mod_add(accum, term, Q);
        }
        let diff = mod_sub(f1_eval[i], accum, Q);
        e_evals.push(diff);
    }
    let unique_e: std::collections::HashSet<u128> = e_evals.iter().copied().collect();
    println!("Derived e_i(1) values: {:?}", e_evals);
    println!("All e_i(1) == 0? {}", e_evals.iter().all(|&v| v == 0));
    println!("Unique e(1) residues: {:?}", unique_e);

    // ---- Inspect ideal traces per column ----
    for j in 0..DIM {
        let mut col_gcd = Q;
        let mut column = Vec::with_capacity(N);
        for i in 0..N {
            let val = f0_eval[i][j];
            column.push(val);
            col_gcd = gcd_u128(col_gcd, val);
        }
        println!(
            "Column {} evals: {:?} → gcd with q = {}",
            j, column, col_gcd
        );
    }

    // ---- Attempt to recover secret by brute forcing e(1) = k·p ----
    println!("\nBrute forcing e(1) candidates (multiples of p)...");
    let mut best_candidate: Option<(u128, Vec<Option<Vec<u128>>>, Vec<(u128, bool)>)> = None;
    let mut best_successes = 0usize;

    for k in 0..=P {
        let e_candidate = (k * P) % Q;
        let rhs: Vec<u128> = f1_eval
            .iter()
            .map(|&value| mod_sub(value, e_candidate, Q))
            .collect();

        let mut factor_solutions: Vec<Option<Vec<u128>>> = Vec::with_capacity(q_factors.len());
        let mut local_success = Vec::with_capacity(q_factors.len());

        for &factor in &q_factors {
            match solve_linear_system_mod(&f0_eval, &rhs, factor) {
                Some(sol) => {
                    let expected: Vec<u128> = secret_eval.iter().map(|&x| x % factor).collect();
                    let recovered: Vec<u128> = sol.iter().map(|&x| x % factor).collect();
                    let success = expected == recovered;
                    factor_solutions.push(Some(sol));
                    local_success.push((factor, success));
                    println!(
                        "k = {:>2}, factor {:>12}: recovered {:?} (expected {:?}) → success = {}",
                        k, factor, recovered, expected, success
                    );
                }
                None => {
                    println!(
                        "k = {:>2}, factor {:>12}: matrix non-invertible, skipping candidate",
                        k, factor
                    );
                    factor_solutions.push(None);
                    local_success.push((factor, false));
                    continue;
                }
            }
        }

        let success_count = local_success.iter().filter(|(_, ok)| *ok).count();
        if success_count >= best_successes {
            best_successes = success_count;
            best_candidate = Some((e_candidate, factor_solutions.clone(), local_success.clone()));
        }
    }

    if let Some((e_guess, solutions, factor_success)) = best_candidate {
        let used_moduli: Vec<u128> = solutions
            .iter()
            .zip(q_factors.iter())
            .filter_map(|(sol, &modulus)| sol.as_ref().map(|_| modulus))
            .collect();

        if used_moduli.is_empty() {
            println!("No factor yielded a solvable system; cannot recover secret evaluations.");
        } else {
            let mut combined = vec![0u128; secret_eval.len()];
            for var_idx in 0..secret_eval.len() {
                let residues: Vec<u128> = solutions
                    .iter()
                    .zip(q_factors.iter())
                    .filter_map(|(sol, &modulus)| sol.as_ref().map(|v| (v[var_idx], modulus)))
                    .map(|(residue, modulus)| residue % modulus)
                    .collect();
                let moduli: Vec<u128> = solutions
                    .iter()
                    .zip(q_factors.iter())
                    .filter_map(|(sol, &modulus)| sol.as_ref().map(|_| modulus))
                    .collect();
                combined[var_idx] = crt_all(&residues, &moduli).unwrap();
            }
            let modulus_product: u128 = used_moduli.iter().product();
            println!("\nBest e(1) candidate: {}", e_guess);
            println!(
                "Recovered secret evaluations modulo {}: {:?}",
                modulus_product, combined
            );
            println!(
                "Matches ground truth modulo {}? {}",
                modulus_product,
                combined
                    .iter()
                    .zip(secret_eval.iter())
                    .all(|(&rec, &truth)| rec % modulus_product == truth % modulus_product)
                    && e_guess == e_evals[0]
            );
        }

        for (factor, success) in factor_success {
            println!(
                "Factor {:>12}: success rate = {}/{}",
                factor,
                if success { 1 } else { 0 },
                1
            );
        }
    } else {
        println!("Failed to recover a consistent secret using the tested e(1) candidates.");
    }
}
