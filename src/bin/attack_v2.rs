//! Diagnostic/attack harness that targets the actual ACES parameters.
//!
//! This program instantiates the scheme with the parameters used by the
//! repository (`p = 32`, `q = 2^40 + 1`, `dim = n = 5`).  It inspects the
//! generated public key `(f0, f1)` and recovers the evaluations of the secret
//! key polynomials at the public evaluation point ω = 1 by solving the induced
//! linear system modulo each prime-power factor of `q` and combining the
//! solutions with the Chinese Remainder Theorem.
//!
//! The code is adapted from `attack.rs`, but instead of using synthetic
//! ciphertext samples it relies on the real key-generation logic from
//! [`Aces::generate_keypair`].

use aces_core::{Aces, ArithChannel};

const P: u128 = 32;
const Q: u128 = (1u128 << 40) + 1;
const DIM: usize = 5;
const N: usize = 5;
const NOISE_LEVEL: u128 = 1; // generate_vanisher(k) with k = 1 in key generation

fn mod_add(a: u128, b: u128, m: u128) -> u128 {
    ((a % m).wrapping_add(b % m)) % m
}

fn mod_sub(a: u128, b: u128, m: u128) -> u128 {
    let a = a % m;
    let b = b % m;
    if a >= b {
        a - b
    } else {
        let diff = (b - a) % m;
        if diff == 0 {
            0
        } else {
            m - diff
        }
    }
}

fn mod_mul(a: u128, b: u128, m: u128) -> u128 {
    ((a % m).wrapping_mul(b % m)) % m
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn mod_inverse(a: u128, m: u128) -> Option<u128> {
    let (g, x, _) = extended_gcd(a % m, m);
    if g != 1 {
        return None;
    }
    let mut x = x % m as i128;
    if x < 0 {
        x += m as i128;
    }
    Some(x as u128)
}

fn extended_gcd(a: u128, b: u128) -> (u128, i128, i128) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x1, y1) = extended_gcd(b, a % b);
        let x = y1;
        let y = x1 - (a / b) as i128 * y1;
        (g, x, y)
    }
}

fn factor_prime_powers(mut n: u128) -> Vec<u128> {
    let mut factors = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            let mut power = 1u128;
            while n % d == 0 {
                n /= d;
                power *= d;
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

fn crt_pair(x1: u128, m1: u128, x2: u128, m2: u128) -> Option<(u128, u128)> {
    let (g, s, t) = extended_gcd(m1, m2);
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

fn solve_linear_system_mod(matrix: &[Vec<u128>], rhs: &[u128], modulus: u128) -> Option<Vec<u128>> {
    let n = matrix.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let m = matrix[0].len();
    assert_eq!(n, m, "square matrix expected");
    assert_eq!(rhs.len(), n);

    let mut aug = vec![vec![0u128; m + 1]; n];
    for i in 0..n {
        for j in 0..m {
            aug[i][j] = matrix[i][j] % modulus;
        }
        aug[i][m] = rhs[i] % modulus;
    }

    for col in 0..m {
        let mut pivot = None;
        for row in col..n {
            let val = aug[row][col] % modulus;
            if val != 0 && gcd(val, modulus) == 1 {
                pivot = Some(row);
                break;
            }
        }
        let pivot_row = pivot?;
        if pivot_row != col {
            aug.swap(pivot_row, col);
        }
        let inv = mod_inverse(aug[col][col], modulus)?;
        for j in col..=m {
            aug[col][j] = mod_mul(aug[col][j], inv, modulus);
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            if factor == 0 {
                continue;
            }
            for j in col..=m {
                let sub = mod_mul(factor, aug[col][j], modulus);
                aug[row][j] = mod_sub(aug[row][j], sub % modulus, modulus);
            }
        }
    }

    let mut solution = vec![0u128; m];
    for i in 0..m {
        solution[i] = aug[i][m] % modulus;
    }
    Some(solution)
}

fn solve_columns(
    matrix: &[Vec<u128>],
    rhs: &[u128],
    columns: &[usize],
    modulus: u128,
) -> Option<Vec<u128>> {
    if columns.is_empty() {
        return Some(Vec::new());
    }
    let target = columns.len();
    if matrix.len() < target {
        return None;
    }

    fn dfs(
        matrix: &[Vec<u128>],
        rhs: &[u128],
        columns: &[usize],
        modulus: u128,
        rows: &mut Vec<usize>,
        start: usize,
    ) -> Option<Vec<u128>> {
        if rows.len() == columns.len() {
            let size = columns.len();
            let mut submatrix = vec![vec![0u128; size]; size];
            let mut subrhs = vec![0u128; size];
            for (i, &row) in rows.iter().enumerate() {
                subrhs[i] = rhs[row] % modulus;
                for (j, &col) in columns.iter().enumerate() {
                    submatrix[i][j] = matrix[row][col] % modulus;
                }
            }
            if let Some(solution) = solve_linear_system_mod(&submatrix, &subrhs, modulus) {
                for (row_idx, row) in matrix.iter().enumerate() {
                    let mut acc = 0u128;
                    for (col_idx, &col) in columns.iter().enumerate() {
                        acc = mod_add(
                            acc,
                            mod_mul(row[col] % modulus, solution[col_idx], modulus),
                            modulus,
                        );
                    }
                    if acc != rhs[row_idx] % modulus {
                        return None;
                    }
                }
                return Some(solution);
            }
            return None;
        }

        for next in start..matrix.len() {
            rows.push(next);
            if let Some(sol) = dfs(matrix, rhs, columns, modulus, rows, next + 1) {
                return Some(sol);
            }
            rows.pop();
        }
        None
    }

    let mut rows = Vec::with_capacity(target);
    dfs(matrix, rhs, columns, modulus, &mut rows, 0)
}

fn main() {
    println!("=== ACES attack v2 (actual parameters) ===");
    let chan = ArithChannel::new(P, Q, DIM, N);
    let (aces, secret_key) = Aces::generate_keypair(&chan);

    println!("Parameters: p = {}, q = {}, dim = {}, n = {}", P, Q, DIM, N);

    let q_factors = factor_prime_powers(Q);
    println!("Prime-power factors of q: {:?}", q_factors);

    let secret_evals: Vec<u128> = secret_key.iter().map(|poly| chan.eval(poly) % Q).collect();
    println!(
        "Ground-truth secret evaluations at ω = 1: {:?}",
        secret_evals
    );

    let mut matrix = vec![vec![0u128; DIM]; N];
    let mut rhs = vec![0u128; N];
    let mut noise_evals = Vec::with_capacity(N);

    for i in 0..N {
        let mut expected = 0u128;
        for j in 0..DIM {
            let a_eval = chan.eval(&aces.f0[i][j]) % Q;
            matrix[i][j] = a_eval;
            expected = mod_add(expected, mod_mul(a_eval, secret_evals[j], Q), Q);
        }
        let b_eval = chan.eval(&aces.f1[i]) % Q;
        rhs[i] = b_eval;
        let noise = mod_sub(b_eval, expected, Q);
        noise_evals.push(noise);
    }

    println!("Noise evaluations (diagnostic): {:?}", noise_evals);

    let mut column_factors = vec![0u128; DIM];
    for j in 0..DIM {
        let mut col_gcd = Q;
        for i in 0..N {
            col_gcd = gcd(col_gcd, matrix[i][j]);
        }
        column_factors[j] = gcd(col_gcd, Q);
    }
    println!(
        "Column gcd factors (expected q_sigma assignments): {:?}",
        column_factors
    );

    let noise_eval = (NOISE_LEVEL * P) % Q;
    println!("Assumed public noise evaluation per row: {}", noise_eval);

    use std::collections::BTreeMap;
    let mut groups: BTreeMap<u128, Vec<usize>> = BTreeMap::new();
    for (idx, &factor) in column_factors.iter().enumerate() {
        groups.entry(factor).or_default().push(idx);
    }

    let rhs_without_noise: Vec<u128> = rhs.iter().map(|&val| mod_sub(val, noise_eval, Q)).collect();

    let mut residues: Vec<Vec<Option<u128>>> = vec![vec![None; DIM]; q_factors.len()];

    for (mod_idx, &qi) in q_factors.iter().enumerate() {
        let rhs_mod: Vec<u128> = rhs_without_noise.iter().map(|&val| val % qi).collect();
        let active_cols: Vec<usize> = (0..DIM)
            .filter(|&j| gcd(column_factors[j], qi) == 1)
            .collect();
        if active_cols.is_empty() {
            println!("Modulo {}: no informative columns (all annihilated).", qi);
            continue;
        }
        match solve_columns(&matrix, &rhs_mod, &active_cols, qi) {
            Some(solution) => {
                println!(
                    "Recovered residues modulo {} for columns {:?}: {:?}",
                    qi, active_cols, solution
                );
                for (idx, &col) in active_cols.iter().enumerate() {
                    residues[mod_idx][col] = Some(solution[idx]);
                }
            }
            None => {
                println!("Failed to solve system modulo {}", qi);
            }
        }
    }

    let idx_small = q_factors.iter().position(|&m| m == 257);
    let idx_large = q_factors.iter().position(|&m| m == 427_825_5361);

    if let (Some(idx_small), Some(_idx_large)) = (idx_small, idx_large) {
        if let (Some(group_small), Some(group_large)) =
            (groups.get(&257), groups.get(&427_825_5361))
        {
            if !group_small.is_empty() && !group_large.is_empty() {
                let can_adjust = group_large
                    .iter()
                    .all(|&col| residues[idx_small][col].is_some());
                if can_adjust {
                    let mut reduced_matrix = vec![vec![0u128; group_small.len()]; N];
                    let mut reduced_rhs = vec![0u128; N];
                    for i in 0..N {
                        let mut rem = rhs_without_noise[i];
                        for &col in group_large {
                            if let Some(residue) = residues[idx_small][col] {
                                let term = mod_mul(matrix[i][col], residue, Q);
                                rem = mod_sub(rem, term, Q);
                            }
                        }
                        if rem % 257 != 0 {
                            println!(
                                "Warning: remainder {} not divisible by 257 at row {}",
                                rem, i
                            );
                        }
                        reduced_rhs[i] = rem / 257;
                        for (pos, &col) in group_small.iter().enumerate() {
                            reduced_matrix[i][pos] = matrix[i][col] / 257;
                        }
                    }
                    let matrix_mod257: Vec<Vec<u128>> = reduced_matrix
                        .iter()
                        .map(|row| row.iter().map(|&val| val % 257).collect())
                        .collect();
                    println!("Reduced matrix mod 257: {:?}", matrix_mod257);
                    println!("Reduced rhs (before mod): {:?}", reduced_rhs);
                    let columns_idx: Vec<usize> = (0..group_small.len()).collect();
                    let rhs_mod_small: Vec<u128> =
                        reduced_rhs.iter().map(|&val| val % 257).collect();
                    if let Some(solution_small) =
                        solve_columns(&reduced_matrix, &rhs_mod_small, &columns_idx, 257)
                    {
                        println!(
                            "Derived additional residues modulo 257 for columns {:?}: {:?}",
                            group_small, solution_small
                        );
                        for (idx, &col) in group_small.iter().enumerate() {
                            residues[idx_small][col] = Some(solution_small[idx]);
                        }
                    } else {
                        println!(
                            "Failed to derive modulo-257 residues for columns {:?}",
                            group_small
                        );
                    }
                } else {
                    println!(
                        "Insufficient data to eliminate large-factor columns when targeting modulo 257."
                    );
                }
            }
        }
    }

    let mut combined = vec![None; DIM];
    for j in 0..DIM {
        let mut residues_list = Vec::new();
        for (mod_idx, &modulus) in q_factors.iter().enumerate() {
            if let Some(value) = residues[mod_idx][j] {
                residues_list.push((value, modulus));
            }
        }
        if residues_list.len() == q_factors.len() {
            let (vals, mods): (Vec<_>, Vec<_>) = residues_list.into_iter().unzip();
            if let Some(value) = crt_all(&vals, &mods) {
                combined[j] = Some(value % Q);
            }
        }
    }

    println!("Recovered residues per modulus: {:?}", residues);
    println!("CRT reconstruction (when possible): {:?}", combined);

    for (idx, maybe_value) in combined.iter().enumerate() {
        match maybe_value {
            Some(value) => println!(
                "coordinate {} ⇒ recovered {} (actual {})",
                idx, value, secret_evals[idx]
            ),
            None => println!(
                "coordinate {} ⇒ partial info only (actual {})",
                idx, secret_evals[idx]
            ),
        }
    }
}
