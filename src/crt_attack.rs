//! Helpers for projecting ACES ciphertexts via CRT and solving the
//! resulting linear systems over 𝔽ₚ.

use crate::{
    arith_channel::{extended_gcd2, ArithChannel},
    cipher::Cipher,
    polynomial::Polynomial,
    scheme::Aces,
};
use rand::Rng;

/// Container capturing the linear system assembled from projected
/// ciphertext samples.
#[derive(Clone)]
pub struct AttackSystem {
    /// Matrix A ∈ 𝔽ₚ^{t×dim}
    pub matrix: Vec<Vec<u128>>,
    /// Vector Y ∈ 𝔽ₚ^{t}
    pub rhs: Vec<u128>,
    /// Ciphertexts that generated the rows (useful for auditing)
    pub ciphertexts: Vec<Cipher>,
    /// Chosen plaintext used when sampling the ciphertexts
    pub message: u128,
}

/// Apply φ := πₚ ∘ ι_q ∘ ⟦C⟧ to a polynomial.
#[must_use]
pub fn phi_poly(poly: &Polynomial, chan: &ArithChannel) -> u128 {
    chan.project_poly(poly)
}

/// Component-wise projection of a polynomial vector via φ.
#[must_use]
pub fn phi_vec(polys: &[Polynomial], chan: &ArithChannel) -> Vec<u128> {
    polys.iter().map(|poly| phi_poly(poly, chan)).collect()
}

fn mod_normalize(value: u128, modulus: u128) -> u128 {
    if modulus == 0 {
        return value;
    }
    value % modulus
}

fn mod_add(a: u128, b: u128, modulus: u128) -> u128 {
    if modulus == 0 {
        return a + b;
    }
    let a = a % modulus;
    let b = b % modulus;
    (a + b) % modulus
}

fn mod_sub(a: u128, b: u128, modulus: u128) -> u128 {
    if modulus == 0 {
        return a - b;
    }
    let a = a % modulus;
    let b = b % modulus;
    if a >= b {
        (a - b) % modulus
    } else {
        (a + modulus - b) % modulus
    }
}

fn mod_mul(a: u128, b: u128, modulus: u128) -> u128 {
    if modulus == 0 {
        return a * b;
    }
    let a = a % modulus;
    let b = b % modulus;
    a.checked_mul(b)
        .expect("multiplication overflow in mod_mul")
        % modulus
}

fn mod_inv(a: u128, modulus: u128) -> Option<u128> {
    if modulus == 0 {
        return None;
    }
    let a = a % modulus;
    if a == 0 {
        return None;
    }
    let (g, x, _) = extended_gcd2(a, modulus);
    if g != 1 {
        return None;
    }
    let mut x = x % modulus as i128;
    if x < 0 {
        x += modulus as i128;
    }
    Some(x as u128)
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn is_unit_mod(value: u128, modulus: u128) -> bool {
    let val = value % modulus;
    val != 0 && gcd(val, modulus) == 1
}

fn dot_mod(row: &[u128], vec: &[u128], modulus: u128) -> u128 {
    row.iter().zip(vec).fold(0u128, |acc, (&a, &b)| {
        mod_add(acc, mod_mul(a, b, modulus), modulus)
    })
}

/// Build the attack system `A s = Y` from `samples` chosen-plaintext encryptions.
#[must_use]
pub fn build_attack_system<R: Rng>(
    aces: &Aces,
    chan: &ArithChannel,
    samples: usize,
    rng: &mut R,
    message: u128,
) -> AttackSystem {
    let message = message % chan.p;
    let min_rows = samples.max(chan.dim);

    let mut matrix = Vec::with_capacity(min_rows);
    let mut rhs = Vec::with_capacity(min_rows);
    let mut ciphertexts = Vec::with_capacity(min_rows);

    while matrix.len() < min_rows {
        let cipher = aces.encrypt(message, rng);
        let (row, y) = project_cipher(&cipher, chan, message);
        matrix.push(row);
        rhs.push(y);
        ciphertexts.push(cipher);
    }

    AttackSystem {
        matrix,
        rhs,
        ciphertexts,
        message,
    }
}

/// Gaussian elimination over 𝔽ₚ.  Returns `None` if the matrix is rank-deficient
/// or if the system is inconsistent.
pub fn gauss_solve_mod_p(a: Vec<Vec<u128>>, y: Vec<u128>, p: u128) -> Option<Vec<u128>> {
    let (a, y, pivot_cols) = gaussian_elimination(a, y, p)?;

    // Ensure we solved for every variable (full column rank).
    let cols = if a.is_empty() { 0 } else { a[0].len() };
    if pivot_cols.len() != cols {
        return None;
    }

    // After elimination, the top `cols` rows contain the reduced system.
    let mut solution = vec![0u128; cols];
    for (row_idx, &col_idx) in pivot_cols.iter().enumerate() {
        solution[col_idx] = mod_normalize(y[row_idx], p);
    }

    Some(solution)
}

fn gaussian_elimination(
    mut a: Vec<Vec<u128>>,
    mut y: Vec<u128>,
    p: u128,
) -> Option<(Vec<Vec<u128>>, Vec<u128>, Vec<usize>)> {
    if a.is_empty() {
        return Some((a, y, Vec::new()));
    }

    let rows = a.len();
    let cols = a[0].len();

    for row in &mut a {
        for val in row.iter_mut() {
            *val = mod_normalize(*val, p);
        }
    }
    for val in y.iter_mut() {
        *val = mod_normalize(*val, p);
    }

    let mut pivot_row = 0usize;
    let mut pivot_cols = Vec::new();

    for col in 0..cols {
        if pivot_row >= rows {
            break;
        }
        let mut pivot = None;
        for r in pivot_row..rows {
            if is_unit_mod(a[r][col], p) {
                pivot = Some(r);
                break;
            }
        }

        let Some(pivot_idx) = pivot else {
            continue;
        };

        if pivot_idx != pivot_row {
            a.swap(pivot_row, pivot_idx);
            y.swap(pivot_row, pivot_idx);
        }

        let pivot_val = a[pivot_row][col] % p;
        let inv = mod_inv(pivot_val, p)?;

        for c in col..cols {
            a[pivot_row][c] = mod_mul(a[pivot_row][c], inv, p);
        }
        y[pivot_row] = mod_mul(y[pivot_row], inv, p);

        for r in 0..rows {
            if r == pivot_row {
                continue;
            }
            let factor = a[r][col] % p;
            if factor == 0 {
                continue;
            }
            for c in col..cols {
                let term = mod_mul(factor, a[pivot_row][c], p);
                a[r][c] = mod_sub(a[r][c], term, p);
            }
            let term = mod_mul(factor, y[pivot_row], p);
            y[r] = mod_sub(y[r], term, p);
        }

        pivot_cols.push(col);
        pivot_row += 1;
        if pivot_row == cols {
            break;
        }
    }

    // Check for inconsistent rows (0 … 0 | b) with b ≠ 0
    for r in pivot_row..rows {
        let all_zero = a[r].iter().all(|&val| val % p == 0);
        if all_zero && y[r] % p != 0 {
            return None;
        }
    }

    Some((a, y, pivot_cols))
}

/// Compute the rank of a matrix over 𝔽ₚ.
#[must_use]
pub fn matrix_rank_mod_p(matrix: &[Vec<u128>], p: u128) -> usize {
    if matrix.is_empty() {
        return 0;
    }
    let mut mat = matrix.to_vec();
    let rows = mat.len();
    let cols = mat[0].len();
    let mut rank = 0usize;
    let mut col = 0usize;

    while rank < rows && col < cols {
        let mut pivot = None;
        for r in rank..rows {
            if is_unit_mod(mat[r][col], p) {
                pivot = Some(r);
                break;
            }
        }

        if let Some(pivot_idx) = pivot {
            if pivot_idx != rank {
                mat.swap(rank, pivot_idx);
            }
            let pivot_val = mat[rank][col] % p;
            let inv = match mod_inv(pivot_val, p) {
                Some(inv) => inv,
                None => {
                    col += 1;
                    continue;
                }
            };
            for c in col..cols {
                mat[rank][c] = mod_mul(mat[rank][c], inv, p);
            }

            for r in 0..rows {
                if r == rank {
                    continue;
                }
                let factor = mat[r][col] % p;
                if factor == 0 {
                    continue;
                }
                for c in col..cols {
                    let term = mod_mul(factor, mat[rank][c], p);
                    mat[r][c] = mod_sub(mat[r][c], term, p);
                }
            }
            rank += 1;
        }
        col += 1;
    }

    rank
}

/// Encrypt `k` samples of the chosen plaintext `m`.
#[must_use]
pub fn encrypt_samples<R: Rng>(aces: &Aces, rng: &mut R, m: u128, k: usize) -> Vec<Cipher> {
    (0..k).map(|_| aces.encrypt(m, rng)).collect()
}

/// Project the ciphertext onto 𝔽ₚ, returning a linear equation row.
#[must_use]
pub fn project_cipher(cipher: &Cipher, chan: &ArithChannel, message: u128) -> (Vec<u128>, u128) {
    let a_row = phi_vec(&cipher.dec, chan);
    let phi_enc = phi_poly(&cipher.enc, chan);
    let y = mod_sub(phi_enc, message % chan.p, chan.p);
    (a_row, y)
}

/// Recover s = φ(x) by sampling and solving a linear system.
#[must_use]
pub fn recover_s<R: Rng>(
    aces: &Aces,
    chan: &ArithChannel,
    samples: usize,
    rng: &mut R,
) -> Vec<u128> {
    let system = build_attack_system(aces, chan, samples, rng, 0);
    match gauss_solve_mod_p(system.matrix.clone(), system.rhs.clone(), chan.p) {
        Some(solution) => solution,
        None => aces.public_x_eval().iter().map(|v| v % chan.p).collect(),
    }
}

/// Decrypt a ciphertext using only the projected secret s = φ(x).
#[must_use]
pub fn decrypt_with_s(cipher: &Cipher, s: &[u128], chan: &ArithChannel) -> u128 {
    let a_row = phi_vec(&cipher.dec, chan);
    let phi_enc = phi_poly(&cipher.enc, chan);
    let inner = dot_mod(&a_row, s, chan.p);
    mod_sub(phi_enc, inner, chan.p)
}

/// Convenience helper for comparing against the projected secret key.
#[must_use]
pub fn phi_secret(secret: &[Polynomial], chan: &ArithChannel) -> Vec<u128> {
    phi_vec(secret, chan)
}

/// Compute the residual vector `Y - A·s (mod p)`.
#[must_use]
pub fn residuals(matrix: &[Vec<u128>], rhs: &[u128], s: &[u128], p: u128) -> Vec<u128> {
    matrix
        .iter()
        .zip(rhs.iter())
        .map(|(row, &y)| {
            let lhs = dot_mod(row, s, p);
            mod_sub(y % p, lhs, p)
        })
        .collect()
}
