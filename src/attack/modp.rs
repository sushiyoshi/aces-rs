//! Lightweight modular arithmetic helpers over ℤ/pℤ used by the attack.

use crate::arith_channel::extended_gcd2;

/// Compute modular inverse of `a` modulo `p` when gcd(a, p) = 1.
#[must_use]
pub fn inv_mod(a: u128, p: u128) -> Option<u128> {
    let a = a % p;
    if a == 0 {
        return None;
    }
    let (g, x, _) = extended_gcd2(a, p);
    if g != 1 {
        None
    } else {
        let mut res = x % p as i128;
        if res < 0 {
            res += p as i128;
        }
        Some(res as u128)
    }
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

fn mul_mod(a: u128, b: u128, p: u128) -> u128 {
    ((a % p) * (b % p)) % p
}

fn sub_mod(a: u128, b: u128, p: u128) -> u128 {
    (a + p - (b % p)) % p
}

/// Perform Gaussian elimination over ℤ/pℤ to solve A·s = y.
pub fn gaussian_elimination_mod_p(
    matrix: &[Vec<u128>],
    rhs: &[u128],
    p: u128,
    dimension: usize,
) -> Option<Vec<u128>> {
    if dimension == 0 || matrix.is_empty() || rhs.len() != matrix.len() {
        return None;
    }

    let rows = matrix.len();
    let cols = dimension;

    if matrix.iter().any(|row| row.len() < cols) {
        return None;
    }

    // Build augmented matrix [A | y] reduced modulo p.
    let mut aug = vec![vec![0u128; cols + 1]; rows];
    for (r, row) in matrix.iter().enumerate() {
        for c in 0..cols {
            aug[r][c] = row[c] % p;
        }
        aug[r][cols] = rhs[r] % p;
    }

    let mut pivot_row = 0usize;
    let mut pivot_cols = Vec::with_capacity(cols);

    for col in 0..cols {
        if pivot_row == rows {
            break;
        }
        let mut pivot = None;
        for r in pivot_row..rows {
            let val = aug[r][col] % p;
            if val != 0 && gcd(val, p) == 1 {
                pivot = Some(r);
                break;
            }
        }

        let Some(pivot_idx) = pivot else { continue };

        if pivot_idx != pivot_row {
            aug.swap(pivot_row, pivot_idx);
        }

        let inv = inv_mod(aug[pivot_row][col], p)?;
        for c in col..=cols {
            aug[pivot_row][c] = mul_mod(aug[pivot_row][c], inv, p);
        }

        for r in 0..rows {
            if r == pivot_row {
                continue;
            }
            let factor = aug[r][col] % p;
            if factor == 0 {
                continue;
            }
            for c in col..=cols {
                let prod = mul_mod(factor, aug[pivot_row][c], p);
                aug[r][c] = sub_mod(aug[r][c], prod, p);
            }
        }

        pivot_cols.push(col);
        pivot_row += 1;
    }

    if pivot_cols.len() < cols {
        return None;
    }

    let mut solution = vec![0u128; cols];
    for (row, &col) in pivot_cols.iter().enumerate() {
        solution[col] = aug[row][cols] % p;
    }

    Some(solution)
}

/// Dot product modulo `p`.
#[must_use]
pub fn dot_mod_p(a: &[u128], b: &[u128], p: u128) -> u128 {
    assert_eq!(a.len(), b.len(), "dimension mismatch in dot_mod_p");
    let mut acc = 0u128;
    for (&ai, &bi) in a.iter().zip(b.iter()) {
        acc = (acc + mul_mod(ai, bi, p)) % p;
    }
    acc
}
