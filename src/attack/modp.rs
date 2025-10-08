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
    matrix: &mut [Vec<u128>],
    rhs: &mut [u128],
    p: u128,
    dimension: usize,
) -> Option<Vec<u128>> {
    if matrix.is_empty() || dimension == 0 {
        return None;
    }
    let rows = matrix.len();
    let cols = dimension;

    for row in matrix.iter_mut() {
        if row.len() < cols {
            return None;
        }
        for val in row.iter_mut() {
            *val %= p;
        }
    }
    for r in rhs.iter_mut() {
        *r %= p;
    }

    let mut pivot_row = 0;
    for col in 0..cols {
        if pivot_row >= rows {
            break;
        }
        let mut pivot = None;
        for r in pivot_row..rows {
            let val = matrix[r][col] % p;
            if val != 0 && gcd(val, p) == 1 {
                pivot = Some(r);
                break;
            }
        }
        let Some(pivot_idx) = pivot else { continue };
        if pivot_idx != pivot_row {
            matrix.swap(pivot_row, pivot_idx);
            rhs.swap(pivot_row, pivot_idx);
        }
        let pivot_val = matrix[pivot_row][col] % p;
        let inv = inv_mod(pivot_val, p)?;
        for c in col..cols {
            matrix[pivot_row][c] = mul_mod(matrix[pivot_row][c], inv, p);
        }
        rhs[pivot_row] = mul_mod(rhs[pivot_row], inv, p);

        for r in 0..rows {
            if r == pivot_row {
                continue;
            }
            let factor = matrix[r][col] % p;
            if factor == 0 {
                continue;
            }
            for c in col..cols {
                let prod = mul_mod(factor, matrix[pivot_row][c], p);
                matrix[r][c] = sub_mod(matrix[r][c], prod, p);
            }
            let prod = mul_mod(factor, rhs[pivot_row], p);
            rhs[r] = sub_mod(rhs[r], prod, p);
        }
        pivot_row += 1;
        if pivot_row == cols {
            break;
        }
    }

    if pivot_row < cols {
        return None;
    }

    let mut solution = vec![0u128; cols];
    for i in 0..cols {
        solution[i] = rhs[i] % p;
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
