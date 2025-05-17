//! rand_iso.rs  ―  Random isomorphism generator for ACES
//! 可逆行列 U ∈ GL(dim, ℤ/intmodℤ) とその逆行列 U⁻¹ をランダム生成する。
//!
//! 秘密鍵ベクトル x を行列 U の列として得るときに利用する。

use rand::{Rng, thread_rng};
use std::cmp::Ordering;

/// 行列型 ― 各行は Vec<u128>
pub type Matrix = Vec<Vec<u128>>;

/// 単なるコンテナ
#[derive(Clone, Debug)]
pub struct RandIso {
    pub intmod: u128,
    pub dim:    usize,
}

impl RandIso {
    /// 新規生成
    #[must_use]
    pub fn new(intmod: u128, dim: usize) -> Self {
        assert!(intmod >= 3, "modulus must be ≥ 3");
        assert!(dim >= 2,    "dimension must be ≥ 2");
        Self { intmod, dim }
    }

    /// ランダムペア (i ≠ j)
    fn generate_pair<R: Rng>(&self, rng: &mut R) -> (usize, usize) {
        let i = rng.gen_range(0..self.dim);
        let mut j = rng.gen_range(0..self.dim);
        while j == i {
            j = rng.gen_range(0..self.dim);
        }
        (i, j)
    }

    /* ---------- 初等変換 3 種 ---------- */

    /// swap 行列 (自己逆)
    fn generate_swap<R: Rng>(&self, rng: &mut R) -> (Matrix, Matrix) {
        let (i, j) = self.generate_pair(rng);
        let mut m = identity(self.dim);
        m[i][i] = 0;
        m[j][j] = 0;
        m[i][j] = 1;
        m[j][i] = 1;
        // swap は自己逆
        (m.clone(), m)
    }

    /// diag( a₀,…,a_{dim-1} ), 逆行列は diag(a₀⁻¹,…)
    fn generate_mult<R: Rng>(&self, rng: &mut R) -> (Matrix, Matrix) {
        let mut m   = identity(self.dim);
        let mut inv = identity(self.dim);
        for r in 0..self.dim {
            let a = loop {
                let cand = rng.gen_range(1..self.intmod); // 0 でなければ OK（q は素数想定）
                if gcd(cand, self.intmod) == 1 {
                    break cand;
                }
            };
            let ainv = mod_inverse(a, self.intmod)
                .expect("modular inverse must exist");
            m[r][r]   = a;
            inv[r][r] = ainv;
        }
        (m, inv)
    }

    /// 単位行列に i 行目←i 行目 + a·j 行目 を加える（自己逆ではない）
    fn generate_line<R: Rng>(&self, rng: &mut R) -> (Matrix, Matrix) {
        let (i, j) = self.generate_pair(rng);
        let a = rng.gen_range(1..self.intmod);
        let mut m   = identity(self.dim);
        let mut inv = identity(self.dim);

        m[i][j]   = a;
        inv[i][j] = (self.intmod - a) % self.intmod; // -a  (mod intmod)

        (m, inv)
    }

    /* ---------- メイン: 合成 ---------- */

    /// 長さ `length` のランダム列で U, U⁻¹ を生成
    #[must_use]
    pub fn generate(
        &self,
        length: usize,
        pswap: usize,
        pmult: usize,
        pline: usize,
    ) -> (Matrix, Matrix) {
        assert!(pswap + pmult + pline > 0, "at least one transformation type");

        let mut rng   = thread_rng();
        let mut u     = identity(self.dim);
        let mut inv_u = identity(self.dim);

        // 重み付き選択用テーブル
        let mut choices = Vec::new();
        choices.extend(std::iter::repeat(Trans::Swap).take(pswap));
        choices.extend(std::iter::repeat(Trans::Mult).take(pmult));
        choices.extend(std::iter::repeat(Trans::Line).take(pline));

        for _ in 0..length {
            let trans = choices[rng.gen_range(0..choices.len())];
            let (m, invm) = match trans {
                Trans::Swap => self.generate_swap(&mut rng),
                Trans::Mult => self.generate_mult(&mut rng),
                Trans::Line => self.generate_line(&mut rng),
            };
            u     = mat_mul_mod(&m, &u, self.intmod);
            inv_u = mat_mul_mod(&inv_u, &invm, self.intmod); // inv は右から掛ける
        }
        (u, inv_u)
    }
}

/* ---------- 内部ユーティリティ ---------- */

#[derive(Copy, Clone)]
enum Trans { Swap, Mult, Line }

/// 単位行列
fn identity(dim: usize) -> Matrix {
    (0..dim)
        .map(|r| {
            (0..dim)
                .map(|c| if r == c { 1 } else { 0 })
                .collect()
        })
        .collect()
}

/// 行列積 (mod m)
fn mat_mul_mod(a: &Matrix, b: &Matrix, m: u128) -> Matrix {
    let dim = a.len();
    assert_eq!(dim, b.len());
    let mut out = vec![vec![0; dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            let mut acc = 0u128;
            for k in 0..dim {
                acc = (acc + a[i][k] * b[k][j]) % m;
            }
            out[i][j] = acc;
        }
    }
    out
}

/// 最大公約数 (ユークリッド互除法)
fn gcd(mut x: u128, mut y: u128) -> u128 {
    while y != 0 {
        let tmp = x % y;
        x = y;
        y = tmp;
    }
    x
}

/// 拡張ユークリッドで a⁻¹ (mod m) を求める
fn mod_inverse(a: u128, modulus: u128) -> Option<u128> {
    // 拡張 Euclid を i128 で
    let (mut t, mut new_t) = (0i128, 1i128);
    let (mut r, mut new_r) = (modulus as i128, a as i128);

    while new_r != 0 {
        let quotient = r / new_r;
        t   -= quotient * new_t;
        r   -= quotient * new_r;
        std::mem::swap(&mut t, &mut new_t);
        std::mem::swap(&mut r, &mut new_r);
    }
    if r.abs() != 1 {
        return None; // gcd != 1
    }
    if t < 0 {
        t += modulus as i128;
    }
    Some(t as u128)
}

/* ---------- テスト ---------- */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rand_iso_invertibility() {
        let iso = RandIso::new(97, 4);
        let (u, invu) = iso.generate(50, 1, 2, 3);

        // U · U⁻¹ = I  (mod intmod)
        let id = mat_mul_mod(&u, &invu, iso.intmod);
        for r in 0..iso.dim {
            for c in 0..iso.dim {
                if r == c {
                    assert_eq!(id[r][c], 1);
                } else {
                    assert_eq!(id[r][c], 0);
                }
            }
        }
    }
}
