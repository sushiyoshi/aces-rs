//! locator.rs
//! ------------------------------------------------------------
//! 秘密鍵を使わない p-Refreshable 判定用ユーティリティ
//! （ACES §5.4 Remark 5.45, Thm 5.46）
//
//! ここでは概念実証として：
//!   * locator として **ゼロベクトル** 1 本だけを採用
//!   * director として     ±e_i  （標準基底ベクトル）を採用
//! をデータベース A = [ℓ0 | ±e₁ | … | ±e_n] に格納し，
//! 線形方程式  A·d = ℓ  (mod q) を Gauss 法で解ければ
//! “locator/director 分解が見つかった” とみなし Refreshable と判定する。
//!
//! （注意）このデータベースは極めてミニマルです。実運用では
//!   * 余裕を大きく取った複数の locator
//!   * ランダム／構造化 director 収集
//!   * 係数 d ∈ {–1,0,1}ⁿ に制限するソルバ
//! などに置き換えてください。

use num_integer::Integer;

/// q 上 Zq-ベクトル
pub type Vector = Vec<u128>;
/// Zq 行列
pub type Matrix = Vec<Vector>;

/// ★ Gauss 消去（mod q）
/// 与えられた A (m×n) と rhs (m) について，A·x = rhs の解を 1 つ返す。
/// （解が存在しない場合は None）
pub fn solve_mod_q(a: &Matrix, rhs: &Vector, q: u128) -> Option<Vector> {
    let m = a.len();
    let n = a[0].len();
    assert_eq!(rhs.len(), m);

    // 拡大行列 [A | rhs]
    let mut aug: Vec<Vec<u128>> = a
        .iter()
        .zip(rhs)
        .map(|(row, &b)| {
            let mut v = row.clone();
            v.push(b % q);
            v
        })
        .collect();

    let mut col = 0usize;
    for row in 0..m {
        // pivot 探索
        let mut pivot = row;
        while pivot < m && aug[pivot][col] == 0 {
            pivot += 1;
        }
        if pivot == m {
            col += 1;
            if col == n {
                break;
            }
            continue;
        }
        aug.swap(row, pivot);

        // 逆元
        let inv = modinv(aug[row][col], q)?;
        for j in col..=n {
            aug[row][j] = (aug[row][j] * inv) % q;
        }

        // 他行を消去
        for i in 0..m {
            if i == row {
                continue;
            }
            let factor = aug[i][col];
            for j in col..=n {
                let t = (factor * aug[row][j]) % q;
                aug[i][j] = (aug[i][j] + q - t) % q;
            }
        }
        col += 1;
        if col == n {
            break;
        }
    }

    // 一貫性チェック & 解抽出
    let mut sol = vec![0u128; n];
    for i in 0..m {
        // 係数が全部ゼロで rhs≠0 → 不整合
        if (0..n).all(|c| aug[i][c] == 0) && aug[i][n] != 0 {
            return None;
        }
        // 基底列（基本変数）を決定
        if let Some(pivot_col) = (0..n).find(|&c| aug[i][c] == 1) {
            sol[pivot_col] = aug[i][n];
        }
    }
    Some(sol)
}

/// 拡張 Euclid – 逆元（q は 64 bit 未満素数想定）
fn modinv(a: u128, q: u128) -> Option<u128> {
    let (mut t, mut new_t) = (0i128, 1i128);
    let (mut r, mut new_r) = (q as i128, a as i128);
    while new_r != 0 {
        let qout = r / new_r;
        (t, new_t) = (new_t, t - qout * new_t);
        (r, new_r) = (new_r, r - qout * new_r);
    }
    if r != 1 {
        return None; // 非可逆
    }
    if t < 0 {
        t += q as i128;
    }
    Some(t as u128)
}

/// ------------------------------------------------------------
/// 公開 DB ＝ [locator | ±director] 行列
pub struct LocatorDB {
    pub a: Matrix, // (n)×(1+2n)
    pub q: u128,
}

impl LocatorDB {
    /// デフォルト生成（locator=0, directors=±e_i）
    pub fn new(dim: usize, q: u128) -> Self {
        let mut a: Matrix = vec![vec![0; 1 + 2 * dim]; dim];
        // locator 列 = 0 ベクトル → 既にゼロ
        // director 列
        for i in 0..dim {
            // +e_i
            a[i][1 + i] = 1 % q;
            // -e_i
            a[i][1 + dim + i] = (q - 1) % q;
        }
        Self { a, q }
    }

    /// 与えられた ℓ が A·d = ℓ (mod q) に分解可能なら true
    pub fn has_decomposition(&self, ell: &Vector) -> bool {
        solve_mod_q(&self.a, ell, self.q).is_some()
    }
}
