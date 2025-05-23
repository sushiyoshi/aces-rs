// ------------------------------------------------------------
// 補助: 素朴 pow / 逆元
// ------------------------------------------------------------
#[inline]
fn mod_pow(mut base: u128, mut exp: u128, modu: u128) -> u128 {
    let mut res = 1u128;
    while exp > 0 {
        if exp & 1 == 1 {
            res = res * base % modu;
        }
        base = base * base % modu;
        exp >>= 1;
    }
    res
}

#[inline]
fn mod_inv(x: u128, modu: u128) -> u128 {
    mod_pow(x, modu - 2, modu) // q は素数仮定
}

// ------------------------------------------------------------
// primitive-root 探索
// ------------------------------------------------------------
fn factorize(mut n: u128) -> Vec<u128> {
    let mut f = Vec::new();
    let mut p = 2u128;
    while p * p <= n {
        if n % p == 0 {
            f.push(p);
            while n % p == 0 {
                n /= p;
            }
        }
        p += if p == 2 { 1 } else { 2 }; // 2 の後は奇数だけ
    }
    if n > 1 {
        f.push(n);
    }
    f
}

fn is_primitive_root(g: u128, modu: u128, factors: &[u128]) -> bool {
    for &p in factors {
        if mod_pow(g, (modu - 1) / p, modu) == 1 {
            return false;
        }
    }
    true
}

/// 任意の “NTT-friendly” 素数に対して原始根 g を返す
fn find_primitive_root(modu: u128) -> u128 {
    let phi = modu - 1;
    let factors = factorize(phi);
    let mut g = 2u128;
    loop {
        if is_primitive_root(g, modu, &factors) {
            return g;
        }
        g += 1;
    }
}

// ------------------------------------------------------------
// bit-reverse & NTT 本体
// ------------------------------------------------------------
fn bit_reverse(vec: &mut [u128]) {
    let n = vec.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            vec.swap(i, j);
        }
    }
}

/// in-place NTT / iNTT
fn ntt(a: &mut [u128], invert: bool, modu: u128, g: u128) {
    let n = a.len();
    bit_reverse(a);

    let mut len = 2;
    while len <= n {
        let w_len = mod_pow(g, (modu - 1) / len as u128, modu);
        let root = if invert { mod_inv(w_len, modu) } else { w_len };

        for i in (0..n).step_by(len) {
            let mut w = 1u128;
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = a[i + j + len / 2] * w % modu;
                a[i + j] = (u + v) % modu;
                a[i + j + len / 2] = (u + modu - v) % modu;
                w = w * root % modu;
            }
        }
        len <<= 1;
    }

    if invert {
        let inv_n = mod_inv(n as u128, modu);
        for x in a.iter_mut() {
            *x = *x * inv_n % modu;
        }
    }
}

// ------------------------------------------------------------
// 公開 API: 多項式畳み込み
// ------------------------------------------------------------
pub fn convolve(a: &[u128], b: &[u128], modu: u128) -> Vec<u128> {
    assert!(modu >= 3, "modulus must be >= 3");
    assert!(modu & 1 == 1, "modulus must be odd prime");

    // 畳み込み長
    let need = a.len() + b.len() - 1;
    let mut ntt_len = 1usize;
    while ntt_len < need {
        ntt_len <<= 1;
    }
    // 2 の冪乗チェック
    assert_eq!(ntt_len & (ntt_len - 1), 0, "length not power-of-two");

    // 原始根探索
    let g = find_primitive_root(modu);

    // コピー + 0 埋め
    let mut fa = vec![0u128; ntt_len];
    let mut fb = vec![0u128; ntt_len];
    fa[..a.len()].copy_from_slice(a);
    fb[..b.len()].copy_from_slice(b);

    // 変換・点乗・逆変換
    ntt(&mut fa, false, modu, g);
    ntt(&mut fb, false, modu, g);
    for i in 0..ntt_len {
        fa[i] = fa[i] * fb[i] % modu;
    }
    ntt(&mut fa, true, modu, g);

    // 余剰ゼロを除去
    while fa.last() == Some(&0) && fa.len() > 1 {
        fa.pop();
    }
    fa
}
