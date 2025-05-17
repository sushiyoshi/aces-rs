//! Polynomial type over Z/qZ with u128 coefficients.

use rand::Rng;
use std::ops::{Add, Mul, Neg, Rem, Sub};
use crate::arith_channel::extended_gcd2;
use crate::ntt;
use std::cmp::max;         // （既にあれば重複削除）

/// f(x) = coeffs[0] + coeffs[1]·x + ...  (always mod `modulus`)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial {
    pub coeffs: Vec<u128>,
    pub modulus: u128,  // Made public since ACES uses it directly
}

impl Polynomial {
    pub fn new(coeffs: Vec<u128>, modulus: u128) -> Self {
        assert!(modulus > 0, "modulus must be positive");
        let mut coeffs = coeffs.into_iter().map(|x| x % modulus).collect::<Vec<_>>();
        Self { coeffs, modulus }
    }
    /// Construct zero polynomial
    pub fn zero(modulus: u128) -> Self {
        Self {
            coeffs: vec![0],
            modulus,
        }
    }

    /// Construct zero polynomial with specified dimension
    pub fn zero_with_dim(modulus: u128, dim: usize) -> Self {
        assert!(dim > 0, "dimension must be positive");
        let mut coeffs = vec![0; dim];
        coeffs[0] = 0; // Ensure constant term is zero
        Self {
            coeffs,
            modulus,
        }
    }

    /// Create a constant polynomial with single coefficient
    pub fn constant(value: u128, modulus: u128) -> Self {
        Self {
            coeffs: vec![value % modulus],
            modulus,
        }
    }

    /// Random polynomial of degree up to `dim-1`
    pub fn random<R: Rng>(modulus: u128, dim: usize, rng: &mut R) -> Self {
        assert!(dim > 0, "dimension must be positive");
        let mut coeffs = Vec::with_capacity(dim);
        for _ in 0..dim {
            coeffs.push(rng.gen_range(0..modulus));
        }
        let mut result = Self { coeffs, modulus };
        result.trim();
        result
    }

    /// Shift constant term so eval_at_one() == target
    // pub fn randshift<R: Rng>(target: u128, modulus: u128, dim: usize, rng: &mut R) -> Self {
    //     assert!(dim > 0, "dimension must be positive");
    //     let mut p = Self::random(modulus, dim, rng);
    //     let eval = p.eval_at_one();
    //     let target = target % modulus;
        
    //     // Calculate shift needed to achieve target evaluation
    //     // Calculate shift needed to achieve target evaluation
    //     let shift = if eval <= target {
    //         // Case 1: need to add positive shift
    //         let diff = target.checked_sub(eval)
    //             .expect("Arithmetic underflow in randshift while calculating target-eval");
    //         diff
    //     } else {
    //         // Case 2: need to subtract and ensure positive result modulo
    //         let diff = eval.checked_sub(target)
    //             .expect("Arithmetic underflow in randshift while calculating eval-target");
    //         modulus.checked_sub(diff % modulus)
    //             .expect("Arithmetic underflow in randshift while calculating modulus-diff")
    //     };

    //     // Add shift to constant term
    //     let new_coeff = p.coeffs[0].checked_add(shift)
    //         .expect("Arithmetic overflow in randshift while adding shift to coefficient");
    //     p.coeffs[0] = new_coeff % modulus;
    //     p
    // }
    // pub fn randshift<R: Rng>(coef: u128, modulus: u128, dim: usize, rng: &mut R) -> Self {
    //     assert!(dim > 0, "dimension must be positive");
        
    //     // ランダムな次数（0〜dim-1）を選ぶ
    //     let degree = rng.gen_range(0..dim);
    
    //     // 全係数を0で初期化し、ランダムな次数の位置にcoef % modulus を代入
    //     let mut coeffs = vec![0; degree + 1];
    //     coeffs[degree] = coef % modulus;
    
    //     Self { coeffs, modulus }
    // }

    pub fn trim(&mut self) {
        while self.coeffs.len() > 1 && *self.coeffs.last().unwrap() == 0 {
            self.coeffs.pop();
        }
    }

    pub fn randshift<R: Rng>(mut coef: u128, modulus: u128, dim: usize, rng: &mut R) -> Self {
        assert!(dim > 0, "dimension must be positive");
    
        let d = rng.gen_range(0..dim);

        // Guarantee a non-zero coefficient modulo `modulus`
        if coef % modulus == 0 {
            coef = rng.gen_range(1..modulus);        // re-sample uniformly
        }

        let mut coeffs = vec![0; d + 1];
        coeffs[d] = coef % modulus;
    
        Self { coeffs, modulus }
    }

    /// Helper function for tests and cases where thread_rng is acceptable
    pub fn randshift_with_thread_rng(target: u128, modulus: u128, dim: usize) -> Self {
        Self::randshift(target, modulus, dim, &mut rand::thread_rng())
    }

    /// Evaluate polynomial at x = value
    pub fn eval(&self, value: u128) -> u128 {
        let mut result = 0;
        for &coeff in self.coeffs.iter().rev() {
            result = (result * value + coeff) % self.modulus;
        }
        result
    }

    /// Fast path evaluation at x = 1
    pub fn eval_at_one(&self) -> u128 {
        self.coeffs.iter().fold(0, |acc, &c| (acc + c) % self.modulus)
    }

    /// Change modulus of coefficients
    pub fn with_modulus(&self, new_modulus: Option<u128>) -> Self {
        match new_modulus {
            Some(m) => Self {
                coeffs: self.coeffs.iter().map(|&x| x % m).collect(),
                modulus: m,
            },
            None => Self {
                coeffs: self.coeffs.clone(),
                modulus: self.modulus,
            }
        }
    }

    /// Polynomial division with remainder
    // pub fn rem_poly(&self, u: &Polynomial) -> Self {
    //     assert_eq!(self.modulus, u.modulus, "modulus mismatch in rem_poly");
    //     let mut r = self.clone();
    //     let deg_u = u.degree();
        
    //     // Handle special cases
    //     if deg_u == 0 {
    //         if u.coeffs[0] == 0 {
    //             panic!("division by zero polynomial");
    //         }
    //         // Division by constant polynomial
    //         let inv = Self::mod_inverse(u.coeffs[0], u.modulus)
    //             .expect("constant term must be invertible");
    //         return Self {
    //             coeffs: vec![(self.eval_at_one() * inv) % self.modulus],
    //             modulus: self.modulus,
    //         };
    //     }

    //     // Perform polynomial long division
    //     while r.degree() >= deg_u {
    //         let k = r.degree() - deg_u;
            
    //         // Calculate quotient term
    //         let leading_r = r.coeffs[r.degree()];
    //         let leading_u = u.coeffs[deg_u];
    //         let inv = Self::mod_inverse(leading_u, u.modulus)
    //             .expect("leading coefficient must be invertible");
    //         let q = (leading_r * inv) % u.modulus;
            
    //         // Subtract q * u * x^k
    //         for i in 0..=deg_u {
    //             let idx = i + k;
    //             let sub = (u.coeffs[i] * q) % r.modulus;
    //             r.coeffs[idx] = (r.coeffs[idx] + r.modulus - sub) % r.modulus;
    //         }
            
    //         // Normalize result
    //         while r.coeffs.len() > 1 && r.coeffs.last().unwrap() == &0 {
    //             r.coeffs.pop();
    //         }
    //     }
        
    //     // Ensure the result has proper dimension
    //     while r.coeffs.len() < deg_u {
    //         r.coeffs.push(0);
    //     }
        
    //     r
    // }
    pub fn rem_poly(&self, u: &Polynomial) -> Self {
        assert_eq!(
            self.modulus, u.modulus,
            "modulus mismatch in rem_poly"
        );
        let mut r = self.clone();
        let deg_u = u.degree();

        if deg_u == 0 {
            if u.coeffs[0] == 0 {
                panic!("division by zero polynomial");
            }
            // Division by constant polynomial
            let inv = Self::mod_inverse(u.coeffs[0], u.modulus)
                .expect("constant term must be invertible");
            // Initialize result with proper dimension
            let mut coeffs = vec![0; u.coeffs.len()];
            coeffs[0] = (self.eval_at_one() * inv) % self.modulus;
            return Self {
                coeffs,
                modulus: self.modulus,
            };
        }

        while r.degree() >= deg_u {
            let k = r.degree() - deg_u;
            
            // Calculate quotient term
            let leading_r = r.coeffs[r.degree()];
            let leading_u = u.coeffs[deg_u];
            let inv = Self::mod_inverse(leading_u, u.modulus)
                .expect("leading coefficient must be invertible");
            let q = (leading_r * inv) % u.modulus;
            
            // Subtract q * u * x^k
            for i in 0..=deg_u {
                let idx = i + k;
                let sub = (u.coeffs[i] * q) % r.modulus;
                r.coeffs[idx] = (r.coeffs[idx] + r.modulus - sub) % r.modulus;
            }
            
            // Normalize result
            while r.coeffs.len() > 1 && r.coeffs.last().unwrap() == &0 {
                r.coeffs.pop();
            }
        }
        r.trim();
        r
    }

    /// Calculate modular multiplicative inverse using extended GCD
    fn mod_inverse(a: u128, modulus: u128) -> Option<u128> {
        let (g, x, _) = extended_gcd2(a, modulus);
        if g != 1 {
            None
        } else {
            // Convert signed x to positive residue modulo m
            let x = x % modulus as i128;
            let x = if x < 0 { x + modulus as i128 } else { x };
            Some(x as u128)
        }
    }

    /// Calculate polynomial degree (highest non-zero coefficient)
    pub fn degree(&self) -> usize {
        let mut d = self.coeffs.len() - 1;
        while d > 0 && self.coeffs[d] == 0 {
            d -= 1;
        }
        d
    }
}
impl Add for &Polynomial {
    type Output = Polynomial;
    fn add(self, rhs: Self) -> Self::Output {
        let len = self.coeffs.len().max(rhs.coeffs.len());
        let mut v = vec![0; len];
        for i in 0..len {
            let a = self.coeffs.get(i).copied().unwrap_or(0);
            let b = rhs.coeffs.get(i).copied().unwrap_or(0);
            v[i] = (a + b) % self.modulus;
        }
        let mut result = Polynomial {
            coeffs: v,
            modulus: self.modulus,
        };
        result.trim();
        result
    }
}

impl Add for Polynomial {
    type Output = Polynomial;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl Add<&Polynomial> for Polynomial {
    type Output = Polynomial;
    fn add(self, rhs: &Polynomial) -> Self::Output {
        &self + rhs
    }
}

impl Add<Polynomial> for &Polynomial {
    type Output = Polynomial;
    fn add(self, rhs: Polynomial) -> Self::Output {
        self + &rhs
    }
}

// impl Mul for &Polynomial {
//     type Output = Polynomial;
//     fn mul(self, rhs: Self) -> Self::Output {
//         // println!("self, rhs: {:?}", (self, rhs));
//         let len = self.coeffs.len() + rhs.coeffs.len() - 1;
//         let mut v = vec![0; len];
//         for (i, &a) in self.coeffs.iter().enumerate() {
//             for (j, &b) in rhs.coeffs.iter().enumerate() {
//                 let idx = i + j;
//                 let prod = a.checked_mul(b).expect("polynomial multiplication overflow");
//                 v[idx] = (v[idx] + prod) % self.modulus;
                
//             }
//         }
//         let mut result = Polynomial {
//             coeffs: v,
//             modulus: self.modulus,
//         };
//         result.trim();
//         result
//     }
// }

// impl Mul for Polynomial {
//     type Output = Polynomial;
//     fn mul(self, rhs: Self) -> Self::Output {
//         &self * &rhs
//     }
// }

// impl Mul<&Polynomial> for Polynomial {
//     type Output = Polynomial;
//     fn mul(self, rhs: &Polynomial) -> Self::Output {
//         &self * rhs
//     }
// }

// impl Mul<Polynomial> for &Polynomial {
//     type Output = Polynomial;
//     fn mul(self, rhs: Polynomial) -> Self::Output {
//         self * &rhs
//     }
// }


impl<'a, 'b> Mul<&'b Polynomial> for &'a Polynomial {
    type Output = Polynomial;

    /// 高速 NTT + フォールバック O(n²) 乗算（参照 × 参照）
    fn mul(self, rhs: &'b Polynomial) -> Polynomial {
        assert_eq!(self.modulus, rhs.modulus, "moduli must match");
        let modu = self.modulus;
        let threshold = 32;

        let coeffs = if max(self.coeffs.len(), rhs.coeffs.len()) > threshold {
            ntt::convolve(&self.coeffs, &rhs.coeffs, modu)
        } else {
            let mut prod = vec![0u128; self.coeffs.len() + rhs.coeffs.len() - 1];
            for i in 0..self.coeffs.len() {
                for j in 0..rhs.coeffs.len() {
                    prod[i + j] =
                        (prod[i + j] + self.coeffs[i] * rhs.coeffs[j]) % modu;
                }
            }
            prod
        };

        // ← *** Self::new ではなく `Polynomial::new` に変更 ***
        Polynomial::new(coeffs, modu)
    }
}

impl<'a> Mul<&'a Polynomial> for Polynomial {
    type Output = Polynomial;
    fn mul(self, rhs: &'a Polynomial) -> Polynomial { (&self).mul(rhs) }
}

impl<'a> Mul<Polynomial> for &'a Polynomial {
    type Output = Polynomial;
    fn mul(self, rhs: Polynomial) -> Polynomial { self.mul(&rhs) }
}

impl Mul<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn mul(self, rhs: Polynomial) -> Polynomial { (&self).mul(&rhs) }
}


impl Neg for Polynomial {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Polynomial {
            coeffs: self.coeffs.iter().map(|&x| (self.modulus - x) % self.modulus).collect(),
            modulus: self.modulus,
        }
    }
}

impl Sub for &Polynomial {
    type Output = Polynomial;
    fn sub(self, rhs: Self) -> Self::Output {
        let len = self.coeffs.len().max(rhs.coeffs.len());
        let mut v = vec![0; len];
        for i in 0..len {
            let a = self.coeffs.get(i).copied().unwrap_or(0);
            let b = rhs.coeffs.get(i).copied().unwrap_or(0);
            v[i] = (self.modulus + a - b) % self.modulus;
        }
        let mut result = Polynomial {
            coeffs: v,
            modulus: self.modulus,
        };
        result.trim();
        result
    }
}

impl Sub for Polynomial {
    type Output = Polynomial;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl Sub<&Polynomial> for Polynomial {
    type Output = Polynomial;
    fn sub(self, rhs: &Polynomial) -> Self::Output {
        &self - rhs
    }
}

impl Sub<Polynomial> for &Polynomial {
    type Output = Polynomial;
    fn sub(self, rhs: Polynomial) -> Self::Output {
        self - &rhs
    }
}

// Basic implementation for references
impl Rem<&Polynomial> for &Polynomial {
    type Output = Polynomial;
    fn rem(self, rhs: &Polynomial) -> Self::Output {
        self.rem_poly(rhs)
    }
}

// All other implementations delegate to the reference implementation
impl Rem<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn rem(self, rhs: Polynomial) -> Self::Output {
        &self % &rhs
    }
}

impl Rem<&Polynomial> for Polynomial {
    type Output = Polynomial;
    fn rem(self, rhs: &Polynomial) -> Self::Output {
        &self % rhs
    }
}

impl Rem<Polynomial> for &Polynomial {
    type Output = Polynomial;
    fn rem(self, rhs: Polynomial) -> Self::Output {
        self % &rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, thread_rng};

    #[test]
    fn test_basic_ops() {
        let m = 16;
        let p1 = Polynomial {
            coeffs: vec![1, 2, 3],
            modulus: m,
        };
        let p2 = Polynomial {
            coeffs: vec![4, 5, 6],
            modulus: m,
        };

        // Add
        let sum = &p1 + &p2;
        assert_eq!(sum.coeffs, vec![5, 7, 9]);

        // Sub
        let diff = &p1 - &p2;
        assert_eq!(diff.coeffs, vec![13, 13, 13]);

        // Mul
        let prod = &p1 * &p2;
        assert_eq!(prod.coeffs, vec![4, 13, 12, 11, 2]);

        // Eval
        assert_eq!(p1.eval(2), (1 + 2 * 2 + 3 * 4) % m);
        assert_eq!(p1.eval_at_one(), 6);
    }

    // #[test]
    // fn test_rem() {
    //     let m = 16;
    //     let p1 = Polynomial {
    //         coeffs: vec![1, 2, 3, 4],
    //         modulus: m,
    //     };
    //     let p2 = Polynomial {
    //         coeffs: vec![1, 1],
    //         modulus: m,
    //     };
    //     let rem = &p1 % &p2;
    //     assert_eq!(rem.coeffs, vec![1]);
    // }

    #[test]
    fn test_random() {
        let m = 17;
        let p = Polynomial::random(m, 5, &mut thread_rng());
        assert!(p.coeffs.len() == 5);
        assert!(p.coeffs.iter().all(|&x| x < m));
    }

    #[test]
    fn test_randshift() {
        let m = 17;
        let target = 5;
        
        // Test with thread_rng helper
        let p = Polynomial::randshift_with_thread_rng(target, m, 3);
        assert_eq!(p.eval_at_one(), target);

        // Test with explicit RNG
        let mut rng = rand::thread_rng();
        let p2 = Polynomial::randshift(target, m, 3, &mut rng);
        assert_eq!(p2.eval_at_one(), target);
    }

    #[test]
    fn test_neg_and_sub() {
        let m = 17;
        let p1 = Polynomial {
            coeffs: vec![1, 2, 3],
            modulus: m,
        };
        let p2 = Polynomial {
            coeffs: vec![4, 5, 6],
            modulus: m,
        };

        let neg_p1 = -p1.clone();
        assert_eq!(neg_p1.coeffs, vec![16, 15, 14]);

        let diff = &p1 - &p2;
        assert_eq!(diff.coeffs, vec![14, 14, 14]);
    }

    #[test]
    #[should_panic(expected = "polynomial multiplication overflow")]
    fn test_mul_overflow() {
        let p1 = Polynomial {
            coeffs: vec![u128::MAX, 2],
            modulus: u128::MAX,
        };
        let p2 = Polynomial {
            coeffs: vec![2, 1],
            modulus: u128::MAX,
        };
        let _ = &p1 * &p2;
    }

    #[test]
    fn test_degree() {
        let m = 17;
        let cases = vec![
            (vec![1], 0),
            (vec![1, 0], 0),
            (vec![1, 2], 1),
            (vec![1, 0, 0], 0),
            (vec![1, 2, 0], 1),
            (vec![1, 2, 3], 2),
            (vec![0, 0, 0, 4], 3),
        ];

        for (coeffs, expected_degree) in cases {
            let p = Polynomial {
                coeffs,
                modulus: m,
            };
            assert_eq!(p.degree(), expected_degree);
        }
    }

    #[test]
    fn test_with_modulus() {
        let p = Polynomial {
            coeffs: vec![15, 20, 25],
            modulus: 17,
        };

        // Change to new modulus
        let p2 = p.with_modulus(Some(7));
        assert_eq!(p2.coeffs, vec![1, 6, 4]);
        assert_eq!(p2.modulus, 7);

        // Remove modulus
        let p3 = p.with_modulus(None);
        assert_eq!(p3.coeffs, vec![15, 20, 25]);
        assert!(p3.modulus == p.modulus);
    }
}
