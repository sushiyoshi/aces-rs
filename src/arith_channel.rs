//! Arithmetic Channel implementation for ACES.

use crate::polynomial::Polynomial;
use rand::{Rng, thread_rng};

/// Arithmetic channel as in ACES §5.1: (p, q, ω, u).
#[derive(Clone, Debug)]
pub struct ArithChannel {
    /// Vanishing modulus p (smaller prime)
    pub p: u128,
    /// Integer modulus q (larger prime)
    pub q: u128,
    /// Dimension for polynomials
    pub dim: usize,
    /// Number of components
    pub n: usize,
    /// Modulus polynomial u
    pub u: Polynomial,
    /// Evaluation point ω
    pub omega: u128,
    
    /// Prime factors of q
    q_factors: Vec<u128>,
    /// n-repartition mapping
    n_repartition: Vec<usize>,
}

impl ArithChannel {
    /// Create new arithmetic channel with p < q and gcd(p,q) = 1
    pub fn new(p: u128, q: u128, dim: usize, n: usize) -> Self {
        assert!(p < q, "p must be smaller than q");
        assert!(gcd(p, q) == 1, "p and q must be coprime");

        // Generate modulus polynomial u with random coefficients
        let u = Self::generate_u(dim, q);

        // Factor q into primes
        let q_factors = factor_intmod(q);
        
        // Generate random n-repartition
        let n_repartition = Self::generate_n_repartition(n, q_factors.len());


        Self {
            p,
            q,
            dim,
            n,
            u,
            omega: 1,
            q_factors,
            n_repartition,
        }
    }

    /// Evaluate polynomial at x = omega
    pub fn eval(&self, poly: &Polynomial) -> u128 {
        assert_eq!(poly.modulus, self.q, "polynomial modulus mismatch");
        poly.eval(self.omega)
    }

    /// Get sigma[q]_{i,j} value
    pub fn sigma_q(&self, i: usize, j: usize) -> u128 {
        let sigma_i = self.n_repartition[i];
        let sigma_j = self.n_repartition[j];
        
        if sigma_i == sigma_j {
            self.q / self.q_factors[sigma_i]
        } else {
            self.q / (self.q_factors[sigma_i] * self.q_factors[sigma_j])
        }
    }

    /// Generate random n-repartition mapping
    fn generate_n_repartition(n: usize, factor_count: usize) -> Vec<usize> {
        let mut rng = rand::thread_rng();
        let mut ret = vec![0; n];
        
        // Ensure each factor has at least one assignment
        let indices: Vec<usize> = (0..n).collect();
        for i in 0..factor_count {
            if let Some(pos) = indices.get(i) {
                ret[*pos] = i;
            }
        }

        // Randomly assign remaining positions
        for i in factor_count..n {
            if let Some(pos) = indices.get(i) {
                ret[*pos] = rng.gen_range(0..factor_count);
            }
        }

        ret
    }

    /// Calculate Bezout coefficients for list of integers using safe modular arithmetic
    pub fn bezout_coefficients(values: &[u128], modulus: u128) -> Vec<u128> {
        let n = values.len();
        if n == 0 {
            return vec![];
        }

        // Initialize with normalized first value
        let mut g = values[0] % modulus;
        let mut x_list = vec![1u128];
        x_list.extend(vec![0u128; n-1]);

        for i in 1..n {
            let v_i = values[i] % modulus;
            let (g_val, m, n_) = extended_gcd2(g, v_i);
            
            // Convert signed coefficients to unsigned modular values
            let m_pos = m.unsigned_abs() as u128;
            let m_unsigned = if m < 0 {
                modulus - (m_pos % modulus)
            } else {
                m_pos % modulus
            };

            let n_pos = n_.unsigned_abs() as u128;
            let n_unsigned = if n_ < 0 {
                modulus - (n_pos % modulus)
            } else {
                n_pos % modulus
            };

            // Update coefficients using safe modular arithmetic
            let mut new_coeffs = vec![0u128; i];
            for j in 0..i {
                new_coeffs[j] = Self::mul(x_list[j], m_unsigned, modulus);
            }
            // Copy back updated coefficients
            for j in 0..i {
                x_list[j] = new_coeffs[j];
            }
            x_list[i] = n_unsigned;
            
            g = g_val % modulus;
        }

        // Final normalization with overflow protection
        for x in x_list.iter_mut() {
            *x = *x % modulus;
        }
        x_list
    }
    /// Safe modular multiplication
    fn mul(a: u128, b: u128, modulus: u128) -> u128 {
        let a = a % modulus;
        let b = b % modulus;
        let result = a.checked_mul(b)
            .expect("Arithmetic overflow occurred in ArithChannel::mul");
        result % modulus
    }

    /// Generate μ list for the given secret key
    pub fn generate_mu_list(&self, secret: &Vec<Polynomial>) -> Vec<u128> {
        let mut divisor_nums = Vec::with_capacity(self.n);
        for i in 0..self.n {
            let q_sigma_i = self.q_factors[self.n_repartition[i]];
            let eval_xi = secret[i].eval(self.omega) % self.q;
            // Safe modular multiplication
            divisor_nums.push(Self::mul(q_sigma_i, eval_xi, self.q));
        }
        
        Self::bezout_coefficients(&divisor_nums, self.q)
    }


    /// Generate initializer f0 such that each row of f0 is a multiple of (at least one) prime factor q_factor
    // pub fn generate_initializer(&self) -> Vec<Vec<Polynomial>> {
    //     let mut f0 = Vec::with_capacity(self.n);
    //     let mut rng = thread_rng();

    //     for _ in 0..self.n {
    //         let mut row = Vec::with_capacity(self.dim);
    //         for j in 0..self.dim {
    //             // Generate random polynomial and evaluate at omega
    //             let k = rng.gen_range(0..self.q);
    //             let rand_poly = Polynomial::random(self.q, self.dim, &mut rng);
    //             let eval = rand_poly.eval(self.omega);
                
    //             // Calculate shift value
    //             let target = (self.p.checked_mul(k)
    //                 .expect("Overflow in p*k") % self.q)
    //                 .checked_sub(eval % self.q)
    //                 .unwrap_or(self.q - (eval % self.q));
                
    //             // Generate shift polynomial and combine with random poly
    //             let shift = Polynomial::randshift(target, self.q, self.dim, &mut rng);
    //             let combined = &(&shift + &rand_poly) % &self.u;
                
    //             // Scale by q_sigma_j
    //             let q_sigma_j = self.q_factors[self.n_repartition[j]];
    //             let scaled = &(&Polynomial::constant(q_sigma_j, self.q) * combined) % &self.u;
                
    //             row.push(scaled);
    //         }
    //         f0.push(row);
    //     }
    //     f0
    // }

    /// Build f₀ exactly like the python reference (no secret mixing).
    pub fn generate_initializer(&self, _secret: &Vec<Polynomial>) -> Vec<Vec<Polynomial>> {
        use rand::Rng;
        let mut rng = thread_rng();
        let mut f0 = Vec::with_capacity(self.n);

        for _ in 0..self.n {
            let mut row = Vec::with_capacity(self.dim);
            for j in 0..self.dim {
                // choose k ∈ [0, q)
                let k = rng.gen_range(0..self.q);

                // fresh random polynomial
                let rand_poly = Polynomial::random(self.q, self.dim, &mut rng);
                let eval = rand_poly.eval_at_one() % self.q;

                // shift so constant term becomes  p·k  (mod q)
                let target = (self.p * k + self.q - eval) % self.q;
                let shift = Polynomial::randshift(target, self.q, self.dim, &mut rng);

                // base polynomial before scaling
                let base = &shift + &rand_poly;   // degree < dim
                // multiply by q_sigma_j  and reduce mod u
                let q_sigma_j = self.q_factors[self.n_repartition[j]];
                let scaled =
                    &Polynomial::constant(q_sigma_j, self.q) * &base;

                row.push(scaled);
            }
            f0.push(row);
        }
        f0
    }
    /// Generate tensor for multiplication
    pub fn generate_tensor(
        &self,
        secret: &Vec<Polynomial>,
        mu_list: &[u128]
    ) -> Vec<Vec<Vec<u128>>> {
        let mut tensor = vec![vec![vec![0; self.n]; self.n]; self.n];
        
        for i in 0..self.n {
            for j in 0..self.n {
                // Calculate x_i * x_j mod u
                let prod = &secret[i] * &secret[j];
                let xi_xj = &prod % &self.u;
                let eval_xixj = xi_xj.eval(self.omega);
                let l_ij = rand::thread_rng().gen_range(0..self.q);
                let sigma_q_ij = self.sigma_q(i, j);

                // Calculate tensor coefficients with careful modular arithmetic
                for k in 0..self.n {
                    let q_sigma_k = self.q_factors[self.n_repartition[k]];
                    let mu_k = mu_list[k];
                    
                    // Compute difference term with safe arithmetic
                    let l_term = (l_ij * sigma_q_ij) % self.q;
                    let diff = if eval_xixj >= l_term {
                        (eval_xixj - l_term) % self.q
                    } else {
                        (self.q + eval_xixj - l_term) % self.q
                    };
                    
                    // Combine terms with safe multiplication
                    let term1 = (q_sigma_k * mu_k) % self.q;
                    tensor[i][j][k] = (term1 * diff) % self.q;
                }
            }
        }

        tensor
    }

    /// Generate error polynomial r(m) following Example 5.26
    /// r(m) = ((m*ω^(-s_m) - Σ(a_m,j*ω^(j-s_m))) mod q)*X^s_m + Σ(a_m,j*X^j)
    pub fn generate_error(&self, m: u128) -> Polynomial {
        let mut rng = rand::thread_rng();
        
        // Choose random s_m in [0, dim-1]
        let s_m = rng.gen_range(0..self.dim);
        
        // Generate random coefficients a_m,j
        let mut a_m = vec![0u128; self.dim];
        for j in 0..self.dim {
            if j != s_m {
                a_m[j] = rng.gen_range(0..self.q);
            }
        }
        
        // Calculate ω^(-s_m) using modular inverse
        let omega_inv = mod_inverse(self.omega, self.q)
            .expect("omega must be invertible");
        let omega_neg_sm = mod_pow(omega_inv, s_m as u128, self.q);
        
        // Calculate Σ(a_m,j*ω^(j-s_m))
        let mut sum = 0u128;
        for j in 0..self.dim {
            if j != s_m {
                let omega_j_sm = if j >= s_m {
                    mod_pow(self.omega, (j - s_m) as u128, self.q)
                } else {
                    mod_pow(omega_inv, (s_m - j) as u128, self.q)
                };
                sum = (sum + a_m[j] * omega_j_sm) % self.q;
            }
        }
        
        // Calculate coefficient of X^s_m
        let coeff_sm = ((m * omega_neg_sm) % self.q + self.q - sum) % self.q;
        
        // Construct the polynomial
        let mut coeffs = a_m;
        coeffs[s_m] = coeff_sm;
        
        Polynomial { coeffs, modulus: self.q }
    }

    /// Generate vanishing noise polynomial e with [C](e) = k*p
    pub fn generate_vanisher(&self, k: u128) -> Polynomial {
        let mut rng = rand::thread_rng();
        let rand_poly = Polynomial::random(self.q, self.dim, &mut rng);
        let eval = rand_poly.eval(self.omega);
        
        // Calculate shift to make [C](e) = k*p
        let target = ((k * self.p) % self.q + self.q - eval) % self.q;
        let shift = Polynomial::randshift(target, self.q, self.dim, &mut rng);
        
        &shift + &rand_poly
    }

    /// Check if polynomial e is in I_k(C)
    pub fn is_vanishing(&self, e: &Polynomial, k: u128) -> bool {
        let eval = e.eval(self.omega);
        
        // First check if evaluation is multiple of p
        if eval % self.p != 0 {
            return false;
        }

        // Ensure the evaluation is within bounds
        let max_val = k.checked_mul(self.p)
            .expect("overflow in vanishing polynomial check");
        eval <= max_val
    }
    /// Generate monic polynomial u of degree `dim` such that
    ///   (1) u(1) = 0   (so ω=1 is a root)
    ///   (2) leading coeff = 1  (monic)
    fn generate_u(dim: usize, q: u128) -> Polynomial {
        use rand::Rng;
        use rand_distr::{Normal, Distribution};
        let mut rng  = rand::thread_rng();

        // ---- 1. ランダムに dim-1 個の係数を作る（定数項は後で調整） ----
        //    u(x) = x^dim + a_{dim-1}x^{dim-1} + ... + a_1 x + a_0
        let mean = (3.0 * dim as f64) / 4.0;
        let std  =  dim as f64 / 4.0;
        let normal = Normal::new(mean, std).unwrap();

        // 個々の係数を出来るだけ散らす
        let mut coeffs: Vec<u128> = Vec::with_capacity(dim + 1);
        coeffs.push(1);                // x^dim の係数 = 1  (monic)

        for _ in 1..dim {
            let c  = rng.gen_range(0..q);
            coeffs.push(c);
        }

        // ---- 2. 定数項を決定して Σ coeffs ≡ 0 (mod q) を満たす ----
        let sum_except_c0: u128 =
            coeffs.iter().fold(0u128, |acc, &x| (acc + x) % q);
        let c0 = (q - sum_except_c0) % q;   // (1) を満たすための定数項
        coeffs.push(c0);
        Polynomial { coeffs: coeffs, modulus: q }
    }
}

/// Extended Euclidean Algorithm with safe handling of negative values
pub fn extended_gcd2(a: u128, b: u128) -> (u128, i128, i128) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x1, y1) = extended_gcd2(b, a % b);
        let x = y1;
        let q = (a / b) as i128;
        let y = x1 - q * y1;
        (g, x, y)
    }
}

/// Calculate GCD using Euclidean algorithm
fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Calculate modular multiplicative inverse
fn mod_inverse(a: u128, m: u128) -> Option<u128> {
    let (g, x, _) = extended_gcd2(a, m);
    if g != 1 {
        None
    } else {
        Some(((x % m as i128 + m as i128) as u128) % m)
    }
}

/// Calculate modular exponentiation
fn mod_pow(base: u128, exponent: u128, modulus: u128) -> u128 {
    if exponent == 0 {
        return 1;
    }
    
    let mut base = base % modulus;
    let mut exponent = exponent;
    let mut result = 1u128;
    
    while exponent > 0 {
        if exponent % 2 == 1 {
            result = (result * base) % modulus;
        }
        base = (base * base) % modulus;
        exponent /= 2;
    }
    
    result
}

/// Factor integer into primes
fn factor_intmod(mut n: u128) -> Vec<u128> {
    let mut factors = Vec::new();
    let mut d = 2;
    
    while d * d <= n {
        while n % d == 0 {
            factors.push(d);
            n /= d;
        }
        d += if d == 2 { 1 } else { 2 };
    }
    
    if n > 1 {
        factors.push(n);
    }
    
    factors
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_channel_creation() {
        let p = 7;
        let q = 31;
        let dim = 4;
        let n = 3;
        
        let chan = ArithChannel::new(p, q, dim, n);
        
        assert_eq!(chan.p, p);
        assert_eq!(chan.q, q);
        assert_eq!(chan.dim, dim);
        assert_eq!(chan.n, n);
        assert_eq!(chan.omega, 1);
        
        // Check modulus polynomial properties
        assert_eq!(chan.u.coeffs.len(), dim + 1);
        assert_eq!(chan.u.modulus, q);
        // Highest coefficient should be 1
        assert_eq!(chan.u.coeffs[dim], 1);
    }

    #[test]
    fn test_factor_intmod() {
        assert_eq!(factor_intmod(12), vec![2, 2, 3]);
        assert_eq!(factor_intmod(17), vec![17]);
        assert_eq!(factor_intmod(100), vec![2, 2, 5, 5]);
    }

    #[test]
    fn test_n_repartition() {
        let n = 5;
        let factor_count = 3;
        
        let rep = ArithChannel::generate_n_repartition(n, factor_count);
        
        assert_eq!(rep.len(), n);
        assert!(rep.iter().all(|&x| x < factor_count));
        
        // Check all factors are used at least once
        let mut used = vec![false; factor_count];
        for &x in &rep {
            used[x] = true;
        }
        assert!(used.iter().all(|&x| x));
    }

    #[test]
    fn test_sigma_q() {
        let chan = ArithChannel::new(7, 15, 4, 3); // 15 = 3 * 5
        
        // Same sigma values
        let same_sigma = chan.sigma_q(0, 0);
        assert_eq!(same_sigma * chan.q_factors[chan.n_repartition[0]], chan.q);
        
        // Different sigma values
        let diff_sigma = chan.sigma_q(0, 1);
        assert_eq!(
            diff_sigma * chan.q_factors[chan.n_repartition[0]] * chan.q_factors[chan.n_repartition[1]],
            chan.q
        );
    }

    /// Test error polynomial generation and evaluation
    #[test]
    fn test_error_polynomial() {
        let chan = ArithChannel::new(7, 31, 4, 3);
        
        // Test multiple messages
        for m in 0..chan.p {
            let r = chan.generate_error(m);
            
            // Check [C](r(m)) = m
            let eval = r.eval(chan.omega) % chan.q;
            assert_eq!(eval, m, "Error polynomial evaluation failed for m = {}", m);
            
            // Check polynomial degree is within bounds
            assert!(r.coeffs.len() <= chan.dim,
                   "Error polynomial degree too large for m = {}", m);
            
            // Check coefficients are in Zq
            assert!(r.coeffs.iter().all(|&c| c < chan.q),
                   "Coefficients not reduced modulo q for m = {}", m);
        }
        
        // Test with omega's inverse exists
        assert!(mod_inverse(chan.omega, chan.q).is_some(),
               "omega must be invertible");
    }
    
    #[test]
    fn test_mod_inverse() {
        // Test known values
        assert_eq!(mod_inverse(3, 11), Some(4)); // 3 * 4 = 12 ≡ 1 (mod 11)
        assert_eq!(mod_inverse(5, 10), None); // 5 and 10 are not coprime
    }
    
    #[test]
    fn test_mod_pow() {
        assert_eq!(mod_pow(2, 5, 13), 6);    // 2^5 = 32 ≡ 6 (mod 13)
        assert_eq!(mod_pow(3, 7, 11), 9);    // 3^7 = 2187 ≡ 9 (mod 11)
    }

    #[test]
    fn test_vanishing_polynomial() {
        let chan = ArithChannel::new(7, 31, 4, 3);
        let k = 2;
        
        let e = chan.generate_vanisher(k);
        assert!(chan.is_vanishing(&e, k));
        
        // Test non-vanishing polynomial
        let other = Polynomial::random(chan.q, chan.dim, &mut rand::thread_rng());
        assert!(!chan.is_vanishing(&other, k));
    }

    #[test]
    #[should_panic(expected = "p must be smaller than q")]
    fn test_invalid_p_q_order() {
        ArithChannel::new(31, 7, 4, 3);
    }

    #[test]
    #[should_panic(expected = "p and q must be coprime")]
    fn test_invalid_p_q_gcd() {
        ArithChannel::new(15, 45, 4, 3);
    }

    #[test]
    fn test_bezout_coefficients() {
        let values = vec![2, 3, 5];
        let modulus = 16;
        let coeffs = ArithChannel::bezout_coefficients(&values, modulus);
        
        // Verify that sum(a_i * x_i) ≡ 1 (mod q)
        let sum = values.iter().zip(coeffs.iter())
            .fold(0u128, |acc, (&v, &c)|
                (acc + (v * c) % modulus) % modulus
            );
        assert_eq!(sum, 1);
    }

    #[test]
    fn test_mu_list() {
        let chan = ArithChannel::new(16, (1u128 << 64) + 1, 5, 5);
        let mut rng = rand::thread_rng();
        
        // Create some secret polynomials
        let secret: Vec<_> = (0..chan.n)
            .map(|_| Polynomial::random(chan.q, chan.dim, &mut rng))
            .collect();
            
        let mu_list = chan.generate_mu_list(&secret);
        
        // Verify length
        assert_eq!(mu_list.len(), chan.n);
        
        // Verify Bezout equation: sum(q_sigma_i * [C](x_i) * mu_i) ≡ 1 (mod q)
        let sum = (0..chan.n).fold(0u128, |acc, i| {
            let q_sigma_i = chan.q_factors[chan.n_repartition[i]];
            let eval_xi = secret[i].eval(chan.omega);
            (acc + (q_sigma_i * eval_xi % chan.q) * mu_list[i]) % chan.q
        });
        assert_eq!(sum, 1);
    }

    #[test]
    fn test_tensor_properties() {
        let chan = ArithChannel::new(16, (1u128 << 64) + 1, 5, 5);
        let mut rng = rand::thread_rng();
        
        // Create secret and generate mu_list
        let secret: Vec<_> = (0..chan.n)
            .map(|_| Polynomial::random(chan.q, chan.dim, &mut rng))
            .collect();
        let mu_list = chan.generate_mu_list(&secret);
        
        // Generate tensor
        let tensor = chan.generate_tensor(&secret, &mu_list);
        
        // Check tensor dimensions
        assert_eq!(tensor.len(), chan.n);
        for i in 0..chan.n {
            assert_eq!(tensor[i].len(), chan.n);
            for j in 0..chan.n {
                assert_eq!(tensor[i][j].len(), chan.n);
            }
        }

        // Verify tensor elements are multiples of q_sigma_k
        for i in 0..chan.n {
            for j in 0..chan.n {
                for k in 0..chan.n {
                    let q_sigma_k = chan.q_factors[chan.n_repartition[k]];
                    assert_eq!(tensor[i][j][k] % q_sigma_k, 0);
                }
            }
        }

        // Optional: Verify tensor equation for x_i * x_j
        for i in 0..chan.n {
            for j in 0..chan.n {
                // Calculate x_i * x_j mod u
                let prod = &(&secret[i] * &secret[j]);
                let xi_xj = &(prod % &chan.u);

                // Calculate sum of lambda_ijk * x_k
                let mut sum = Polynomial::zero(chan.q);
                for k in 0..chan.n {
                    let lambda = tensor[i][j][k];
                    let rand_term = &Polynomial::randshift_with_thread_rng(lambda, chan.q, chan.dim);
                    let term = &(rand_term * &secret[k]);
                    let term_mod = &(term % &chan.u);
                    sum = sum + term_mod;
                }

                // Final reduction modulo u
                let sum_reduced = &sum % &chan.u;
                
                // Check that difference is a multiple of sigma_q_ij
                let sigma_q_ij = chan.sigma_q(i, j);
                let diff = &(xi_xj - sum_reduced);
                assert!(diff.eval_at_one() % sigma_q_ij == 0,
                    "Tensor equation failed for i = {}, j = {}", i, j);
            }
        }
    }
}

