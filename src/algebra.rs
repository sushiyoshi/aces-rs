//! Homomorphic add / mult.  λ‑tensor is given externally.

use crate::{cipher::Cipher, polynomial::Polynomial};

/// Holds tensor & parameters.
pub struct AcesAlgebra {
    vanmod: u128,
    tensor: Vec<Vec<Vec<u128>>>, // λ_ijk  (dim³)
    dim: usize,
    u: Polynomial,
}

impl AcesAlgebra {
    /// Create new algebra with tensor from ArithChannel
    pub fn new(chan: &crate::arith_channel::ArithChannel, secret_key: &Vec<Polynomial>) -> Self {
        // Generate μ list and tensor for multiplication
        let mu_list = chan.generate_mu_list(secret_key);
        let tensor = chan.generate_tensor(secret_key, &mu_list);
        
        Self {
            vanmod: chan.p,
            tensor,
            dim: chan.dim,
            u: chan.u.clone(),
        }
    }

    pub fn add(&self, a: &Cipher, b: &Cipher) -> Cipher {
        let dec: Vec<Polynomial> = (0..self.dim)
            .map(|k| &(&a.dec[k] + &b.dec[k]) % &self.u)
            .collect();
        let enc = &(&a.enc + &b.enc) % &self.u;
        let level = a.level.checked_add(b.level)
            .expect("Arithmetic overflow occurred in add operation while adding noise levels");
        Cipher {
            dec,
            enc,
            level,
        }
    }

    pub fn mult(&self, a: &Cipher, b: &Cipher) -> Cipher {
        // t_k = Σ_i Σ_j λ_ijk · c0_i · d0_j
        let mut t = vec![Polynomial::zero(self.u.modulus); self.dim];
        for k in 0..self.dim {
            for i in 0..self.dim {
                for j in 0..self.dim {
                    let coeff = self.tensor[i][j][k];
                    let coeff_poly = Polynomial {
                        coeffs: vec![coeff],
                        modulus: self.u.modulus,
                    };
                    let term = &(&(&coeff_poly * &a.dec[i]) % &self.u * &b.dec[j]) % &self.u;
                    t[k] = &(&t[k] + term) % &self.u;
                }
            }
        }
        let dec: Vec<Polynomial> = (0..self.dim)
            .map(|k| {
                let term1 = &(&b.enc * &a.dec[k]) % &self.u;
                let term2 = &(&a.enc * &b.dec[k]) % &self.u;
                let sum = &(&term1 + &term2) % &self.u;
                &(sum - &t[k]) % &self.u
            })
            .collect();
        let enc = &(&a.enc * &b.enc) % &self.u;
        let level_sum = a.level.checked_add(b.level)
            .expect("Arithmetic overflow occurred in mult operation while adding noise levels");
        let level_prod = a.level.checked_mul(b.level)
            .expect("Arithmetic overflow occurred in mult operation while multiplying noise levels");
        let intermediate_level = level_sum.checked_add(level_prod)
            .expect("Arithmetic overflow occurred in mult operation while combining noise levels");
        let final_level = intermediate_level.checked_mul(self.vanmod as u128)
            .expect("Arithmetic overflow occurred in mult operation while scaling noise level");
        // println!("level_sum, level_prod, final_level: {}, {}, {}", level_sum, level_prod, final_level);
        Cipher {
            dec,
            enc,
            level: final_level,
        }
    }
}
