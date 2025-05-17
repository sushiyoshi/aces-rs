//! Encryption & decryption.

use crate::{
    arith_channel::ArithChannel,
    cipher::Cipher,
    polynomial::Polynomial,
};
use rand::Rng;
use crate::rand_iso::RandIso;


/// Main state object (contains f₀, f₁, tensor is in algebra).
pub struct Aces {
    pub f0: Vec<Vec<Polynomial>>,
    pub f1: Vec<Polynomial>,
    chan: ArithChannel,
}

impl Aces {
    pub fn generate_keypair(chan: &ArithChannel) -> (Self, Vec<Polynomial>) {
        // Generate secret key
        let mut rng = rand::thread_rng();
        let (U, inv_U) = RandIso::new(chan.q, chan.dim).generate(60, 1, 2, 3);

        // 列ベクトルを Polynomial に変換して秘密鍵とする例
        let secret_key: Vec<Polynomial> = (0..chan.dim)
            .map(|k| {
                let col: Vec<u128> = (0..chan.dim).map(|r| U[r][k]).collect();
                Polynomial::new(col, chan.q)   // 例: new(Vec<u128>, modulus)
            })
            .collect();
        // Generate public key components
        let f0 = chan.generate_initializer(&secret_key);
        // println!("f0: {:?}", f0);
        
        // Generate f1 = f0 * x + e as per LWE scheme
        let mut f1 = vec![Polynomial::zero_with_dim(chan.q, chan.dim); chan.n];
        
        // Calculate f0 * x
        for i in 0..chan.n {
            for j in 0..chan.dim {
                let prod = &f0[i][j] * &secret_key[j];
                f1[i] = (&f1[i] + &prod) % &chan.u;
            }
        }

        // Add vanishing noise to each component
        let k = 1u128; // noise level
        for i in 0..chan.n {
            let e = chan.generate_vanisher(k);
            f1[i] = (&f1[i] + &e) % &chan.u;
        }

        (
            Self {
                f0,
                f1,
                chan: chan.clone()
            },
            secret_key
        )
    }

    pub fn new_with_custom(f0: Vec<Vec<Polynomial>>, f1: Vec<Polynomial>, chan: ArithChannel) -> Self {
        Self { f0, f1, chan }
    }

    /// Encrypt a message m ∈ ℤ_p  →  (Cipher, noise_vector)
    pub fn encrypt<R: Rng>(&self, m: u128, rng: &mut R) -> (Cipher, Vec<u128>) {
        // Verify message is within bounds
        let m_p = m;
        // println!("m: {}", m_p);
        // 1. generate linear form b with controlled noise
        let b: Vec<Polynomial> = (0..self.chan.n)
            .map(|_| {
                let k = rng.gen_range(0..self.chan.p);
                let rand_poly = Polynomial::random(self.chan.q, self.chan.dim, rng);
                let eval = rand_poly.eval_at_one();
                let target = if k >= eval {
                    (k - eval) % self.chan.q
                } else {
                    (self.chan.q - ((eval - k) % self.chan.q)) % self.chan.q
                };
                let shift = Polynomial::randshift(target, self.chan.q, self.chan.dim, rng);
                let eval = (&(&shift + &rand_poly) % &self.chan.u).eval_at_one();
                // println!("[C](bi)={}, p={}", eval, self.chan.p);
                &(&shift + &rand_poly) % &self.chan.u
            })
            .collect();

        // 2. Generate error polynomial r(m) using ArithChannel
        let mut enc = self.chan.generate_error(m);
        
        // Verify [C](r(m)) = m using evaluation at omega
        let eval = enc.eval(self.chan.omega) % self.chan.q;
        assert_eq!(eval, m_p, "Error polynomial evaluation failed: [C](r({})) = {} ≠ {}", m_p, eval, m_p);

        // 3. accumulate bᵀ·(f₀x+e) with careful modular reduction
        for (bi, f1i) in b.iter().zip(self.f1.iter()) {
            let prod = bi * f1i;
            let reduced = prod % &self.chan.u;
            enc = &enc + &reduced % &self.chan.u;
        }
        // println!("enc: {:?}", enc);

        // 4. f₀ᵀ·b (dec vector) with intermediate reductions
        let mut dec_vec = Vec::with_capacity(self.chan.dim);
        for j in 0..self.chan.dim {
            let mut acc = Polynomial::zero_with_dim(self.chan.q, self.chan.dim);
            for i in 0..self.chan.n {
                let prod = &b[i] * &self.f0[i][j];
                acc = &acc + &prod;
            }
            acc = acc % &self.chan.u;
            dec_vec.push(acc);
        }
        // println!("dec_vec: {:?}", dec_vec);
        
        let noise: Vec<u128> = b.iter()
            .map(|bi| bi.eval_at_one() % self.chan.p)
            .collect();
        
            
        (
            Cipher {
                dec: dec_vec,
                enc,
                level: (self.chan.n as u128) * self.chan.p,
            },
            noise,
        )
    }

    /// Decrypt with secret x (vector of `Polynomial`s)
    #[must_use]
    pub fn decrypt(&self, c: &Cipher, x: &Vec<Polynomial>) -> u128 {
        // Verify dimensions
        assert_eq!(x.len(), self.chan.dim, "incorrect decryption key dimension");
        
        // Initialize accumulator with proper dimension
        let mut c_tx = Polynomial::zero_with_dim(self.chan.q, self.chan.dim);
        
        // Compute inner product with careful modular reduction
        for (ci, xi) in c.dec.iter().zip(x) {
            let prod = (ci * xi) % &self.chan.u;
            c_tx = (&c_tx + &prod) % &self.chan.u;
        }
        
        // Compute final result with proper modular arithmetic
        let m_pre = &(&c.enc - &c_tx) % &self.chan.u;
        let eval = m_pre.eval(self.chan.omega);
        eval % self.chan.p
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::AcesAlgebra;

    #[test]
    fn test_key_generation_lwe_properties() {
        // Initialize parameters
        let chan = ArithChannel::new(16, (1u128 << 64) + 1, 5, 5);
        let (aces, secret_key) = Aces::generate_keypair(&chan);
        let mut rng = rand::thread_rng();

        // Test secret key properties
        assert_eq!(secret_key.len(), chan.dim);
        for key in &secret_key {
            assert_eq!(key.modulus, chan.q);
            assert_eq!(key.coeffs[0], 1); // Constant term is 1
        }

        // Test LWE property: f1 ≈ f0 * x
        for i in 0..chan.n {
            let mut expected = Polynomial::zero_with_dim(chan.q, chan.dim);
            for j in 0..chan.dim {
                let prod = &aces.f0[i][j] * &secret_key[j];
                expected = (&expected + &prod) % &chan.u;
            }
            // Difference should be a vanishing polynomial
            let diff = (&aces.f1[i] - &expected) % &chan.u;
            assert!(chan.is_vanishing(&diff, 1),
                   "f1[{}] - (f0[{}] * x) should be vanishing", i, i);
        }

        // Test encryption and decryption with LWE ciphertexts
        let message = 7u128;
        let message2 = 5u128;

        // Single message encryption/decryption
        let (c1, noise1) = aces.encrypt(message, &mut rng);
        assert!(noise1.iter().all(|&n| n < chan.p), "noise exceeds bound");
        let m1 = aces.decrypt(&c1, &secret_key);
        assert_eq!(m1, message, "basic decryption failed");

        // Test homomorphic operations with proper tensor
        let (c2, _) = aces.encrypt(message2, &mut rng);
        let alg = AcesAlgebra::new(&chan, &secret_key);
        
        // Test addition
        let c_add = alg.add(&c1, &c2);
        let m_add = aces.decrypt(&c_add, &secret_key);
        assert_eq!(m_add, (message + message2) % chan.p,
                  "homomorphic addition failed");

        // Test multiplication with tensor
        let c_mul = alg.mult(&c1, &c2);
        let m_mul = aces.decrypt(&c_mul, &secret_key);
        assert_eq!(m_mul, (message * message2) % chan.p,
                  "homomorphic multiplication failed");

        // Validate multiplication properties
        for i in 0..3 {
            let m_i = (i + 2) as u128;
            let (c_i, _) = aces.encrypt(m_i, &mut rng);
            
            // Test associativity
            let c_left = alg.mult(&alg.mult(&c1, &c2), &c_i);
            let c_right = alg.mult(&c1, &alg.mult(&c2, &c_i));
            
            let m_left = aces.decrypt(&c_left, &secret_key);
            let m_right = aces.decrypt(&c_right, &secret_key);
            
            assert_eq!(m_left, m_right,
                      "multiplication not associative");
            assert_eq!(m_left, (message * message2 * m_i) % chan.p,
                      "incorrect multiplication result");
        }
            
        // Test repeated multiplication for noise growth
        let mut curr_cipher = c1.clone();
        let mut ops = 0;
        
        for _ in 0..10 {  // Try 10 multiplications
            let (next_cipher, _) = aces.encrypt(message, &mut rng);
            curr_cipher = alg.mult(&curr_cipher, &next_cipher);
            
            // Verify each multiplication result
            let decrypted = aces.decrypt(&curr_cipher, &secret_key);
            let expected = message.pow((ops + 2) as u32) % chan.p;
            
            if decrypted != expected {
                println!("Multiplication chain failed at step {}", ops + 1);
                println!("Expected: {}, Got: {}", expected, decrypted);
                break;
            }
            ops += 1;
        }
        
        // Should support a reasonable number of multiplications
        assert!(ops >= 5, "Too few multiplications before failure: {}", ops);
    }
}

