//! Cheap level‑reset (Refresh) operation.
use core::panic;

use rand::{Rng, thread_rng};
use crate::{
    algebra::AcesAlgebra,
    arith_channel::ArithChannel,
    cipher::Cipher,
    polynomial::Polynomial,
    scheme::Aces,
};

/// Stateless helper that owns references to ACES scheme + algebra.
pub struct Refresher<'a> {
    scheme: &'a Aces,
    alg: &'a AcesAlgebra,
    chan: &'a ArithChannel,
}

impl<'a> Refresher<'a> {
    pub fn new(scheme: &'a Aces, alg: &'a AcesAlgebra, chan: &'a ArithChannel) -> Self {
        Self { scheme, alg, chan }
    }
    pub fn is_refreshable(&self, c: &Cipher, secret_x: &Vec<Polynomial>) -> bool {
        let ps_c0: Vec<u128> = c.dec
            .iter()
            .map(|ci| {
                let eval = ci.eval(self.chan.omega) % self.chan.q;
                let (diff, overflow1) = self.chan.q.overflowing_sub(eval);
                if overflow1 {
                    println!("Warning: Underflow detected at line 30: q ({}) - eval ({})", self.chan.q, eval);
                }
                let (rem, overflow2) = diff.overflowing_rem(self.chan.q);
                if overflow2 {
                    println!("Warning: Overflow detected at line 30 in modulo operation");
                }
                rem
                // ((self.chan.q as i128 - eval as i128) % self.chan.q as i128) as u128
            })
            .collect();
        let ps_c1 = c.enc.eval(self.chan.omega) % self.chan.q;

        let x_evals: Vec<u128> = secret_x.iter()
            .map(|xi| {
                let xi_u = xi % &self.chan.u;
                xi_u.eval(self.chan.omega) % self.chan.q
            })
            .collect();

        let mut c0x_sum = 0u128;
        for (ps_c0_i, x_evals_i) in ps_c0.iter().zip(x_evals.iter()) {
            // println!("c0_x_sum: {}, ps_c0_i: {}, x_evals_i: {}", c0x_sum, ps_c0_i, x_evals_i);
            let (prod, overflow1) = ps_c0_i.overflowing_mul(*x_evals_i);
            if overflow1 {
                println!("Warning: Overflow detected at line 46 in multiplication: ps_c0_i ({}) * x_evals_i ({})", ps_c0_i, x_evals_i);
            }
            let (sum, overflow2) = c0x_sum.overflowing_add(prod);
            if overflow2 {
                println!("Warning: Overflow detected at line 46 in addition: c0x_sum ({}) + prod ({})", c0x_sum, prod);
            }
            c0x_sum = sum;
        }
        let (sum, overflow1) = ps_c1.overflowing_add(c0x_sum);
        if overflow1 {
            println!("Warning: Overflow detected at line 48: ps_c1 ({}) + c0x_sum ({})", ps_c1, c0x_sum);
        }
        let (rem1, overflow2) = sum.overflowing_rem(self.chan.q);
        if overflow2 {
            println!("Warning: Overflow detected at line 48 in first modulo operation");
        }
        let (right, overflow3) = rem1.overflowing_rem(self.chan.p);
        if overflow3 {
            println!("Warning: Overflow detected at line 48 in second modulo operation");
        }

        let (ps_c1_mod_p, overflow4) = ps_c1.overflowing_rem(self.chan.p);
        if overflow4 {
            println!("Warning: Overflow detected at line 49 in ps_c1 modulo operation");
        }
        let (c0x_sum_mod_p, overflow5) = c0x_sum.overflowing_rem(self.chan.p);
        if overflow5 {
            println!("Warning: Overflow detected at line 49 in c0x_sum modulo operation");
        }
        let (sum2, overflow6) = ps_c1_mod_p.overflowing_add(c0x_sum_mod_p);
        if overflow6 {
            println!("Warning: Overflow detected at line 49: ps_c1_mod_p ({}) + c0x_sum_mod_p ({})",
                ps_c1_mod_p, c0x_sum_mod_p);
        }
        let (left, overflow7) = sum2.overflowing_rem(self.chan.p);
        if overflow7 {
            println!("Warning: Overflow detected at line 49 in final modulo operation");
        }
        // println!("ps_c1, {}, c0x_sum: {},ps_c1 + c0x_sum: {}", ps_c1, c0x_sum, ps_c1 + c0x_sum);
        // println!("Left: {} Right: {}", left, right);
        let (_left2, overflow) = ps_c1.overflowing_add(c0x_sum);
        if overflow {
            println!("Warning: Overflow detected at line 103: ps_c1 ({}) + c0x_sum ({})", ps_c1, c0x_sum);
        }
        left == right
        // let mut count = 0;
        // if ps_c1 + c0x_sum > self.chan.q {
        //     let max_loop = 10;
        //     for i in 0..max_loop {
        //        if left2 < self.chan.q {
        //            break;
        //        }
        //        left2 = left2 - self.chan.q;
        //        count = count + 1;
        //     }
        // }
        // println!("Left: {}, left-count: {},right: {}", left, left2 % self.chan.p, right);
        // println!("count: {}", count);
        // left-count == right

    }
    // pub fn is_refreshable(&self, c: &Cipher, _secret_x: &Vec<Polynomial>) -> bool {
    //     use crate::locator::{LocatorDB, Vector};

    //     // --- 0) 基本レベル判定（Theorem 5.46 左辺）
    //     let level_bound = (self.chan.q + 1) / self.chan.p - 1;
    //     if c.level <= level_bound {
    //         return true;
    //     }

    //     // --- 1) 〚C〛((c)) = (-c₀)   を Z_q ベクトルで取得
    //     let ell: Vector = c
    //         .dec
    //         .iter()
    //         .map(|ci| {
    //             let eval = ci.eval(self.chan.omega) % self.chan.q;
    //             (self.chan.q + self.chan.q - eval) % self.chan.q // = −c₀ (mod q)
    //         })
    //         .collect();

    //     // --- 2) locator/director DB で分解を試みる
    //     let db = LocatorDB::new(self.chan.dim, self.chan.q);
    //     db.has_decomposition(&ell)
    // }

    /// Make a ciphertext refreshable by repeatedly adding encrypted zeros
    /// Will panic after MAX_ATTEMPTS (default: 10) unsuccessful attempts
    pub fn make_refreshable<R: Rng>(
        &self,
        cipher: &Cipher,
        secret_x: &Vec<Polynomial>,
        rng: &mut R,
    ) -> (Option<Cipher>, u32) {
        const MAX_ATTEMPTS: u32 = 10000;

        // First check if already refreshable
        if self.is_refreshable(cipher, secret_x) {
            return (Some(cipher.clone()), 0);
        }

        // Check initial level bound
        let level_bound = (self.chan.q + 1) / self.chan.p - 1;
        let mut current = cipher.clone();
        let mut attempts = 0;

        // Keep trying new zero ciphers until we get a refreshable result
        // while !c&current, secret_x) {
        loop {
            attempts += 1;
            if attempts > MAX_ATTEMPTS {
                panic!("Failed to make cipher refreshable after {} attempts", MAX_ATTEMPTS);
            }

            // Generate a fresh encryption of zero
            let (zero_cipher, _) = self.scheme.encrypt(0, rng);
            
            // Add it to our current cipher
            let next = self.alg.add(&current, &zero_cipher);
            // println!("nextlevel: {}, level_bound: {}", next.level, level_bound);
            // Check if level would exceed bound
            if next.level >= level_bound {
                println!("Level would exceed bound: {} >= {}",next.level, level_bound);
                panic!("Level exceeded bound after {} attempts", attempts);
            }
            
            // current = next;
            if self.is_refreshable(&next, secret_x) {
                current = next;
                break;
            }

            // println!("current ciphertext: {:?}", current);
            // println!("Attempt {}/{}: Added zero cipher", attempts, MAX_ATTEMPTS);
        }

        // println!("Successfully made refreshable after {} attempts", attempts);
        (Some(current),attempts)
    }
    /// Core refresh algorithm (ACES §5.5).
    #[must_use]
    pub fn refresh(&self, c: &Cipher, secret_x: &Vec<Polynomial>) -> Cipher {
        // 1. Encrypt every x_i (refresher vector ρ)
        let mut rng = thread_rng();

        // 1. Create refresher vector from secret key evaluations
        let (rho_cipher, rho_noise): (Vec<_>, Vec<_>) = secret_x
            .iter()
            .map(|xi| {
                let xi_u = xi % &self.chan.u;
                let eval = xi_u.eval(self.chan.omega) % self.chan.q % self.chan.p;
                self.scheme.encrypt(eval, &mut rng)
            })
            .unzip();

        // 2. Create negated c₀ components
        let c0_int: Vec<u128> = c.dec
            .iter()
            .map(|ci| {
                let ci_u = ci % &self.chan.u;
                let eval = ci_u.eval(self.chan.omega) % self.chan.q;
                (self.chan.q - eval) % self.chan.q
            })
            .collect();
        let (c0_enc, c0_noise): (Vec<_>, Vec<_>) = c0_int
            .iter()
            .map(|&m| self.scheme.encrypt(m % self.chan.p, &mut rng))
            .unzip();
            
        // 3. Get evaluation of c₁ = [C](c)
        let c1_eval = c.enc.eval(self.chan.omega) % self.chan.q;
        let c1_mod_p = c1_eval % self.chan.p;
        let (c1_enc, c1_noise) = self.scheme.encrypt(c1_mod_p, &mut rng);

        // 5. Scalar product ⟨c0_enc, ρ⟩
        let mut sp = self.alg.mult(&c0_enc[0], &rho_cipher[0]);
        for i in 1..rho_cipher.len() {
            sp = self.alg.add(&sp, &self.alg.mult(&c0_enc[i], &rho_cipher[i]));
        }
        
        // 6. Add c1_enc
        let result = self.alg.add(&c1_enc, &sp);

        // 7. recompute noise level (κ₀+κ₁)
        
        // Calculate new noise level
        // Calculate κ₀ = p * (k₂ + Σ(κᵢ + k₁ᵢ + κᵢk₁ᵢ))
        let kappa0 = {
            let mut c1_sum: u128 = 0;
            for &x in c1_noise.iter() {
                let (sum, overflow) = c1_sum.overflowing_add(x);
                if overflow {
                    println!("Warning: Overflow detected at line 171 in c1_noise sum: acc ({}) + x ({})", c1_sum, x);
                }
                c1_sum = sum;
            }

            let mut rho_sum: u128 = 0;
            for v in rho_noise.iter() {
                for &x in v.iter() {
                    let (sum, overflow) = rho_sum.overflowing_add(x);
                    if overflow {
                        println!("Warning: Overflow detected at line 172 in rho_noise sum");
                    }
                    rho_sum = sum;
                }
            }

            let mut c0_sum: u128 = 0;
            for v in c0_noise.iter() {
                for &x in v.iter() {
                    let (sum, overflow) = c0_sum.overflowing_add(x);
                    if overflow {
                        println!("Warning: Overflow detected at line 173 in c0_noise sum");
                    }
                    c0_sum = sum;
                }
            }

            let mut product_sum: u128 = 0;
            for (rf_noise, c0_noise) in rho_noise.iter().zip(c0_noise.iter()) {
                let mut rf_sum = 0u128;
                for &x in rf_noise.iter() {
                    let (sum, overflow) = rf_sum.overflowing_add(x);
                    if overflow {
                        println!("Warning: Overflow detected in rf_noise sum");
                    }
                    rf_sum = sum;
                }
                
                let mut c0_sum = 0u128;
                for &x in c0_noise.iter() {
                    let (sum, overflow) = c0_sum.overflowing_add(x);
                    if overflow {
                        println!("Warning: Overflow detected in c0_noise sum");
                    }
                    c0_sum = sum;
                }
                
                let (prod, overflow1) = rf_sum.overflowing_mul(c0_sum);
                if overflow1 {
                    println!("Warning: Overflow detected in product: rf_sum ({}) * c0_sum ({})", rf_sum, c0_sum);
                }
                let (sum, overflow2) = product_sum.overflowing_add(prod);
                if overflow2 {
                    println!("Warning: Overflow detected in sum accumulation: acc ({}) + prod ({})", product_sum, prod);
                }
                product_sum = sum;
            }

            // Final calculations with overflow detection
            let (sum1, overflow1) = c1_sum.overflowing_add(rho_sum);
            if overflow1 {
                println!("Warning: Overflow detected in final sum step 1");
            }
            let (sum2, overflow2) = sum1.overflowing_add(c0_sum);
            if overflow2 {
                println!("Warning: Overflow detected in final sum step 2");
            }
            let (sum3, overflow3) = sum2.overflowing_add(product_sum);
            if overflow3 {
                println!("Warning: Overflow detected in final sum step 3");
            }
            let (result, overflow4) = self.chan.p.overflowing_mul(sum3);
            if overflow4 {
                println!("Warning: Overflow detected in final multiplication: p ({}) * sum3 ({})",
                    self.chan.p, sum3);
            }
            result
        };

        let kappa1 = {
            let p = self.chan.p;
            let p_minus_1 = p.checked_sub(1)
                .expect("Arithmetic underflow in p-1");
            let p_minus_1_squared = p_minus_1.checked_mul(p_minus_1)
                .expect("Arithmetic overflow in (p-1)²");
            let dim_term = (self.chan.dim as u128).checked_mul(p_minus_1_squared)
                .expect("Arithmetic overflow in n(p-1)²");
            
            ((p_minus_1 + dim_term) / p) as u128
        };

        Cipher {
            dec: result.dec,
            enc: result.enc,
            level: kappa0.checked_add(kappa1)
                .expect("Arithmetic overflow in final level calculation"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn test_refresh_behavior() {
        // Use minimal parameters for clear testing
        let chan = ArithChannel::new(5, 17, 3, 2);
        let (aces, secret_key) = Aces::generate_keypair(&chan);
        let alg = AcesAlgebra::new(&chan, &secret_key);
        let refresher = Refresher::new(&aces, &alg, &chan);
        let mut rng = thread_rng();

        println!("\nTesting refresh with parameters:");
        println!("p={}, q={}, dim={}", chan.p, chan.q, chan.dim);
        println!("refresh threshold: {}", (chan.q + 1) / chan.p - 1);

        // Create and encrypt small message
        let m = 2u128;
        let (cipher, _) = aces.encrypt(m, &mut rng);
        
        // Record initial state
        let init_level = cipher.level;
        let init_eval = cipher.enc.eval(chan.omega) % chan.p;
        println!("\nInitial state:");
        println!("  Level: {}", init_level);
        println!("  [C](c) mod p: {}", init_eval);
        println!("  Dec(c): {}", aces.decrypt(&cipher, &secret_key));

        // Make refreshable if needed
        let cipher = refresher.make_refreshable(
            &cipher,
            &secret_key,
            &mut rng
        ).expect("Failed to make refreshable");

        // Perform refresh
        let refreshed = refresher.refresh(&cipher, &secret_key);
        
        // Check post-refresh state
        println!("\nPost-refresh state:");
        println!("  Level: {} -> {}", cipher.level, refreshed.level);
        println!("  [C](c) mod p: {} -> {}",
            cipher.enc.eval(chan.omega) % chan.p,
            refreshed.enc.eval(chan.omega) % chan.p);
        println!("  Dec(c): {} -> {}",
            aces.decrypt(&cipher, &secret_key),
            aces.decrypt(&refreshed, &secret_key));

        // Verify all properties
        assert_eq!(
            aces.decrypt(&refreshed, &secret_key), m,
            "Refresh changed the message"
        );
        assert!(
            refreshed.level <= cipher.level,
            "Level increased: {} -> {}",
            cipher.level, refreshed.level
        );
        assert_eq!(
            cipher.enc.eval(chan.omega) % chan.p,
            refreshed.enc.eval(chan.omega) % chan.p,
            "Evaluation mod p changed"
        );
        assert!(
            refresher.is_refreshable(&refreshed, &secret_key),
            "Result not refreshable"
        );
    }

    #[test]
    fn test_refresh_chain() {
        // Use minimal parameters
        let chan = ArithChannel::new(3, 11, 2, 2);
        let (aces, secret_key) = Aces::generate_keypair(&chan);
        let alg = AcesAlgebra::new(&chan, &secret_key);
        let refresher = Refresher::new(&aces, &alg, &chan);
        let mut rng = thread_rng();

        println!("\nRefresh Chain Test");
        println!("Parameters: p={}, q={}, dim={}", chan.p, chan.q, chan.dim);
        println!("Refresh threshold: {}", (chan.q + 1) / chan.p - 1);

        let m = 1u128;
        let (mut cipher, _) = aces.encrypt(m, &mut rng);
        let mut refresh_count = 0;

        // Try multiple refresh cycles
        for i in 0..3 {
            // Record pre-state
            let pre_level = cipher.level;
            let pre_eval = cipher.enc.eval(chan.omega) % chan.p;
            let pre_value = aces.decrypt(&cipher, &secret_key);

            println!("\nIteration {}:", i);
            println!("Pre-refresh:");
            println!("  Level: {}", pre_level);
            println!("  [C](c) mod p: {}", pre_eval);
            println!("  Dec(c): {}", pre_value);

            // Verify pre-state
            assert_eq!(pre_value, m, "Value changed before refresh");

            // Make refreshable if needed
            cipher = refresher.make_refreshable(
                &cipher,
                &secret_key,
                &mut rng
            ).expect("Failed to make refreshable");

            // Do refresh
            cipher = refresher.refresh(&cipher, &secret_key);
            refresh_count += 1;

            // Check post-state
            println!("Post-refresh:");
            println!("  Level: {}", cipher.level);
            println!("  [C](c) mod p: {}", cipher.enc.eval(chan.omega) % chan.p);
            println!("  Dec(c): {}", aces.decrypt(&cipher, &secret_key));

            // Verify refresh properties
            assert_eq!(
                aces.decrypt(&cipher, &secret_key),
                m,
                "Refresh changed value at iter {}", i
            );
            assert!(
                cipher.level <= pre_level,
                "Level increased at iter {}", i
            );
            assert_eq!(
                cipher.enc.eval(chan.omega) % chan.p,
                pre_eval,
                "Evaluation changed at iter {}", i
            );
            assert!(
                refresher.is_refreshable(&cipher, &secret_key),
                "Result not refreshable at iter {}", i
            );
        }

        println!("\nCompleted {} successful refreshes", refresh_count);
        assert!(
            refresh_count >= 2,
            "Should support at least 2 refreshes, got {}", refresh_count
        );
    }

    

    #[test]
    fn test_noise_control() {
        // Use parameters that allow noise growth
        let chan = ArithChannel::new(5, 31, 3, 2);
        let (aces, secret_key) = Aces::generate_keypair(&chan);
        let alg = AcesAlgebra::new(&chan, &secret_key);
        let refresher = Refresher::new(&aces, &alg, &chan);
        let mut rng = thread_rng();

        println!("\nTesting noise control");
        println!("Parameters: p={}, q={}, dim={}", chan.p, chan.q, chan.dim);
        println!("Refresh threshold: {}", (chan.q + 1) / chan.p - 1);

        // Start with a simple message
        let m = 2u128;
        let (mut cipher, _) = aces.encrypt(m, &mut rng);
        
        // Do multiple multiplications to increase noise
        for i in 0..3 {
            let (tmp, _) = aces.encrypt(1, &mut rng);
            cipher = alg.mult(&cipher, &tmp);
            
            println!("\nAfter {} multiplications:", i+1);
            println!("Level: {}", cipher.level);
            
            // Try refresh if possible
            if refresher.is_refreshable(&cipher, &secret_key) {
                let pre_level = cipher.level;
                cipher = refresher.refresh(&cipher, &secret_key);
                println!("Refreshed: {} -> {}", pre_level, cipher.level);
                
                // Verify refresh properties
                assert!(cipher.level < pre_level,
                    "Level did not decrease at step {}", i);
                assert_eq!(
                    aces.decrypt(&cipher, &secret_key),
                    m,
                    "Value changed at step {}", i
                );
            }
        }
    }


}
