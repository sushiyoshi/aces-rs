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
    
    pub fn fractional_part_f64(a: u128, b: u128) -> f64 {
        // (1) 整数部を取り除いた余り
        let mut rem = a % b;
    
        // (2) 結果と、現在のビット重み
        let mut result = 0.0f64;
        let mut factor = 0.5f64;
    
        // 2*rem >= b の判定を、溢れを起こさずに行うための閾値
        // b + 1 がオーバーフローする場合も検知してパニック
        let threshold = b
            .checked_add(1)
            .expect("overflow computing threshold for (b + 1)")
            / 2;
    
        // 二進法の小数ビットを最大 64 桁抽出
        for _ in 0..53 {
            if rem == 0 {
                break;
            }
            if rem >= threshold {
                // bit = 1
                // rem = 2*rem - b を checked_* で
                rem = rem
                    .checked_shl(1)
                    .expect("overflow on rem << 1");
                rem = rem
                    .checked_sub(b)
                    .expect("overflow on (2*rem - b)");
                result += factor;
            } else {
                // bit = 0, rem = 2*rem
                rem = rem
                    .checked_shl(1)
                    .expect("overflow on rem << 1");
            }
            factor *= 0.5;
        }
    
        result
    }
    
    pub fn is_refreshable(&self, ct: &Cipher) -> bool {
        let p     = self.chan.p  as u128;
        let q     = self.chan.q  as u128;
        let omega = self.chan.omega as u128;
        if ct.level > q/p {
            println!("Ciphertext level {} exceeds threshold {}", ct.level, q/p);
            return false;
        }
        // (1) ℓ := 〚C〛((c))
        let ell: Vec<u128> = ct.dec
            .iter()
            .map(|poly| (poly % &self.chan.u).eval(omega) % q)
            .collect();

           
        // (2) ω_q ∘ 〚C〛((x)) 
        let x_eval = &self.scheme.x_eval;

        // (3) Definition 5.41: p-locator 判定
        let dot: u128 = ell.iter()
        .zip(x_eval.iter())
        .fold(0u128, |acc, (e, x)| {
            let (prod, overflow1) = e.overflowing_mul(*x);
            if overflow1 {
                panic!("Overflow occurred during multiplication: {} * {}", e, x);
            }
            let (sum, overflow2) = acc.overflowing_add(prod);
            if overflow2 {
                panic!("Overflow occurred during addition: {} + {}", acc, prod);
            }
            sum
        });
    
        let x_eval_sum = self.scheme.x_eval_sum;
        let dot_div = dot / q;
        if x_eval_sum <= dot_div {
            println!("x_eval_sum {} < dot_div {}", x_eval_sum, dot_div);
            return false;
        }
        let diff = x_eval_sum - dot_div;
        if diff % p != 0 {
            return false;
        }
        let margin = Refresher::fractional_part_f64(dot, q);
        // let margin = (dot as f64 / q as f64) - (dot_div as f64);
        if margin < 0.0 {
            panic!("Margin is negative: {}", margin);
        } else {
            println!("Margin: {}", margin);
        }
        println!("diff/p: {}", diff / p);
        
        // (4) Definition 5.43: margin ∈ [0,1)
        let margin = Refresher::fractional_part_f64(dot, q);
        if margin < 0.0 {
            panic!("Margin is negative: {}", margin);
        } else {
            println!("Margin: {}", margin);
        }
        let lhs: f64 = (p * (ct.level + 1) - 1) as f64 / q as f64;
        println!("lhs: {}, rhs: {}", lhs, 1.0 - margin);
        lhs < 1.0 - margin
    }

    /// Make a ciphertext refreshable by repeatedly adding encrypted zeros
    /// Will panic after MAX_ATTEMPTS (default: 10) unsuccessful attempts
    pub fn make_refreshable<R: Rng>(
        &self,
        cipher: &Cipher,
        rng: &mut R,
    ) -> (Option<Cipher>, u32) {
        const MAX_ATTEMPTS: u32 = 10000;

        // First check if already refreshable
        if self.is_refreshable(cipher) {
            return (Some(cipher.clone()), 0);
        }

        // Check initial level bound
        let level_bound = (self.chan.q + 1) / self.chan.p - 1;
        let mut current = cipher.clone();
        let mut attempts = 0;

        // Keep trying new zero ciphers until we get a refreshable result
        loop {
            attempts += 1;
            if attempts > MAX_ATTEMPTS {
                panic!("Failed to make cipher refreshable after {} attempts", MAX_ATTEMPTS);
            }

            // Generate a fresh encryption of zero
            let zero_cipher = self.scheme.encrypt(0, rng);
            
            // Add it to our current cipher
            let next = self.alg.add(&current, &zero_cipher);
            // println!("nextlevel: {}, level_bound: {}", next.level, level_bound);
            // Check if level would exceed bound
            if next.level >= level_bound {
                println!("Level would exceed bound: {} >= {}",next.level, level_bound);
                panic!("Level exceeded bound after {} attempts", attempts);
            }
            
            // current = next;
            if self.is_refreshable(&next) {
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
        let rho_cipher = secret_x
            .iter()
            .map(|xi| {
                let xi_u = xi % &self.chan.u;
                let eval = xi_u.eval(self.chan.omega) % self.chan.q % self.chan.p;
                self.scheme.encrypt(eval, &mut rng)
            })
            .collect::<Vec<_>>();
            

        // 2. Create negated c₀ components
        let c0_int: Vec<u128> = c.dec
            .iter()
            .map(|ci| {
                let ci_u = ci % &self.chan.u;
                let eval = ci_u.eval(self.chan.omega) % self.chan.q;
                (self.chan.q - eval) % self.chan.q
            })
            .collect();
        
        let c0_enc = c0_int
            .iter()
            .map(|&m| self.scheme.encrypt(m % self.chan.p, &mut rng))
            .collect::<Vec<_>>();
            
            
        // 3. Get evaluation of c₁ = [C](c)
        let c1_eval = c.enc.eval(self.chan.omega) % self.chan.q;
        let c1_mod_p = c1_eval % self.chan.p;
        let c1_enc = self.scheme.encrypt(c1_mod_p, &mut rng);

        // 5. Scalar product ⟨c0_enc, ρ⟩
        let mut sp = self.alg.mult(&c0_enc[0], &rho_cipher[0]);
        for i in 1..rho_cipher.len() {
            sp = self.alg.add(&sp, &self.alg.mult(&c0_enc[i], &rho_cipher[i]));
        }
        
        // 6. Add c1_enc
        let result = self.alg.add(&c1_enc, &sp);

        // 7. recompute noise level (κ₀+κ₁)
        
        // Calculate new noise level
        let kappa_asterisk_star = {
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
            level: result.level.checked_add(kappa_asterisk_star)
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
        let cipher = aces.encrypt(m, &mut rng);
        
        // Record initial state
        let init_level = cipher.level;
        let init_eval = cipher.enc.eval(chan.omega) % chan.p;
        println!("\nInitial state:");
        println!("  Level: {}", init_level);
        println!("  [C](c) mod p: {}", init_eval);
        println!("  Dec(c): {}", aces.decrypt(&cipher, &secret_key));

        // Make refreshable if needed
        let cip = refresher.make_refreshable(
            &cipher,
            &mut rng
        );
        let cipher = cip.0.expect("Failed to make refreshable");

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
            refresher.is_refreshable(&refreshed),
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
        let mut cipher = aces.encrypt(m, &mut rng);
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
            let cip = refresher.make_refreshable(
                &cipher,
                &mut rng
            );
            
            cipher = cip.0.expect("Failed to make refreshable");

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
                refresher.is_refreshable(&cipher),
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
        let mut cipher = aces.encrypt(m, &mut rng);
        
        // Do multiple multiplications to increase noise
        for i in 0..3 {
            let tmp = aces.encrypt(1, &mut rng);
            cipher = alg.mult(&cipher, &tmp);
            
            println!("\nAfter {} multiplications:", i+1);
            println!("Level: {}", cipher.level);
            
            // Try refresh if possible
            if refresher.is_refreshable(&cipher) {
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
