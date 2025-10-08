//! Tools for projecting ACES instances onto the CRT component where the vanishing noise disappears.

use anyhow::{ensure, Context, Result};
use rand::Rng;
use serde::Serialize;

use crate::{arith_channel::ArithChannel, cipher::Cipher, polynomial::Polynomial, scheme::Aces};

pub mod modp;

use modp::{dot_mod_p, gaussian_elimination_mod_p};

/// Result of projecting a ciphertext via the public evaluation map.
#[derive(Debug, Clone, Serialize)]
pub struct ProjectedSample {
    /// Projected vector part Aᵢ = φ(cᵢ).
    pub a: Vec<u128>,
    /// Projected scalar part yᵢ = φ(cᵢ') − mᵢ (mod p).
    pub y: u128,
}

/// Output of the linear reconstruction attack.
#[derive(Debug, Clone, Serialize)]
pub struct AttackResult {
    /// Recovered public image of the secret vector s = φ(x).
    pub secret_image: Vec<u128>,
    /// Samples used to build the linear system (useful for reproducibility).
    pub samples: Vec<ProjectedSample>,
}

/// Apply the evaluation-and-reduction homomorphism φ to a polynomial.
#[must_use]
pub fn phi_polynomial(chan: &ArithChannel, poly: &Polynomial) -> u128 {
    chan.eval(poly) % chan.p
}

/// Project the vector component of a ciphertext via φ.
#[must_use]
pub fn project_dec_vector(chan: &ArithChannel, dec: &[Polynomial]) -> Vec<u128> {
    dec.iter().map(|poly| phi_polynomial(chan, poly)).collect()
}

/// Project an entire ciphertext, subtracting the chosen plaintext.
#[must_use]
pub fn project_cipher(chan: &ArithChannel, cipher: &Cipher, message: u128) -> ProjectedSample {
    let a = project_dec_vector(chan, &cipher.dec);
    let c_eval = phi_polynomial(chan, &cipher.enc);
    let m_mod = message % chan.p;
    let y = (c_eval + chan.p - m_mod) % chan.p;
    ProjectedSample { a, y }
}

/// Recover φ(x) from chosen ciphertext samples of encryption.
pub fn run_attack<R: Rng>(
    aces: &Aces,
    rng: &mut R,
    sample_count: usize,
    message: u128,
) -> Result<AttackResult> {
    let chan = aces.channel();
    ensure!(
        sample_count >= chan.dim,
        "need at least {} samples to determine φ(x)",
        chan.dim
    );

    let msg_mod = message % chan.p;
    let mut samples = Vec::with_capacity(sample_count);
    for _ in 0..sample_count {
        let ct = aces.encrypt(msg_mod, rng);
        samples.push(project_cipher(chan, &ct, msg_mod));
    }

    let mut matrix: Vec<Vec<u128>> = samples.iter().map(|s| s.a.clone()).collect();
    let mut rhs: Vec<u128> = samples.iter().map(|s| s.y).collect();
    let secret_image = gaussian_elimination_mod_p(&mut matrix, &mut rhs, chan.p, chan.dim)
        .context("linear system is rank deficient; try more samples")?;

    Ok(AttackResult {
        secret_image,
        samples,
    })
}

/// Decrypt a ciphertext using only the recovered φ(x).
#[must_use]
pub fn decrypt_with_projected_secret(
    chan: &ArithChannel,
    cipher: &Cipher,
    secret_image: &[u128],
) -> u128 {
    let a = project_dec_vector(chan, &cipher.dec);
    let c_eval = phi_polynomial(chan, &cipher.enc);
    let dot = dot_mod_p(&a, secret_image, chan.p);
    (c_eval + chan.p - dot) % chan.p
}
