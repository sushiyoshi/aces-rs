//! ACES ciphertext container.

use crate::polynomial::Polynomial;

/// (c₀, c₁, level)
#[derive(Clone, Debug)]
pub struct Cipher {
    pub dec: Vec<Polynomial>, // vector part: length = dim
    pub enc: Polynomial,      // scalar part
    pub level: u128,          // noise level tracker (fits in u128 for demo)
}
