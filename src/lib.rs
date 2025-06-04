//! ACES ― Arithmetic Channel Encryption Scheme  (research prototype)

#![forbid(unsafe_code)]
#![warn(clippy::pedantic, missing_docs)]

extern crate rand;

pub mod polynomial;
pub mod cipher;
pub mod arith_channel;
pub mod scheme;
pub mod algebra;
pub mod refresher;
pub mod rand_iso; 
pub mod ntt;          

pub use arith_channel::ArithChannel;
pub use cipher::Cipher;
pub use polynomial::Polynomial;
pub use scheme::Aces;
pub use algebra::AcesAlgebra;
pub use refresher::Refresher;
