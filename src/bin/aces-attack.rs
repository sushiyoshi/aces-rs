//! Proof-of-concept attack exploiting CRT projection to recover φ(x).

use std::fmt::Write as _;

use anyhow::{Context, Result};
use clap::Parser;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde_json::json;

use aces_core::{
    attack::{decrypt_with_projected_secret, phi_polynomial, run_attack},
    scheme::Aces,
    ArithChannel,
};

#[derive(Parser, Debug)]
#[command(author, version, about = "Recover φ(x) via CRT projection and decrypt", long_about = None)]
struct Args {
    /// Number of plaintext/ciphertext samples to collect.
    #[arg(long)]
    samples: Option<usize>,
    /// Override the prime p of the arithmetic channel.
    #[arg(long, value_name = "P")]
    prime_p: Option<u128>,
    /// Seed for deterministic randomness (default 0).
    #[arg(long, default_value_t = 0)]
    seed: u64,
    /// Emit JSON transcript of the attack.
    #[arg(long, default_value_t = false)]
    json: bool,
}

const DEFAULT_P: u128 = 32;
const DEFAULT_Q: u128 = (1u128 << 40) + 1;
const DEFAULT_DIM: usize = 5;
const DEFAULT_N: usize = 5;
const DEFAULT_MESSAGE: u128 = 0;
const VERIFY_SAMPLES: usize = 10;

fn main() -> Result<()> {
    let args = Args::parse();
    let p = args.prime_p.unwrap_or(DEFAULT_P);
    let chan = ArithChannel::new(p, DEFAULT_Q, DEFAULT_DIM, DEFAULT_N);
    let (aces, secret_key) = Aces::generate_keypair(&chan);

    let mut rng = StdRng::seed_from_u64(args.seed);
    let min_samples = aces.channel().dim + 12;
    let sample_target = args
        .samples
        .map(|s| s.max(aces.channel().dim))
        .unwrap_or(min_samples);

    let attack = run_attack(&aces, &mut rng, sample_target, DEFAULT_MESSAGE)
        .context("attack failed: insufficient rank over F_p")?;

    let actual_phi: Vec<u128> = secret_key
        .iter()
        .map(|poly| phi_polynomial(aces.channel(), poly))
        .collect();

    let mut success = 0usize;
    let mut failures = Vec::new();
    for _ in 0..VERIFY_SAMPLES {
        let message = rng.gen::<u128>() % aces.channel().p;
        let cipher = aces.encrypt(message, &mut rng);
        let recovered =
            decrypt_with_projected_secret(aces.channel(), &cipher, &attack.secret_image);
        if recovered == message % aces.channel().p {
            success += 1;
        } else {
            failures.push((message, recovered));
        }
    }

    if args.json {
        let report = json!({
            "parameters": {
                "p": chan.p,
                "q": chan.q,
                "dim": chan.dim,
                "n": chan.n,
            },
            "collected_samples": attack.samples.len(),
            "recovered_phi_x": attack.secret_image,
            "actual_phi_x": actual_phi,
            "verification": {
                "trials": VERIFY_SAMPLES,
                "successes": success,
                "failures": failures,
            },
            "samples": attack.samples,
        });
        println!("{}", serde_json::to_string_pretty(&report)?);
    } else {
        let mut summary = String::new();
        writeln!(
            &mut summary,
            "Recovered φ(x) (mod p): {:?}",
            attack.secret_image
        )?;
        writeln!(&mut summary, "Actual φ(x) (mod p):    {:?}", actual_phi)?;
        writeln!(
            &mut summary,
            "Verification: {}/{} ciphertexts decrypted correctly",
            success, VERIFY_SAMPLES
        )?;
        if !failures.is_empty() {
            writeln!(&mut summary, "Failures: {:?}", failures)?;
        }
        print!("{}", summary);
    }

    Ok(())
}
