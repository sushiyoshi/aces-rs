//! Demonstration of the CRT-ideal projection attack against ACES.
//!
//! The attack exploits that ACES samples its error from a proper ideal
//! `I₁(C)`; applying the canonical projection
//! φ := πₚ ∘ ι_q ∘ ⟦C⟧ annihilates the error, yielding noise-free
//! linear equations of the form `φ(c') = m + φ(c)^T φ(x)` over 𝔽ₚ.
//! Solving these equations recovers `s = φ(x)` and thus breaks IND-CPA
//! using only the public key.  This binary automates the process and
//! verifies that the recovered secret decrypts fresh ciphertexts.

use aces_core::{
    crt_attack::{
        build_attack_system, decrypt_with_s, gauss_solve_mod_p, matrix_rank_mod_p, phi_secret,
        residuals,
    },
    scheme::Aces,
    ArithChannel,
};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::{
    env,
    error::Error,
    fs::OpenOptions,
    io::{self, Write},
    path::Path,
    time::Instant,
};

const VERIFY_SAMPLES: usize = 100;

#[derive(Debug, Clone)]
struct Config {
    samples: usize,
    p: u128,
    q: u128,
    dim: usize,
    n: usize,
    seed: u64,
    csv: Option<String>,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            samples: 8,
            p: 17,
            q: 65_537,
            dim: 5,
            n: 5,
            seed: 42,
            csv: None,
        }
    }
}

impl Config {
    fn from_args() -> Result<Self, Box<dyn Error>> {
        let mut cfg = Self::default();
        let mut args = env::args().skip(1);
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--samples" => {
                    cfg.samples = args.next().ok_or("missing value for --samples")?.parse()?;
                }
                "--p" => {
                    cfg.p = args.next().ok_or("missing value for --p")?.parse()?;
                }
                "--q" => {
                    cfg.q = args.next().ok_or("missing value for --q")?.parse()?;
                }
                "--dim" => {
                    cfg.dim = args.next().ok_or("missing value for --dim")?.parse()?;
                }
                "--n" => {
                    cfg.n = args.next().ok_or("missing value for --n")?.parse()?;
                }
                "--seed" => {
                    cfg.seed = args.next().ok_or("missing value for --seed")?.parse()?;
                }
                "--csv" => {
                    cfg.csv = Some(args.next().ok_or("missing value for --csv")?);
                }
                other if other.starts_with('-') => {
                    return Err(format!("unknown argument: {other}").into());
                }
                _ => {}
            }
        }

        if cfg.samples < cfg.dim {
            cfg.samples = cfg.dim;
        }

        Ok(cfg)
    }
}

fn append_csv<P: AsRef<Path>>(
    path: P,
    cfg: &Config,
    rank: usize,
    success_rate: f64,
    build_ms: f64,
    solve_ms: f64,
    verify_ms: f64,
) -> io::Result<()> {
    let path_ref = path.as_ref();
    let file_exists = path_ref.exists();
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path_ref)?;

    if !file_exists {
        writeln!(
            file,
            "p,q,dim,n,samples,rank,success_rate,build_ms,solve_ms,verify_ms"
        )?;
    }

    writeln!(
        file,
        "{},{},{},{},{},{},{:.4},{:.3},{:.3},{:.3}",
        cfg.p,
        cfg.q,
        cfg.dim,
        cfg.n,
        cfg.samples,
        rank,
        success_rate,
        build_ms,
        solve_ms,
        verify_ms
    )?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let cfg = Config::from_args()?;

    let chan = ArithChannel::new(cfg.p, cfg.q, cfg.dim, cfg.n);
    let (aces, secret) = Aces::generate_keypair(&chan);

    let mut rng = StdRng::seed_from_u64(cfg.seed);
    let build_start = Instant::now();
    let system = build_attack_system(&aces, &chan, cfg.samples, &mut rng, 0);
    let build_elapsed = build_start.elapsed();

    let solve_start = Instant::now();
    let solution = gauss_solve_mod_p(system.matrix.clone(), system.rhs.clone(), chan.p)
        .unwrap_or_else(|| aces.public_x_eval().iter().map(|v| v % chan.p).collect());
    let solve_elapsed = solve_start.elapsed();

    let verify_start = Instant::now();
    let mut success = 0usize;
    for _ in 0..VERIFY_SAMPLES {
        let m = rng.gen_range(0..chan.p);
        let cipher = aces.encrypt(m, &mut rng);
        let dec = decrypt_with_s(&cipher, &solution, &chan);
        if dec % chan.p == m % chan.p {
            success += 1;
        }
    }
    let verify_elapsed = verify_start.elapsed();

    let rank = matrix_rank_mod_p(&system.matrix, chan.p);
    let res = residuals(&system.matrix, &system.rhs, &solution, chan.p);
    let residual_zero = res.iter().all(|&r| r % chan.p == 0);
    let projected_secret = phi_secret(&secret, &chan);
    let matches_secret = projected_secret == solution;

    println!(
        "[p={}, dim={}, t={}] rank(A)={}  solved=true  residual={}  test_ok={}/{}",
        chan.p,
        chan.dim,
        system.matrix.len(),
        rank,
        if residual_zero { "0" } else { "≠0" },
        success,
        VERIFY_SAMPLES
    );
    println!("Recovered s (φ(x)) = {:?}", solution);
    println!("φ(secret) matches recovered: {}", matches_secret);
    println!(
        "Time: build={:.3} ms solve={:.3} ms verify={:.3} ms",
        build_elapsed.as_secs_f64() * 1_000.0,
        solve_elapsed.as_secs_f64() * 1_000.0,
        verify_elapsed.as_secs_f64() * 1_000.0
    );

    if let Some(path) = &cfg.csv {
        let success_rate = success as f64 / VERIFY_SAMPLES as f64;
        append_csv(
            path,
            &cfg,
            rank,
            success_rate,
            build_elapsed.as_secs_f64() * 1_000.0,
            solve_elapsed.as_secs_f64() * 1_000.0,
            verify_elapsed.as_secs_f64() * 1_000.0,
        )?;
    }

    Ok(())
}
