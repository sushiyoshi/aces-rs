# aces-rs

Rust implementation of **ACES (Arithmetic Channel Encryption Scheme)**, based on the original Python prototype by \[[arXiv:2401.13255](https://arxiv.org/abs/2401.13255)].  This repository targets **research‑grade benchmarking** of ACES on modern x86‑64 CPUs (macOS M2 under Rosetta 2, Intel/AMD Linux, etc.).

---

## Features

* **Fully typed, idiomatic Rust** (`edition = 2021`)
* **AVX / SSE optimized** polynomial arithmetic via `RUSTFLAGS=-C target-feature=+avx`
* Refreshability check and noise tracking helpers
* CSV benchmark logger for large test campaigns
* Easy parameter tweaking (`p`, `q`, `dim`, `n`)
---

## Building

```bash
# clone your own fork
git clone https://github.com/sushiyoshi/aces-rs.git
cd aces-rs

# debug build
cargo run

# release build with AVX acceleration
RUSTFLAGS='-C target-feature=+avx' cargo run --release
```

---

## Benchmark results

All timings collected on a **MacBook Pro 14‑inch (M2 Pro, 12‑core CPU, 32 GB RAM)**
(Host macOS 14.x; Rust 1.78; default `--release` unless noted).

### p = 16 (4‑bit modulus)

| Build                           | Mul count | Wall‑clock   | Avg refresh  | Max refresh | Min refresh |
| ------------------------------- | --------- | ------------ | ------------ | ----------- | ----------- |
| `cargo run` (debug)             | 1 000 000 | 2940.93 s    | 4.879 ms     | 94.673 ms   | 4.498 ms    |
| `cargo run` (debug, second run) | 1 000 000 | 3188.88 s    | 5.399 ms     | 145.201 ms  | 4.911 ms    |
| `+avx --release`                | 1 000 000 | **704.57 s** | **1.279 ms** | 54.653 ms   | 1.169 ms    |

### p = 32 (5‑bit modulus)

| Build       | Mul count | Wall‑clock | Avg refresh | Max refresh | Min refresh |
| ----------- | --------- | ---------- | ----------- | ----------- | ----------- |
| `cargo run` | 100 000   | 898.97 s   | 3.946 ms    | 34.328 ms   | 3.672 ms    |

---

## Quick start

```rust
use aces_core::{Aces, AcesAlgebra, ArithChannel, Polynomial, Refresher};
use rand::{Rng, thread_rng};
use serde::de;
use std::time::{Duration, Instant};
use std::error::Error;
fn main() -> Result<(), Box<dyn Error>> {

    let max_mults = 10000;
    // Initialize optimal parameters from Python implementation
    let p: u128 = 16;  // small prime
    let q: u128 = 16_u128.pow(15) + 1; // small prime
    let dim = 5;       // optimal polynomial dimension
    let n = 5;         // optimal number of components

    println!("Initializing ACES with optimal parameters:");
    println!("p = {}, q = {}, dim = {}, n = {}", p, q, dim, n);

    // Create arithmetic channel and ACES instance
    let chan = ArithChannel::new(p, q, dim, n);
    // Generate keypair for encryption and decryption
    let (aces, secret_key) = Aces::generate_keypair(&chan);
    // Create algebra with proper tensor for multiplication
    let alg = AcesAlgebra::new(&chan, &secret_key);
    
    // Prepare for multiplication test
    let mut rng = thread_rng();
    
    let mut refresh_times = Vec::new();
    let mut successful_ops = 0;
    let mut refresh_success = 0;
    
    // println!("\nStarting multiplication chain test ({} operations):", max_mults);
    let start = Instant::now();

    // Initial message
    let m: u128 = rng.gen_range(1..p);
    let mut cur_value = m;
    let (mut cipher, _) = aces.encrypt(m, &mut rng);
    
    println!("Initial message: {}", m);
    
    // Use generated secret key for decryption
    let decrypted = aces.decrypt(&cipher, &secret_key);
    println!("Decrypted initial message: {}", decrypted);
    assert_eq!(decrypted, m, "Initial decryption failed");
    let noise_max = (q + 1)/p -1;

    // Test loop
    for i in 0..max_mults {
        // Generate random value in [1, p-1]
        let mi = rng.gen::<u128>() % (p - 1) + 1;
        // println!("mi: {}", mi);
        let (cipher_i, _) = aces.encrypt(mi, &mut rng);
        
        // Perform multiplication (cur_value == 0: add, else: multiply)
        cipher = if cur_value == 0 {
            alg.add(&cipher, &cipher_i)
        } else {
            alg.mult(&cipher, &cipher_i)
        };
        cur_value = if cur_value == 0 {
            (cur_value + mi) % p
        } else {
            (cur_value * mi) % p
        };
        let initial_level = (n as u128) * (p as u128);
        
        let next_level = (cipher.level + initial_level + cipher.level * initial_level) * p;
        let refresh_frag = next_level > noise_max;

        let mut zero_adds = 0;
        // Check if refresh is needed
        let refresher = Refresher::new(&aces, &alg, &chan);
        if refresh_frag {
            if !refresher.is_refreshable(&cipher, &secret_key) {
                // println!("Making cipher refreshable at operation {}", i + 1);
                let cipher_tuple = refresher.make_refreshable(
                    &cipher,
                    &secret_key,
                    &mut rng
                );
                cipher = cipher_tuple.0.expect("Failed to make refreshable");
                zero_adds = cipher_tuple.1;
                // let decrypted = aces.decrypt(&cipher, &secret_key);
                // println!("Decrypted after making refreshable: {}", decrypted);
            }
            
            // println!("Performing refresh at operation {}", i + 1);
            let refresh_start = Instant::now();
            cipher = refresher.refresh(&cipher, &secret_key);
            refresh_times.push(refresh_start.elapsed());

            let elapsed_ms = refresh_start.elapsed().as_secs_f64() * 1000.0;
            println!("Refresh time: {:.6} ms", elapsed_ms);
            println!("zero_adds: {}", zero_adds);
           
            }

        // Verify result
        let decrypted = aces.decrypt(&cipher, &secret_key);
        if decrypted == cur_value {
            successful_ops += 1;
            // println!("Operation {} successful: expected {}, got {}",i + 1, cur_value, decrypted);
            if refresh_frag {
                refresh_success += 1;
                // println!("Refresh successful at operation {}", i + 1);
            }
        } else {
            println!("Operation {} failed: expected {}, got {}",i + 1, cur_value, decrypted);
            break;
        }
    }

    let total_time = start.elapsed();
    
    // Print results
    println!("\nTest completed:");
    println!("Total time: {:?}", total_time);
    println!("Successful operations: {}/{}", successful_ops, max_mults);
    println!("Refresh operations: {}/{}", refresh_success, successful_ops);
    
    if !refresh_times.is_empty() {
        let avg_refresh = refresh_times.iter().sum::<Duration>() / refresh_times.len() as u32;
        let max_refresh = refresh_times.iter().max().unwrap();
        let min_refresh = refresh_times.iter().min().unwrap();
        println!("\nRefresh statistics:");
        println!("Count: {}", refresh_times.len());
        println!("Average time: {:?}", avg_refresh);
        println!("Maximum time: {:?}", max_refresh);
        println!("Minimum time: {:?}", min_refresh);
    }
    Ok(())
}


```

Compile and run with

```bash
RUSTFLAGS='-C target-feature=+avx' cargo run --release --example quickstart
```

---
## Note on Refreshability Checker

> **Warning**: The current implementation of the *refreshability check* is **incomplete**.
> At present, `is_refreshable` uses the **secret key** directly to determine whether a ciphertext is refreshable.
> This shortcut was adopted for early benchmarking purposes and **does not represent a fully secure or key-independent check**.
>
> In future updates, a proper key-independent heuristic or estimator will be implemented to ensure refreshability can be evaluated **without access to the secret key**.

Developers using `Refresher::is_refreshable` should be aware that it is only a temporary placeholder.
Use with caution in security-sensitive scenarios.

---

## License

`aces-rs` is distributed under the **MIT License**.

Portions of the code are derivative works of the reference Python prototype by  
Rémy Tuyéras (<https://github.com/remytuyeras/aces>) and are reused here under the
same MIT terms.  Both copyright notices are included in `LICENSE`, and
third-party attributions are summarized in `NOTICE`.

