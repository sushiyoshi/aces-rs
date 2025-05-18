# aces-rs

Rust implementation of **ACES (Arithmetic Channel Encryption Scheme)**, based on the original Python prototype by Tuyer & Rotaru \[[arXiv:2401.13255](https://arxiv.org/abs/2401.13255)].  This repository targets **research‑grade benchmarking** of ACES on modern x86‑64 CPUs (macOS M2 under Rosetta 2, Intel/AMD Linux, etc.).

---

## Features

* **Fully typed, idiomatic Rust** (`edition = 2021`)
* **AVX / SSE optimized** polynomial arithmetic via `RUSTFLAGS=-C target-feature=+avx`
* Refreshability check and noise tracking helpers
* CSV benchmark logger for large test campaigns
* Easy parameter tweaking (`p`, `q`, `dim`, `n`)

---

## Prerequisites

| Tool           | Minimum version | Notes                            |
| -------------- | --------------- | -------------------------------- |
| Rust toolchain | 1.77 (stable)   | Install with `rustup`            |
| CMake          | 3.25            | Needed by some transitive C deps |
| NASM / yasm    | latest          | Optional, for hand‑rolled ASM    |

> **macOS on Apple Silicon**
> Runs fine under Rosetta 2. Native aarch64 support is planned (ARM NEON intrinsics).

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
use aces_core::{Aces, ArithChannel};

fn main() {
    // Parameters
    let p  = 16u128;
    let q  = 1_152_921_504_606_846_977u128; // 2^60 + 17
    let dim = 5;
    let n   = 5;

    // Initialize
    let mut rng  = rand::thread_rng();
    let mut aces = Aces::new(p, q, dim, n, &mut rng);

    // Encrypt two plaintexts
    let ct1 = aces.encrypt(3u128);
    let ct2 = aces.encrypt(29u128);

    // Multiply, refresh if needed, decrypt
    let mut prod = aces.mul(&ct1, &ct2);
    if aces.is_refreshable(&prod) {
        prod = aces.refresh(prod);
    }
    let m = aces.decrypt(&prod);
    println!("3 × 29 mod {p} = {m}");
}
```

Compile and run with

```bash
RUSTFLAGS='-C target-feature=+avx' cargo run --release --example quickstart
```

---

## License

`aces-rs` is distributed under the **MIT License**.

Portions of the code are derivative works of the reference Python prototype by  
Rémy Tuyéras (<https://github.com/remytuyeras/aces>) and are reused here under the
same MIT terms.  Both copyright notices are included in `LICENSE`, and
third-party attributions are summarized in `NOTICE`.
