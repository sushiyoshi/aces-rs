//! ACES demo implementation with optimal parameters

use aces_core::{Aces, AcesAlgebra, ArithChannel, Polynomial, Refresher};
use rand::{Rng, thread_rng};
use serde::de;
use std::time::{Duration, Instant};
use csv::Writer;         
use std::error::Error;
fn main() -> Result<(), Box<dyn Error>> {
    // -------------- 計測用 CSV ライタ --------------
    let mut wtr = Writer::from_path("mul_stats.csv")?;
    wtr.write_record(&["index", "time_ms", "zero_adds"])?;

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
    
    // println!("Initial message: {}", m);
    
    // Use generated secret key for decryption
    let decrypted = aces.decrypt(&cipher, &secret_key);
    // println!("Decrypted initial message: {}", decrypted);
    // assert_eq!(decrypted, m, "Initial decryption failed");
    let noise_max = (q + 1)/p -1;

    // Test loop
    for i in 0..max_mults {
        // Generate random value in [1, p-1]
        let mi = rng.gen::<u128>() % (p - 1) + 1;
        // println!("mi: {}", mi);
        let (cipher_i, _) = aces.encrypt(mi, &mut rng);
        
        // // Perform multiplication (cur_value == 0: add, else: multiply)
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

        // cipher = alg.add(&cipher, &cipher_i);
        // cur_value = (cur_value + mi) % p;

        // decrypt the result
        // let decrypted = aces.decrypt(&cipher, &secret_key);
        // println!("Decrypted operations: {}", decrypted);
        // println!("Current value: {}", cur_value);
        // assert_eq!(decrypted, cur_value, "Decryption failed at operation {}", i + 1);
        
        // Check if noise is within bounds
        // calculate next noise level
        // let next_uplvl = (ciph.uplvl + initial_uplvl + ciph.uplvl * initial_uplvl) * p
        let initial_level = (n as u128) * (p as u128);
        
        let next_level = (cipher.level + initial_level + cipher.level * initial_level) * p;
        let refresh_frag = next_level > noise_max;

        let mut zero_adds = 0;
        // Check if refresh is needed
        let refresher = Refresher::new(&aces, &alg, &chan);
        if refresh_frag {
            // println!("Refresh start");
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
            wtr.write_record(&[
                i.to_string(),
                format!("{:.6}", elapsed_ms),
                zero_adds.to_string(),
            ])?;
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
    wtr.flush()?;           // 忘れずフラッシュ
    Ok(())
}
