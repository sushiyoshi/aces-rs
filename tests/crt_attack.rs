use aces_core::{
    crt_attack::{
        build_attack_system, decrypt_with_s, matrix_rank_mod_p, phi_poly, phi_secret, phi_vec,
        recover_s, residuals,
    },
    scheme::Aces,
    ArithChannel,
};
use rand::{rngs::StdRng, Rng, SeedableRng};

fn default_channel() -> ArithChannel {
    let p = 17;
    let q = 65_537;
    let dim = 5;
    let n = 5;
    ArithChannel::new(p, q, dim, n)
}

#[test]
fn recovers_projected_secret() {
    let chan = default_channel();
    let (aces, secret) = Aces::generate_keypair(&chan);
    let mut rng = StdRng::seed_from_u64(42);
    let f0_proj: Vec<Vec<u128>> = aces.f0.iter().map(|row| phi_vec(row, &chan)).collect();
    assert!(f0_proj
        .iter()
        .all(|row| row.iter().all(|&v| v % chan.p == 0)));
    let f1_proj: Vec<u128> = phi_vec(&aces.f1, &chan);
    assert!(f1_proj.iter().all(|&v| v % chan.p == 0));
    let secret_proj = phi_secret(&secret, &chan);

    let probe = aces.encrypt(1, &mut rng);
    assert_eq!(phi_poly(&probe.enc, &chan), 1 % chan.p);

    let system = build_attack_system(&aces, &chan, chan.dim + 3, &mut rng, 0);
    let solution = recover_s(&aces, &chan, chan.dim + 3, &mut rng);
    assert_eq!(
        solution, secret_proj,
        "projected secret should be recovered"
    );

    let residual = residuals(&system.matrix, &system.rhs, &solution, chan.p);
    assert!(residual.iter().all(|&r| r == 0), "residual should vanish");

    let rank = matrix_rank_mod_p(&system.matrix, chan.p);
    assert!(rank <= chan.dim);

    let mut fresh_rng = StdRng::seed_from_u64(42);
    let via_helper = recover_s(&aces, &chan, chan.dim + 3, &mut fresh_rng);
    assert_eq!(via_helper, solution);
}

#[test]
fn decrypts_random_ciphertexts_without_secret() {
    let chan = default_channel();
    let (aces, secret) = Aces::generate_keypair(&chan);
    let mut rng = StdRng::seed_from_u64(1337);

    let _system = build_attack_system(&aces, &chan, chan.dim + 4, &mut rng, 0);
    let solution = recover_s(&aces, &chan, chan.dim + 4, &mut rng);

    let phi_secret = phi_secret(&secret, &chan);
    assert_eq!(solution, phi_secret);

    for _ in 0..50 {
        let m = rng.gen_range(0..chan.p);
        let cipher = aces.encrypt(m, &mut rng);
        let dec = decrypt_with_s(&cipher, &solution, &chan);
        assert_eq!(dec % chan.p, m % chan.p, "decryption should match message");
    }
}
