use aces_core::{
    attack::{decrypt_with_projected_secret, phi_polynomial, run_attack},
    scheme::Aces,
    ArithChannel,
};
use rand::{rngs::StdRng, Rng, SeedableRng};

#[test]
fn recover_secret_image_via_projection() {
    let p: u128 = 32;
    let q: u128 = (1u128 << 40) + 1;
    let dim = 5usize;
    let n = 5usize;

    let chan = ArithChannel::new(p, q, dim, n);
    let (aces, secret_key) = Aces::generate_keypair(&chan);

    let mut rng = StdRng::seed_from_u64(42);
    let samples = dim + 12;
    let attack = run_attack(&aces, &mut rng, samples, 0).expect("attack should succeed");

    let expected: Vec<u128> = secret_key
        .iter()
        .map(|poly| phi_polynomial(aces.channel(), poly))
        .collect();
    assert_eq!(attack.secret_image, expected);

    for _ in 0..10 {
        let message = rng.gen::<u128>() % chan.p;
        let cipher = aces.encrypt(message, &mut rng);
        let recovered =
            decrypt_with_projected_secret(aces.channel(), &cipher, &attack.secret_image);
        assert_eq!(recovered, message % chan.p);
    }
}
