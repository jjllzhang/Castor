use castor::{BrakedownQaPcs, ChallengeField as F, FriConfig, PcsParams, QuasarQaPcs, Transcript};
use p3_binary_dft::EncodableLevel;
use p3_binary_field::{BinaryField64, TowerLevel};
use p3_binary_qa::QaField;
use p3_field::Algebra;
use rand::{SeedableRng, rngs::StdRng};

fn assert_legacy<T: QaField + EncodableLevel>(root_hex: &str, brakedown_hex: &str, quasar_hex: &str)
where
    F: Algebra<T>,
{
    let params = PcsParams::new(7, 2, 2, 20, 1).unwrap();
    let fri = FriConfig {
        log_blowup: 2,
        fold_step: 2,
        queries: 4,
        cap_log: 1,
    };
    let witness: Vec<T> = (0u64..128)
        .map(|i| T::from_le_byte_iter((131 * i + 7).to_le_bytes().into_iter().chain([0; 8])))
        .collect();
    let point: Vec<_> = (0..7u128)
        .map(|i| F::from_repr((89 * i + 45).pow(3) << 57))
        .collect();
    let transcript = || Transcript::new(b"castor/v1-compatibility");
    let brakedown = BrakedownQaPcs::<T>::setup(params, &mut StdRng::seed_from_u64(203)).unwrap();
    let quasar = QuasarQaPcs::<T>::setup(params, fri, &mut StdRng::seed_from_u64(203)).unwrap();
    let (root, data) = brakedown.commit_full(witness.clone()).unwrap();
    let y = brakedown.evaluate(&data, &point).unwrap();
    let (proof, _) = brakedown.eval(&data, &point, y, &mut transcript()).unwrap();
    brakedown
        .verify(&root, &point, y, &proof, &mut transcript())
        .unwrap();
    assert_eq!(root, *blake3::Hash::from_hex(root_hex).unwrap().as_bytes());
    assert_eq!(
        blake3::hash(&proof.to_bytes()),
        blake3::Hash::from_hex(brakedown_hex).unwrap()
    );
    let (root, data) = quasar.commit_full(witness).unwrap();
    let (y, proof, _) = quasar.eval(&data, &point, &mut transcript()).unwrap();
    quasar
        .verify(&root, &point, y, &proof, &mut transcript())
        .unwrap();
    assert_eq!(root, *blake3::Hash::from_hex(root_hex).unwrap().as_bytes());
    assert_eq!(
        blake3::hash(&proof.to_bytes()),
        blake3::Hash::from_hex(quasar_hex).unwrap()
    );
}

#[test]
fn full_width_commitments_and_proofs_match_v1() {
    // Roots and proof hashes generated with Castor e501780 on these same fixtures.
    assert_legacy::<BinaryField64>(
        "cfd6ee81faa06986fe0f25a4148fa6f13c8c97b74d0adeee14be56b583dec3e8",
        "e2e5a8a8df24460e8bd85dca304efe2339b7a86126542060e135d7f7058d4aac",
        "469ed80c807d88b9e6650d7cf828074a3927bfc3c9f42343a8fab3ed4a849f88",
    );
    assert_legacy::<F>(
        "da73721ba1c504e160b695e17c1ca12a72e2cbd29ffdeb25c128c4238b409a35",
        "47de7dbf2cbfd551c0fae4a191cb1b8c95cad051ad11ca77950a0ece2507d7c3",
        "6ed10122355c216946da4187635f8b5af54dc7c7b99698c952d481ab024750bf",
    );
}
