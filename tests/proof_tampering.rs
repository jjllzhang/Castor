use castor::{
    BrakedownProof, BrakedownQaPcs, CastorError, ChallengeField as F, FriConfig, MainOpening,
    MultiOpening, PcsParams, QuasarProof, QuasarQaPcs, Transcript, VerifyTimings,
};
use p3_binary_field::{BinaryField16, BinaryField64, TowerLevel};
use rand::{SeedableRng, rngs::StdRng};

// Walk the typed proof in wire order and change one byte in its containing field or digest.
struct ByteFlip(Option<usize>);

impl ByteFlip {
    fn offset_in(&mut self, len: usize) -> Option<usize> {
        let offset = self.0?;
        self.0 = offset.checked_sub(len);
        self.0.is_none().then_some(offset)
    }

    fn bytes(&mut self, bytes: &mut [u8]) {
        if let Some(offset) = self.offset_in(bytes.len()) {
            bytes[offset] ^= 1;
        }
    }

    fn fields<T: TowerLevel>(&mut self, fields: &mut [T]) {
        if let Some(offset) = self.offset_in(fields.len() * T::NUM_BYTES) {
            let field = &mut fields[offset / T::NUM_BYTES];
            let mut bytes: Vec<_> = field.into_bytes().into_iter().collect();
            bytes[offset % T::NUM_BYTES] ^= 1;
            *field = T::from_le_byte_iter(bytes.into_iter());
        }
    }

    fn main_opening<W: TowerLevel, E: TowerLevel>(&mut self, opening: &mut MainOpening<W, E>) {
        self.fields(&mut opening.systematic);
        self.fields(&mut opening.parity);
        self.bytes(opening.siblings.as_flattened_mut());
        self.bytes(opening.cap.as_flattened_mut());
    }

    fn opening<T: TowerLevel>(&mut self, opening: &mut MultiOpening<T>) {
        self.fields(&mut opening.values);
        self.bytes(opening.siblings.as_flattened_mut());
        self.bytes(opening.cap.as_flattened_mut());
    }
}

fn assert_rejects_each_serialized_byte<P: Clone>(
    proof: &P,
    serialize: impl Fn(&P) -> Vec<u8>,
    flip: impl Fn(&mut P, &mut ByteFlip),
    verify: impl Fn(&P) -> Result<VerifyTimings, CastorError>,
) {
    verify(proof).unwrap();
    let mut expected = serialize(proof);
    assert!(!expected.is_empty());
    for offset in 0..expected.len() {
        let mut changed = proof.clone();
        flip(&mut changed, &mut ByteFlip(Some(offset)));
        expected[offset] ^= 1;
        assert_eq!(serialize(&changed), expected, "wire byte {offset}");
        assert!(
            matches!(verify(&changed), Err(CastorError::InvalidProof(_))),
            "expected verification to reject the flip at wire byte {offset}"
        );
        expected[offset] ^= 1;
    }
}

fn fixture_inputs(params: PcsParams) -> (Vec<BinaryField16>, Vec<F>) {
    let witness = (0..params.num_values())
        .map(|i| BinaryField16::from_repr(131 * i as u16 + 7))
        .collect();
    let point = (0..params.ell() as u128)
        .map(|i| F::from_repr(((89 * i + 45).pow(3) << 57) | (i + 2)))
        .collect();
    (witness, point)
}

fn assert_siblings_and_cap<T>(opening: &MultiOpening<T>) {
    assert!(!opening.siblings.is_empty());
    assert!(!opening.cap.is_empty());
}

#[test]
fn brakedown_rejects_a_flip_at_every_serialized_proof_byte() {
    let params = PcsParams::new(5, 2, 2, 8, 1).unwrap();
    let pcs =
        BrakedownQaPcs::<BinaryField64>::setup(params, &mut StdRng::seed_from_u64(38)).unwrap();
    let (witness, point) = fixture_inputs(params);
    let (root, data) = pcs.commit(witness).unwrap();
    let y = pcs.evaluate(&data, &point).unwrap();
    let transcript = || Transcript::new(b"brakedown/byte-tampering");
    let (proof, _) = pcs.eval(&data, &point, y, &mut transcript()).unwrap();
    assert!(!proof.columns.systematic.is_empty());
    assert!(!proof.columns.parity.is_empty());
    assert!(!proof.columns.siblings.is_empty());
    assert!(!proof.columns.cap.is_empty());

    assert_rejects_each_serialized_byte(
        &proof,
        BrakedownProof::to_bytes,
        |proof, flip| {
            flip.fields(&mut proof.proximity_row);
            flip.fields(&mut proof.evaluation_row);
            flip.main_opening(&mut proof.columns);
        },
        |proof| pcs.verify(&root, &point, y, proof, &mut transcript()),
    );
}

#[test]
fn quasar_rejects_a_flip_at_every_serialized_proof_byte() {
    let params = PcsParams::new(7, 2, 2, 8, 1).unwrap();
    // The 2/2/1 fold schedule and rate 1/4 leave room for a sibling and cap at every layer.
    let fri = FriConfig {
        log_blowup: 2,
        fold_step: 2,
        queries: 1,
        cap_log: 1,
    };
    let pcs =
        QuasarQaPcs::<BinaryField64>::setup(params, fri, &mut StdRng::seed_from_u64(39)).unwrap();
    let (witness, point) = fixture_inputs(params);
    let (root, data) = pcs.commit(witness).unwrap();
    let transcript = || Transcript::new(b"quasar/byte-tampering");
    let (y, proof, _) = pcs.eval(&data, &point, &mut transcript()).unwrap();
    let verifier = pcs.verifier();
    assert!(!proof.columns.systematic.is_empty());
    assert!(!proof.columns.parity.is_empty());
    assert!(!proof.columns.siblings.is_empty());
    assert!(!proof.columns.cap.is_empty());
    assert_siblings_and_cap(&proof.fri.witness_packets);
    assert_siblings_and_cap(&proof.fri.parameter_packets);
    assert_eq!(proof.fri.layers.len(), 2);
    for layer in &proof.fri.layers {
        assert_siblings_and_cap(&layer.packets);
    }

    assert_rejects_each_serialized_byte(
        &proof,
        QuasarProof::to_bytes,
        |proof, flip| {
            flip.fields(&mut proof.row_evaluations);
            flip.bytes(&mut proof.spectrum_root);
            flip.main_opening(&mut proof.columns);
            for round in &mut proof.sumcheck {
                flip.fields(std::slice::from_mut(&mut round.at_zero));
                flip.fields(std::slice::from_mut(&mut round.at_omega));
                flip.fields(std::slice::from_mut(&mut round.leading));
            }
            flip.fields(&mut proof.endpoints);
            let fri = &mut proof.fri;
            flip.fields(fri.sumcheck.as_flattened_mut());
            for layer in &mut fri.layers {
                flip.bytes(&mut layer.root);
            }
            flip.fields(std::slice::from_mut(&mut fri.final_value));
            flip.opening(&mut fri.witness_packets);
            flip.opening(&mut fri.parameter_packets);
            for layer in &mut fri.layers {
                flip.opening(&mut layer.packets);
            }
        },
        |proof| verifier.verify(&root, &point, y, proof, &mut transcript()),
    );
}
