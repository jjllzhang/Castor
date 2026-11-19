use std::{hint::black_box, time::Duration};

use castor::{
    BrakedownQaPcs, ChallengeField as F, FriConfig, PcsParams, ProximityBoundMode, ProximityBudget,
    QuasarQaPcs, Transcript,
};
use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};
use p3_binary_field::{BinaryField16, BinaryField32, BinaryField64, TowerLevel};
use p3_binary_qa::QaField;
use p3_field::Algebra;
use rand::{Rng, SeedableRng, rngs::StdRng};

const SEED: u64 = 20261119;

fn fixture<W: QaField>(m: usize) -> (Vec<W>, Vec<F>) {
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut bytes = [0; 16];
    let witness = (0..1 << m)
        .map(|_| {
            rng.fill_bytes(&mut bytes[..W::NUM_BYTES]);
            W::from_le_byte_iter(bytes.into_iter())
        })
        .collect();
    let point = (0..m)
        .map(|_| {
            rng.fill_bytes(&mut bytes);
            F::from_repr(u128::from_le_bytes(bytes))
        })
        .collect();
    (witness, point)
}

fn params(m: usize, r: usize, c: usize) -> PcsParams {
    let budget = ProximityBudget::derive(
        m - r,
        c,
        128,
        ProximityBoundMode::UniqueDecodingCA,
        103,
        None,
    )
    .unwrap();
    PcsParams::new(m, r, c.ilog2() as usize, budget.queries, 0).unwrap()
}

fn transcript() -> Transcript {
    Transcript::new(b"castor/criterion")
}

fn brakedown(criterion: &mut Criterion, m: usize, c: usize) {
    let (witness, point) = fixture::<BinaryField32>(m);
    let pcs = BrakedownQaPcs::<F>::setup(params(m, m / 2, c), &mut StdRng::seed_from_u64(SEED ^ 1))
        .unwrap();
    let (root, data) = pcs.commit(witness.clone()).unwrap();
    let value = pcs.evaluate(&data, &point).unwrap();
    let (proof, _) = pcs.eval(&data, &point, value, &mut transcript()).unwrap();
    pcs.verify(&root, &point, value, &proof, &mut transcript())
        .unwrap();
    let mut group = criterion.benchmark_group(format!("brakedown/m{m}/c{c}/bf32"));
    group.throughput(Throughput::Elements(1 << m));
    group.bench_function("commit", |b| {
        b.iter_batched(
            || witness.clone(),
            |table| black_box(pcs.commit(table).unwrap()),
            BatchSize::LargeInput,
        )
    });
    group.bench_function("eval", |b| {
        b.iter(|| black_box(pcs.eval(&data, &point, value, &mut transcript()).unwrap()))
    });
    group.bench_function("verify", |b| {
        b.iter(|| {
            black_box(
                pcs.verify(&root, &point, value, &proof, &mut transcript())
                    .unwrap(),
            )
        })
    });
    group.finish();
}

fn quasar<W: QaField>(criterion: &mut Criterion, m: usize, c: usize)
where
    F: p3_field::ExtensionField<W> + Algebra<W>,
{
    let (witness, point) = fixture::<W>(m);
    let pcs = QuasarQaPcs::<F>::setup(
        params(m, 5, c),
        FriConfig::default(),
        &mut StdRng::seed_from_u64(SEED ^ 1),
    )
    .unwrap();
    let (root, data) = pcs.commit::<W>(witness.clone()).unwrap();
    let (value, proof, _) = pcs.eval::<W>(&data, &point, &mut transcript()).unwrap();
    let verifier = pcs.verifier();
    verifier
        .verify(&root, &point, value, &proof, &mut transcript())
        .unwrap();
    let mut group = criterion.benchmark_group(format!("quasar/m{m}/c{c}/bf{}", W::NUM_BYTES * 8));
    group.throughput(Throughput::Elements(1 << m));
    group.bench_function("commit", |b| {
        b.iter_batched(
            || witness.clone(),
            |table| black_box(pcs.commit::<W>(table).unwrap()),
            BatchSize::LargeInput,
        )
    });
    group.bench_function("eval", |b| {
        b.iter(|| black_box(pcs.eval::<W>(&data, &point, &mut transcript()).unwrap()))
    });
    group.bench_function("verify", |b| {
        b.iter(|| {
            black_box(
                verifier
                    .verify(&root, &point, value, &proof, &mut transcript())
                    .unwrap(),
            )
        })
    });
    group.finish();
}

fn pcs_comparison(criterion: &mut Criterion) {
    for m in [16, 20] {
        for c in [2, 4] {
            brakedown(criterion, m, c);
            quasar::<BinaryField32>(criterion, m, c);
        }
    }
    quasar::<BinaryField16>(criterion, 20, 2);
    quasar::<BinaryField64>(criterion, 20, 2);
    quasar::<F>(criterion, 20, 2);
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(2));
    targets = pcs_comparison
}
criterion_main!(benches);
