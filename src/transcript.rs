//! Framed BLAKE3 Fiat-Shamir transcript for the Castor prototypes.

use std::collections::HashMap;

use p3_binary_field::TowerLevel;

use crate::{ChallengeField, merkle::update_fields};

#[derive(Clone)]
pub struct Transcript {
    state: blake3::Hasher,
}

impl Transcript {
    pub fn new(label: &[u8]) -> Self {
        let mut transcript = Self {
            state: blake3::Hasher::new(),
        };
        transcript.absorb(b"castor/transcript/v1", label);
        transcript
    }

    pub fn absorb(&mut self, label: &[u8], bytes: &[u8]) {
        self.frame(label, bytes.len());
        self.state.update(bytes);
    }

    pub(crate) fn absorb_fields<F: TowerLevel>(&mut self, label: &[u8], values: &[F]) {
        self.frame(label, values.len() * F::NUM_BYTES);
        update_fields(&mut self.state, values);
    }

    pub(crate) fn absorb_main_columns<W: TowerLevel, E: TowerLevel>(
        &mut self,
        columns: &crate::MainOpening<W, E>,
    ) {
        self.frame(b"opened-columns", columns.values_bytes());
        update_fields(&mut self.state, &columns.systematic);
        update_fields(&mut self.state, &columns.parity);
    }

    fn frame(&mut self, label: &[u8], len: usize) {
        self.state.update(&(label.len() as u64).to_le_bytes());
        self.state.update(label);
        self.state.update(&(len as u64).to_le_bytes());
    }

    fn squeeze(&mut self, label: &[u8]) -> [u8; 32] {
        self.absorb(b"sample", label);
        let output = *self.state.finalize().as_bytes();
        self.absorb(b"sample-result", &output);
        output
    }

    pub(crate) fn field(&mut self, label: &[u8]) -> ChallengeField {
        let bytes = self.squeeze(label);
        ChallengeField::from_repr(u128::from_le_bytes(bytes[..16].try_into().unwrap()))
    }

    pub(crate) fn fields(&mut self, label: &[u8], count: usize) -> Vec<ChallengeField> {
        (0..count).map(|_| self.field(label)).collect()
    }

    pub(crate) fn sorted_columns(
        &mut self,
        label: &[u8],
        count: usize,
        width: usize,
    ) -> Vec<usize> {
        let mut indices = self.columns(label, count, width);
        indices.sort_unstable();
        indices
    }

    fn below(&mut self, label: &[u8], bound: usize) -> usize {
        let bound = bound as u64;
        let threshold = bound.wrapping_neg() % bound;
        loop {
            let bytes = self.squeeze(label);
            let candidate = u64::from_le_bytes(bytes[..8].try_into().unwrap());
            if candidate >= threshold {
                return (candidate % bound) as usize;
            }
        }
    }

    /// Sample uniformly without replacement, retaining the order of the draws.
    pub(crate) fn columns(&mut self, label: &[u8], count: usize, width: usize) -> Vec<usize> {
        assert!(count <= width);
        let mut positions = HashMap::with_capacity(count);
        (width + 1 - count..=width)
            .rev()
            .map(|remaining| {
                let draw = self.below(label, remaining);
                let picked = positions.remove(&draw).unwrap_or(draw);
                let last = remaining - 1;
                if draw != last {
                    let replacement = positions.remove(&last).unwrap_or(last);
                    positions.insert(draw, replacement);
                }
                picked
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn framing_and_sampling_are_reproducible() {
        let mut a = Transcript::new(b"brakedown");
        let mut b = a.clone();
        a.absorb(b"commitment", b"test");
        b.absorb(b"commitment", b"test");
        let first = a.columns(b"columns", 8, 8);
        assert_eq!(first, b.columns(b"columns", 8, 8));
        let mut sorted = first;
        sorted.sort_unstable();
        assert_eq!(sorted, (0..8).collect::<Vec<_>>());
        assert_eq!(a.field(b"challenge"), b.field(b"challenge"));
        a.absorb(b"a", b"bc");
        b.absorb(b"ab", b"c");
        assert_ne!(a.field(b"challenge"), b.field(b"challenge"));
    }
}
