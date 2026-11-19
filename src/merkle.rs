//! BLAKE3 commitments to fixed-width packets with shared authentication paths.

use p3_binary_field::TowerLevel;
use p3_maybe_rayon::prelude::*;

pub type Digest = [u8; 32];

pub(crate) fn stream_fields<F: TowerLevel>(fields: &[F], mut write: impl FnMut(&[u8])) {
    let width = F::NUM_BYTES;
    let mut buffer = [0_u8; 1024];
    for block in fields.chunks(buffer.len() / width) {
        for (slot, &field) in buffer.chunks_exact_mut(width).zip(block) {
            for (byte, value) in slot.iter_mut().zip(field.into_bytes()) {
                *byte = value;
            }
        }
        write(&buffer[..block.len() * width]);
    }
}

pub(crate) fn update_fields<F: TowerLevel>(hasher: &mut blake3::Hasher, fields: &[F]) {
    stream_fields(fields, |bytes| {
        hasher.update(bytes);
    });
}

fn leaf<F: TowerLevel>(packet: &[F]) -> Digest {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"castor/packet/v1");
    hasher.update(&(F::NUM_BYTES as u64).to_le_bytes());
    hasher.update(&(packet.len() as u64).to_le_bytes());
    update_fields(&mut hasher, packet);
    *hasher.finalize().as_bytes()
}

fn node(left: &Digest, right: &Digest) -> Digest {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"castor/node/v1");
    hasher.update(left);
    hasher.update(right);
    *hasher.finalize().as_bytes()
}

fn parents(level: &[Digest]) -> Vec<Digest> {
    level
        .par_chunks_exact(2)
        .with_min_len(64)
        .map(|pair| node(&pair[0], &pair[1]))
        .collect()
}

pub(crate) struct MerkleTree {
    levels: Vec<Vec<Digest>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MultiOpening<F> {
    pub values: Vec<F>,
    pub siblings: Vec<Digest>,
    pub cap: Vec<Digest>,
}

impl<F: TowerLevel> MultiOpening<F> {
    pub fn values_bytes(&self) -> usize {
        self.values.len() * F::NUM_BYTES
    }

    pub fn authentication_bytes(&self) -> usize {
        (self.siblings.len() + self.cap.len()) * size_of::<Digest>()
    }

    pub(crate) fn write_bytes(&self, out: &mut Vec<u8>) {
        stream_fields(&self.values, |bytes| out.extend_from_slice(bytes));
        out.extend_from_slice(self.siblings.as_flattened());
        out.extend_from_slice(self.cap.as_flattened());
    }
}

impl MerkleTree {
    pub(crate) fn new<F: TowerLevel>(values: &[F], width: usize) -> Self {
        assert_eq!(values.len() % width, 0);
        assert!((values.len() / width).is_power_of_two());
        let leaves = values
            .par_chunks_exact(width)
            .with_min_len(8)
            .map(leaf)
            .collect();
        let mut levels: Vec<Vec<Digest>> = vec![leaves];
        while levels.last().unwrap().len() > 1 {
            levels.push(parents(levels.last().unwrap()));
        }
        Self { levels }
    }

    pub(crate) fn root(&self) -> Digest {
        self.levels.last().unwrap()[0]
    }

    /// Indices are unique and increasing. No positions are transmitted on the wire.
    pub(crate) fn open<F: TowerLevel>(
        &self,
        values: &[F],
        width: usize,
        indices: &[usize],
        cap_log: usize,
    ) -> MultiOpening<F> {
        let stop_level = self.levels.len() - 1 - cap_log;
        let mut active = indices.to_vec();
        let mut siblings = Vec::new();
        for level in &self.levels[..stop_level] {
            let mut next = Vec::with_capacity(active.len());
            let mut cursor = 0;
            while cursor < active.len() {
                let index = active[cursor];
                cursor += 1;
                if index & 1 == 0 && active.get(cursor) == Some(&(index + 1)) {
                    cursor += 1;
                } else {
                    siblings.push(level[index ^ 1]);
                }
                next.push(index >> 1);
            }
            active = next;
        }
        let cap = if cap_log == 0 {
            Vec::new()
        } else {
            self.levels[stop_level].clone()
        };
        MultiOpening {
            values: indices
                .iter()
                .flat_map(|&index| values[index * width..(index + 1) * width].iter().copied())
                .collect(),
            siblings,
            cap,
        }
    }
}

pub(crate) fn verify<F: TowerLevel>(
    root: &Digest,
    indices: &[usize],
    opening: &MultiOpening<F>,
    width: usize,
    log_leaves: usize,
    cap_log: usize,
) -> bool {
    if indices.is_empty()
        || indices.windows(2).any(|w| w[0] >= w[1])
        || *indices.last().unwrap() >= 1 << log_leaves
        || opening.values.len() != indices.len() * width
        || opening.cap.len() != if cap_log == 0 { 0 } else { 1 << cap_log }
    {
        return false;
    }
    let mut active: Vec<_> = indices
        .iter()
        .copied()
        .zip(opening.values.chunks_exact(width).map(leaf))
        .collect();
    let mut siblings = opening.siblings.iter();
    for _ in 0..log_leaves - cap_log {
        let mut next = Vec::with_capacity(active.len());
        let mut items = active.into_iter().peekable();
        while let Some((index, digest)) = items.next() {
            let sibling =
                if index & 1 == 0 && items.peek().is_some_and(|&(other, _)| other == index + 1) {
                    items.next().unwrap().1
                } else if let Some(&sibling) = siblings.next() {
                    sibling
                } else {
                    return false;
                };
            let parent = if index & 1 == 0 {
                node(&digest, &sibling)
            } else {
                node(&sibling, &digest)
            };
            next.push((index >> 1, parent));
        }
        active = next;
    }
    if siblings.next().is_some() {
        return false;
    }
    if cap_log == 0 {
        return active == [(0, *root)];
    }
    if active
        .iter()
        .any(|&(index, digest)| opening.cap[index] != digest)
    {
        return false;
    }
    let mut level = opening.cap.clone();
    while level.len() > 1 {
        level = parents(&level);
    }
    level[0] == *root
}

#[cfg(test)]
mod tests {
    use super::*;
    use p3_binary_field::BinaryField64;

    #[test]
    fn opens_every_nonempty_subset_with_or_without_cap() {
        let values: Vec<_> = (0..32).map(BinaryField64::from_repr).collect();
        let tree = MerkleTree::new(&values, 4);
        for mask in 1_usize..256 {
            let indices: Vec<_> = (0..8).filter(|i| mask & (1 << i) != 0).collect();
            for cap in 0..=3 {
                let mut proof = tree.open(&values, 4, &indices, cap);
                assert!(verify(&tree.root(), &indices, &proof, 4, 3, cap));
                assert_eq!(proof.values_bytes(), indices.len() * 32);
                let mut bytes = Vec::new();
                proof.write_bytes(&mut bytes);
                assert_eq!(
                    bytes.len(),
                    proof.values_bytes() + proof.authentication_bytes()
                );
                if cap == 0 && indices.len() > 1 {
                    assert!(proof.authentication_bytes() < indices.len() * 3 * 32);
                }
                if let Some(digest) = proof.siblings.first_mut().or_else(|| proof.cap.first_mut()) {
                    digest[0] ^= 1;
                    assert!(!verify(&tree.root(), &indices, &proof, 4, 3, cap));
                    proof = tree.open(&values, 4, &indices, cap);
                }
                proof.values[0] += BinaryField64::from_repr(1);
                assert!(!verify(&tree.root(), &indices, &proof, 4, 3, cap));
            }
        }
    }
}
