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

/// Main commitment columns in query order, with systematic entries in the witness field.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MainOpening<W, E = W> {
    pub systematic: Vec<W>,
    pub parity: Vec<E>,
    pub siblings: Vec<Digest>,
    pub cap: Vec<Digest>,
}

impl<W: TowerLevel, E: TowerLevel> MainOpening<W, E> {
    pub fn values_bytes(&self) -> usize {
        self.systematic.len() * W::NUM_BYTES + self.parity.len() * E::NUM_BYTES
    }

    pub fn authentication_bytes(&self) -> usize {
        (self.siblings.len() + self.cap.len()) * size_of::<Digest>()
    }

    pub(crate) fn write_bytes(&self, out: &mut Vec<u8>) {
        stream_fields(&self.systematic, |bytes| out.extend_from_slice(bytes));
        stream_fields(&self.parity, |bytes| out.extend_from_slice(bytes));
        out.extend_from_slice(self.siblings.as_flattened());
        out.extend_from_slice(self.cap.as_flattened());
    }
}

fn packets<F: TowerLevel>(values: &[F], width: usize, indices: &[usize], first: usize) -> Vec<F> {
    indices
        .iter()
        .flat_map(|&index| {
            let start = (index - first) * width;
            values[start..start + width].iter().copied()
        })
        .collect()
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
        Self::from_leaves(leaves)
    }

    pub(crate) fn new_split<W: TowerLevel, E: TowerLevel>(
        systematic: &[W],
        parity: &[E],
        width: usize,
    ) -> Self {
        assert_eq!(systematic.len() % width, 0);
        assert_eq!(parity.len() % width, 0);
        assert!(((systematic.len() + parity.len()) / width).is_power_of_two());
        let leaves = systematic
            .par_chunks_exact(width)
            .with_min_len(8)
            .map(leaf)
            .chain(parity.par_chunks_exact(width).with_min_len(8).map(leaf))
            .collect();
        Self::from_leaves(leaves)
    }

    fn from_leaves(leaves: Vec<Digest>) -> Self {
        let mut levels = vec![leaves];
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
        let (siblings, cap) = self.open_path(indices, cap_log);
        MultiOpening {
            values: packets(values, width, indices, 0),
            siblings,
            cap,
        }
    }

    pub(crate) fn open_split<W: TowerLevel, E: TowerLevel>(
        &self,
        systematic: &[W],
        parity: &[E],
        width: usize,
        indices: &[usize],
        cap_log: usize,
    ) -> MainOpening<W, E> {
        let systematic_cols = systematic.len() / width;
        let split = indices.partition_point(|&j| j < systematic_cols);
        let (siblings, cap) = self.open_path(indices, cap_log);
        MainOpening {
            systematic: packets(systematic, width, &indices[..split], 0),
            parity: packets(parity, width, &indices[split..], systematic_cols),
            siblings,
            cap,
        }
    }

    fn open_path(&self, indices: &[usize], cap_log: usize) -> (Vec<Digest>, Vec<Digest>) {
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
        (siblings, cap)
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
    if !valid_indices(indices, log_leaves) || opening.values.len() != indices.len() * width {
        return false;
    }
    let active = indices
        .iter()
        .copied()
        .zip(opening.values.chunks_exact(width).map(leaf))
        .collect();
    verify_path(
        root,
        active,
        &opening.siblings,
        &opening.cap,
        log_leaves,
        cap_log,
    )
}

pub(crate) fn verify_split<W: TowerLevel, E: TowerLevel>(
    root: &Digest,
    indices: &[usize],
    opening: &MainOpening<W, E>,
    width: usize,
    systematic_cols: usize,
    log_leaves: usize,
    cap_log: usize,
) -> bool {
    let split = indices.partition_point(|&j| j < systematic_cols);
    if !valid_indices(indices, log_leaves)
        || opening.systematic.len() != split * width
        || opening.parity.len() != (indices.len() - split) * width
    {
        return false;
    }
    let leaves = opening
        .systematic
        .chunks_exact(width)
        .map(leaf)
        .chain(opening.parity.chunks_exact(width).map(leaf));
    let active = indices.iter().copied().zip(leaves).collect();
    verify_path(
        root,
        active,
        &opening.siblings,
        &opening.cap,
        log_leaves,
        cap_log,
    )
}

fn valid_indices(indices: &[usize], log_leaves: usize) -> bool {
    indices.last().is_some_and(|&last| last < 1 << log_leaves)
        && indices.windows(2).all(|pair| pair[0] < pair[1])
}

fn verify_path(
    root: &Digest,
    mut active: Vec<(usize, Digest)>,
    siblings: &[Digest],
    cap: &[Digest],
    log_leaves: usize,
    cap_log: usize,
) -> bool {
    if cap.len() != if cap_log == 0 { 0 } else { 1 << cap_log } {
        return false;
    }
    let mut siblings = siblings.iter();
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
    if active.iter().any(|&(index, digest)| cap[index] != digest) {
        return false;
    }
    let mut level = parents(cap);
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
                // Compare paths truncated at the same height, with the shared cap sent once.
                let independent = proof.values_bytes()
                    + (indices.len() * (3 - cap) + proof.cap.len()) * size_of::<Digest>();
                assert!(bytes.len() <= independent);
                if indices.len() == 1 {
                    assert_eq!(bytes.len(), independent);
                } else if cap == 0 {
                    assert!(bytes.len() < independent);
                }

                for i in 0..proof.values.len() {
                    let value = proof.values[i];
                    proof.values[i] += BinaryField64::from_repr(1);
                    assert!(!verify(&tree.root(), &indices, &proof, 4, 3, cap));
                    proof.values[i] = value;
                }
                for i in 0..proof.siblings.len() {
                    proof.siblings[i][0] ^= 1;
                    assert!(!verify(&tree.root(), &indices, &proof, 4, 3, cap));
                    proof.siblings[i][0] ^= 1;
                }
                for i in 0..proof.cap.len() {
                    proof.cap[i][0] ^= 1;
                    assert!(!verify(&tree.root(), &indices, &proof, 4, 3, cap));
                    proof.cap[i][0] ^= 1;
                }
            }
        }
    }
}
