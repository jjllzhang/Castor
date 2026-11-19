//! Multilinear tables in increasing index order, with the first variable in the low bit.

use p3_field::{Algebra, PrimeCharacteristicRing};
use p3_maybe_rayon::prelude::*;

use crate::ChallengeField as F;

pub(crate) fn eq_table(point: &[F]) -> Vec<F> {
    let mut weights = Vec::with_capacity(1 << point.len());
    weights.push(F::ONE);
    for &coordinate in point {
        let width = weights.len();
        weights.resize(2 * width, F::ZERO);
        let (low, high) = weights.split_at_mut(width);
        for (lo, hi) in low.iter_mut().zip(high) {
            *hi = *lo * coordinate;
            *lo += *hi;
        }
    }
    weights
}

#[cfg(test)]
pub(crate) fn eq_index(point: &[F], index: usize) -> F {
    point
        .iter()
        .enumerate()
        .map(|(var, &x)| {
            if (index >> var) & 1 == 0 {
                F::ONE + x
            } else {
                x
            }
        })
        .product()
}

pub(crate) fn eq_points(left: &[F], right: &[F]) -> F {
    assert_eq!(left.len(), right.len());
    left.iter()
        .zip(right)
        .map(|(&x, &y)| F::ONE + x + y)
        .product()
}

pub(crate) fn dot<S: Copy>(weights: &[F], values: &[S]) -> F
where
    F: Algebra<S>,
{
    assert_eq!(weights.len(), values.len());
    weights.iter().zip(values).map(|(&w, &x)| w * x).sum()
}

/// Fold the contiguous columns of a column-interleaved matrix.
pub(crate) fn fold_row<S: Copy + Sync>(matrix: &[S], weights: &[F]) -> Vec<F>
where
    F: Algebra<S>,
{
    matrix
        .par_chunks_exact(weights.len())
        .with_min_len(32)
        .map(|column| dot(weights, column))
        .collect()
}

/// Apply several row combinations in one pass over the witness.
pub(crate) fn fold_rows<S: Copy + Sync>(matrix: &[S], rho: &[F], chi: &[F]) -> [Vec<F>; 2]
where
    F: Algebra<S>,
{
    let height = rho.len();
    let width = matrix.len() / height;
    let mut v = F::zero_vec(width);
    let mut p = F::zero_vec(width);
    v.par_iter_mut()
        .zip(p.par_iter_mut())
        .zip(matrix.par_chunks_exact(height))
        .with_min_len(32)
        .for_each(|((v, p), column)| {
            for ((&x, &a), &b) in column.iter().zip(rho).zip(chi) {
                *v += a * x;
                *p += b * x;
            }
        });
    [v, p]
}

/// Evaluate each message row at the column point while reading contiguous columns.
pub(crate) fn evaluate_rows<S: Copy + Sync>(
    matrix: &[S],
    height: usize,
    column_weights: &[F],
) -> Vec<F>
where
    F: Algebra<S>,
{
    assert_eq!(matrix.len(), height * column_weights.len());
    const TILE: usize = 64;
    matrix
        .par_chunks(TILE * height)
        .zip(column_weights.par_chunks(TILE))
        .with_min_len(8)
        .par_fold_reduce(
            || vec![F::ZERO; height],
            |mut sums, (tile, weights)| {
                for (column, &weight) in tile.chunks_exact(height).zip(weights) {
                    for (sum, &value) in sums.iter_mut().zip(column) {
                        *sum += weight * value;
                    }
                }
                sums
            },
            |mut left, right| {
                for (a, b) in left.iter_mut().zip(right) {
                    *a += b;
                }
                left
            },
        )
}

pub(crate) fn evaluate<S: Copy + Sync>(matrix: &[S], row: &[F], column: &[F]) -> F
where
    F: Algebra<S>,
{
    assert_eq!(matrix.len(), row.len() * column.len());
    matrix
        .par_chunks_exact(row.len())
        .zip(column.par_iter())
        .with_min_len(32)
        .map(|(values, &weight)| weight * dot(row, values))
        .sum()
}

#[cfg(test)]
mod tests {
    use p3_binary_field::{BinaryField16, TowerLevel};

    use super::*;

    #[test]
    fn equality_tables_match_scalar_products() {
        let point = [F::from_repr(2), F::from_repr(9), F::from_repr(11)];
        let weights = eq_table(&point);
        assert_eq!(weights.len(), 8);
        for (index, weight) in weights.into_iter().enumerate() {
            assert_eq!(weight, eq_index(&point, index));
            let boolean: Vec<_> = (0..3)
                .map(|i| F::from_bool((index >> i) & 1 == 1))
                .collect();
            assert_eq!(weight, eq_points(&point, &boolean));
        }
        let [x, y, _] = point;
        assert_eq!(eq_points(&[x], &[y]), x * y + (F::ONE + x) * (F::ONE + y));
    }

    #[test]
    fn mixed_fold_matches_naive() {
        let witness: Vec<_> = (0..32).map(|i| BinaryField16::from_repr(17 * i)).collect();
        let z = [F::from_repr(29), F::from_repr(13)];
        let row = eq_table(&z);
        let reverse: Vec<_> = row.iter().copied().rev().collect();
        let column = eq_table(&[F::from_repr(12), F::from_repr(3), F::from_repr(5)]);
        let expected: Vec<F> = witness
            .chunks_exact(4)
            .map(|col| col.iter().zip(&row).map(|(&x, &r)| F::from(x) * r).sum())
            .collect();
        assert_eq!(fold_row(&witness, &row), expected);
        assert_eq!(
            fold_rows(&witness, &row, &reverse),
            [expected.clone(), fold_row(&witness, &reverse)]
        );
        assert_eq!(evaluate(&witness, &row, &column), dot(&column, &expected));
    }
}
