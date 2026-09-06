//! Small analysis helpers for finite-size recipes; no plotting dependency.
use crate::{Error, Result};
use rand::RngExt;
fn curve(data: &[(f64, f64)]) -> Result<()> {
    if data.len() < 2
        || data.iter().any(|(x, y)| !x.is_finite() || !y.is_finite())
        || data.windows(2).any(|p| p[0].0 >= p[1].0)
    {
        Err(Error("curve needs increasing finite coordinates"))
    } else {
        Ok(())
    }
}
/// All linear-interpolated crossings on a shared grid. Empty/multiple results
/// remain explicit. Coincident segments are ambiguous and return an error.
pub fn crossings(a: &[(f64, f64)], b: &[(f64, f64)]) -> Result<Vec<f64>> {
    curve(a)?;
    curve(b)?;
    if a.len() != b.len() || a.iter().zip(b).any(|(a, b)| a.0 != b.0) {
        return Err(Error("crossing grids differ"));
    }
    let mut out = Vec::new();
    for i in 0..a.len() - 1 {
        let d0 = a[i].1 - b[i].1;
        let d1 = a[i + 1].1 - b[i + 1].1;
        if !d0.is_finite() || !d1.is_finite() {
            return Err(Error("curve difference overflow"));
        }
        if d0 == 0.0 && d1 == 0.0 {
            return Err(Error("coincident curves"));
        }
        if d0 == 0.0 {
            out.push(a[i].0);
        } else if d0.is_sign_positive() != d1.is_sign_positive() && d1 != 0.0 {
            let scale = d0.abs().max(d1.abs());
            let f = (d0 / scale) / ((d0 / scale) - (d1 / scale));
            out.push((1.0 - f) * a[i].0 + f * a[i + 1].0);
        }
    }
    if a.last().unwrap().1 == b.last().unwrap().1 {
        out.push(a.last().unwrap().0);
    }
    Ok(out)
}
/// An interior unique maximum on the supplied grid; no unreported extrapolation.
pub fn peak(data: &[(f64, f64)]) -> Result<Option<(f64, f64)>> {
    curve(data)?;
    let max = data.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
    let indexes: Vec<_> = data
        .iter()
        .enumerate()
        .filter(|(_, p)| p.1 == max)
        .map(|(i, _)| i)
        .collect();
    Ok(
        if indexes.len() == 1 && indexes[0] > 0 && indexes[0] + 1 < data.len() {
            Some(data[indexes[0]])
        } else {
            None
        },
    )
}
/// Nonoverlapping block bootstrap of an arbitrary scalar statistic. Resamples
/// complete blocks; excludes trailing incomplete samples. Caller chooses blocks
/// longer than the correlation scale. Returns the bootstrap standard deviation.
pub fn block_bootstrap(
    values: &[f64],
    block: usize,
    repetitions: usize,
    rng: &mut impl RngExt,
    statistic: impl Fn(&[f64]) -> Result<f64>,
) -> Result<f64> {
    if block == 0
        || repetitions < 2
        || values.len() / block < 8
        || values.iter().any(|x| !x.is_finite())
    {
        return Err(Error("insufficient or invalid bootstrap data"));
    }
    let n = values.len() / block;
    let mut moments = crate::statistics::Moments::default();
    let mut resampled = Vec::with_capacity(n * block);
    for _ in 0..repetitions {
        resampled.clear();
        for _ in 0..n {
            let j = rng.random_range(0..n);
            resampled.extend_from_slice(&values[j * block..(j + 1) * block]);
        }
        moments.push(statistic(&resampled)?)?;
    }
    Ok(moments.variance().unwrap().sqrt())
}
#[cfg(test)]
mod tests {
    use super::*;
    use rand::{SeedableRng, rngs::StdRng};
    #[test]
    fn curves_and_bootstrap() {
        let a = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)];
        let b = [(0.0, 0.5), (1.0, 0.5), (2.0, 0.5)];
        assert_eq!(crossings(&a, &b).unwrap(), vec![0.5, 1.5]);
        assert!(crossings(&a, &a).is_err());
        assert_eq!(peak(&a).unwrap(), Some((1.0, 1.0)));
        assert!(peak(&b).unwrap().is_none());
        assert!(crossings(&a, &b[..2]).is_err());
        assert!(crossings(&[(0.0, 0.0)], &[(0.0, 0.0)]).is_err());
        let values: Vec<_> = (0..100).map(|i| i as f64).collect();
        let mut rng = StdRng::seed_from_u64(3);
        let se = block_bootstrap(&values, 1, 500, &mut rng, |x| {
            Ok(x.iter().sum::<f64>() / x.len() as f64)
        })
        .unwrap();
        assert!((se - 2.9).abs() < 0.4);
        assert!(block_bootstrap(&values, 0, 500, &mut rng, |_| Ok(0.0)).is_err());
    }
}
#[cfg(test)]
mod range_tests {
    use super::*;
    #[test]
    fn ranges_without_unwarranted_crossings() {
        let a = [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)];
        let b = [(0.0, 1.5), (1.0, 1.5), (2.0, 1.5)];
        assert_eq!(crossings(&a, &b).unwrap(), vec![1.5]);
        assert!(crossings(&a[..2], &b[..2]).unwrap().is_empty());
        assert!(peak(&a).unwrap().is_none());
        assert!(crossings(&[(0.0, f64::NAN), (1.0, 1.0)], &a[..2]).is_err());
    }
}
