//! Orthorhombic periodic geometry. Half-box ties map to the negative side.
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[cfg_attr(
    feature = "serde",
    serde(bound(
        serialize = "[f64; D]: serde::Serialize",
        deserialize = "[f64; D]: serde::Deserialize<'de>"
    ))
)]
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicBox<const D: usize> {
    lengths: [f64; D],
}
impl<const D: usize> PeriodicBox<D> {
    pub fn new(lengths: [f64; D]) -> Result<Self, &'static str> {
        if D == 0 || lengths.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err("box lengths must be finite and positive");
        }
        Ok(Self { lengths })
    }
    pub fn lengths(&self) -> &[f64; D] {
        &self.lengths
    }
    pub fn wrap(&self, position: [f64; D]) -> Result<[f64; D], &'static str> {
        Self::new(self.lengths)?;
        if position.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite position");
        }
        Ok(std::array::from_fn(|i| {
            let x = position[i].rem_euclid(self.lengths[i]);
            if x == self.lengths[i] { 0.0 } else { x }
        }))
    }
    pub fn displacement(&self, from: [f64; D], to: [f64; D]) -> Result<[f64; D], &'static str> {
        let a = self.wrap(from)?;
        let b = self.wrap(to)?;
        Ok(std::array::from_fn(|i| {
            let mut d = b[i] - a[i];
            let half = self.lengths[i] * 0.5;
            if d >= half {
                d -= self.lengths[i];
            } else if d < -half {
                d += self.lengths[i];
            }
            d
        }))
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
    Open,
    Periodic,
}

/// Rectangular lattice; periodic dimensions must be >= 3 to avoid ambiguous
/// multiple bonds. Each undirected nearest-neighbor bond appears once in `bonds`.
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[derive(Clone, Debug)]
pub struct Lattice2D {
    width: usize,
    height: usize,
    boundary: Boundary,
}
impl Lattice2D {
    pub fn new(width: usize, height: usize, boundary: Boundary) -> Result<Self, &'static str> {
        let n = width.checked_mul(height).ok_or("lattice size overflow")?;
        if n == 0
            || n > isize::MAX as usize
            || (boundary == Boundary::Periodic && (width < 3 || height < 3))
        {
            return Err("invalid lattice dimensions");
        }
        Ok(Self {
            width,
            height,
            boundary,
        })
    }
    pub fn valid(&self) -> bool {
        Self::new(self.width, self.height, self.boundary).is_ok()
    }
    pub fn len(&self) -> usize {
        self.width * self.height
    }
    pub fn is_empty(&self) -> bool {
        false
    }
    pub fn neighbors(&self, site: usize) -> Result<Vec<usize>, &'static str> {
        if site >= self.len() {
            return Err("site outside lattice");
        }
        let (x, y) = (site % self.width, site / self.width);
        let mut out = Vec::with_capacity(4);
        if x > 0 {
            out.push(site - 1);
        } else if self.boundary == Boundary::Periodic {
            out.push(site + self.width - 1);
        }
        if x + 1 < self.width {
            out.push(site + 1);
        } else if self.boundary == Boundary::Periodic {
            out.push(site + 1 - self.width);
        }
        if y > 0 {
            out.push(site - self.width);
        } else if self.boundary == Boundary::Periodic {
            out.push(site + (self.height - 1) * self.width);
        }
        if y + 1 < self.height {
            out.push(site + self.width);
        } else if self.boundary == Boundary::Periodic {
            out.push(x);
        }
        Ok(out)
    }
    pub fn bonds(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.len()).flat_map(|i| {
            self.neighbors(i)
                .expect("valid site")
                .into_iter()
                .filter(move |j| i < *j)
                .map(move |j| (i, j))
        })
    }
}

#[cfg(feature = "serde")]
impl<'de, const D: usize> serde::Deserialize<'de> for PeriodicBox<D>
where
    [f64; D]: serde::Deserialize<'de>,
{
    fn deserialize<T: serde::Deserializer<'de>>(d: T) -> Result<Self, T::Error> {
        #[derive(serde::Deserialize)]
        #[serde(bound(deserialize = "[f64; D]: serde::Deserialize<'de>"))]
        struct Raw<const D: usize> {
            lengths: [f64; D],
        }
        let raw = Raw::<D>::deserialize(d)?;
        Self::new(raw.lengths).map_err(serde::de::Error::custom)
    }
}
#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for Lattice2D {
    fn deserialize<T: serde::Deserializer<'de>>(d: T) -> Result<Self, T::Error> {
        #[derive(serde::Deserialize)]
        struct Raw {
            width: usize,
            height: usize,
            boundary: Boundary,
        }
        let raw = Raw::deserialize(d)?;
        Self::new(raw.width, raw.height, raw.boundary).map_err(serde::de::Error::custom)
    }
}
