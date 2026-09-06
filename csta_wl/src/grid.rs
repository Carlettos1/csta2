use crate::{Error, Result};

/// Shared sampling and thermodynamic energy coordinates.
/// Interval DOS values are bin masses, not densities per unit energy.
#[derive(Clone, Debug, PartialEq)]
pub struct EnergyGrid {
    energies: Vec<f64>,
    edges: Option<Vec<f64>>,
}

impl EnergyGrid {
    /// Half-open bins, except the final bin includes `max`.
    pub fn continuous(min: f64, max: f64, bins: usize) -> Result<Self> {
        if bins == 0 || !min.is_finite() || !max.is_finite() || min >= max {
            return Err(Error::Invalid("nonempty grid needs finite min < max"));
        }
        let width = (max - min) / bins as f64;
        if !width.is_finite() || width <= 0.0 {
            return Err(Error::Invalid("bin width must be finite and positive"));
        }
        if bins
            .checked_add(1)
            .is_none_or(|n| n > isize::MAX as usize / std::mem::size_of::<f64>())
        {
            return Err(Error::Invalid(
                "grid allocation exceeds addressable capacity",
            ));
        }
        let edges: Vec<_> = (0..=bins)
            .map(|i| {
                if i == bins {
                    max
                } else {
                    min + i as f64 * width
                }
            })
            .collect();
        if edges.windows(2).any(|v| v[0] >= v[1]) {
            return Err(Error::Invalid("bin edges are not representable"));
        }
        let energies: Vec<_> = edges
            .windows(2)
            .map(|v| v[0] + (v[1] - v[0]) / 2.0)
            .collect();
        let mut grid = Self::discrete(energies)?;
        grid.edges = Some(edges);
        Ok(grid)
    }

    /// Exact finite, strictly increasing levels. Off-grid energies are rejected;
    /// models should compute energies in the same units as these levels.
    pub fn discrete(energies: Vec<f64>) -> Result<Self> {
        if energies.is_empty()
            || energies.iter().any(|e| !e.is_finite())
            || energies
                .windows(2)
                .any(|v| v[0] >= v[1] || !(v[1] - v[0]).is_finite())
        {
            return Err(Error::Invalid(
                "energy levels must be finite, nonempty and strictly increasing",
            ));
        }
        Ok(Self {
            energies,
            edges: None,
        })
    }

    pub fn energies(&self) -> &[f64] {
        &self.energies
    }
    pub fn len(&self) -> usize {
        self.energies.len()
    }
    pub fn is_empty(&self) -> bool {
        self.energies.is_empty()
    }
    pub fn bin(&self, energy: f64) -> Result<Option<usize>> {
        if !energy.is_finite() {
            return Err(Error::InvalidEnergy);
        }
        if let Some(edges) = &self.edges {
            if energy < edges[0] || energy > edges[edges.len() - 1] {
                return Ok(None);
            }
            Ok(Some(
                edges
                    .partition_point(|e| *e <= energy)
                    .saturating_sub(1)
                    .min(self.len() - 1),
            ))
        } else {
            Ok(self
                .energies
                .binary_search_by(|e| e.partial_cmp(&energy).unwrap())
                .ok())
        }
    }
}
