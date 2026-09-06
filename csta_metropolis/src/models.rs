//! Small validated models. Hamiltonians count each undirected bond once.
use crate::{Error, Result, State};
pub use csta_core::geometry::Boundary;
use csta_core::geometry::Lattice2D;
use rand::RngExt;

#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug)]
pub struct Ising2D {
    lattice: Lattice2D,
    spins: Vec<i8>,
    j: f64,
    h: f64,
    energy: f64,
}
impl Ising2D {
    pub fn new(
        width: usize,
        height: usize,
        boundary: Boundary,
        j: f64,
        h: f64,
        spins: Vec<i8>,
    ) -> Result<Self> {
        let lattice = Lattice2D::new(width, height, boundary).map_err(Error)?;
        if spins.len() != lattice.len()
            || spins.iter().any(|s| *s != 1 && *s != -1)
            || !j.is_finite()
            || !h.is_finite()
        {
            return Err(Error("invalid spins or Hamiltonian"));
        }
        let mut out = Self {
            lattice,
            spins,
            j,
            h,
            energy: 0.0,
        };
        out.energy = out.recompute_energy();
        if !out.energy.is_finite() {
            return Err(Error("energy overflow"));
        }
        Ok(out)
    }
    pub fn aligned(size: usize) -> Result<Self> {
        let lattice = Lattice2D::new(size, size, Boundary::Periodic).map_err(Error)?;
        let mut spins = Vec::new();
        spins
            .try_reserve_exact(lattice.len())
            .map_err(|_| Error("lattice allocation failed"))?;
        spins.resize(lattice.len(), 1);
        Self::new(size, size, Boundary::Periodic, 1.0, 0.0, spins)
    }
    pub fn spins(&self) -> &[i8] {
        &self.spins
    }
    pub fn magnetization(&self) -> f64 {
        self.spins.iter().map(|s| *s as f64).sum()
    }
    pub fn recompute_energy(&self) -> f64 {
        -self.j
            * self
                .lattice
                .bonds()
                .map(|(i, j)| self.spins[i] as f64 * self.spins[j] as f64)
                .sum::<f64>()
            - self.h * self.magnetization()
    }
    pub fn flip_energy(&self, i: usize) -> Result<f64> {
        let neighbors = self.lattice.neighbors(i).map_err(Error)?;
        Ok(2.0
            * self.spins[i] as f64
            * (self.j * neighbors.iter().map(|j| self.spins[*j] as f64).sum::<f64>() + self.h))
    }
    /// One Wolff cluster, only finite beta >= 0, J >= 0, h == 0.
    /// Returns flipped spins; this is a cluster update, not a single-spin sweep.
    pub fn wolff_step(&mut self, beta: f64, rng: &mut impl RngExt) -> Result<usize> {
        if !beta.is_finite() || beta < 0.0 || self.j < 0.0 || self.h != 0.0 {
            return Err(Error("Wolff requires beta >= 0, J >= 0 and zero field"));
        }
        let p = if beta == 0.0 || self.j == 0.0 {
            0.0
        } else {
            -(-2.0 * beta * self.j).exp_m1()
        };
        let root = rng.random_range(0..self.spins.len());
        let spin = self.spins[root];
        let mut marked = vec![false; self.spins.len()];
        marked[root] = true;
        let mut cluster = vec![root];
        let mut pos = 0;
        while pos < cluster.len() {
            let site = cluster[pos];
            pos += 1;
            for j in self.lattice.neighbors(site).map_err(Error)? {
                if !marked[j] && self.spins[j] == spin && rng.random::<f64>() < p {
                    marked[j] = true;
                    cluster.push(j);
                }
            }
        }
        for i in &cluster {
            self.spins[*i] *= -1;
        }
        let next = self.recompute_energy();
        if !next.is_finite() {
            for i in &cluster {
                self.spins[*i] *= -1;
            }
            return Err(Error("cluster energy overflow"));
        }
        self.energy = next;
        Ok(cluster.len())
    }
}
impl State for Ising2D {
    type Params = ();
    type Change = (usize, f64, f64);
    fn valid_state(&self, _: &()) -> bool {
        self.lattice.valid()
            && self.spins.len() == self.lattice.len()
            && self.spins.iter().all(|s| *s == 1 || *s == -1)
            && self.j.is_finite()
            && self.h.is_finite()
            && self.energy.is_finite()
            && (self.energy - self.recompute_energy()).abs() <= 1e-9 * (1.0 + self.energy.abs())
    }
    fn energy(&self, _: &mut ()) -> f64 {
        self.energy
    }
    fn propose_change(&self, rng: &mut impl RngExt) -> Self::Change {
        let i = rng.random_range(0..self.spins.len());
        (i, self.energy, self.flip_energy(i).expect("valid site"))
    }
    fn apply_change(&mut self, (i, old, d): Self::Change) {
        self.spins[i] *= -1;
        self.energy = old + d;
    }
    fn revert_change(&mut self, (i, old, _): Self::Change) {
        self.spins[i] *= -1;
        self.energy = old;
    }
    fn delta_energy(&self, c: &Self::Change, _: &()) -> Option<f64> {
        Some(c.2)
    }
}

/// Lattice gas H=-J sum n_i n_j, sampled using effective energy H-mu*N.
/// A uniformly chosen site is toggled, so insertion/deletion proposals are symmetric.
#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug)]
pub struct LatticeGas {
    lattice: Lattice2D,
    occupation: Vec<bool>,
    j: f64,
    mu: f64,
    energy: f64,
}
impl LatticeGas {
    pub fn new(
        width: usize,
        height: usize,
        boundary: Boundary,
        j: f64,
        mu: f64,
        occupation: Vec<bool>,
    ) -> Result<Self> {
        let lattice = Lattice2D::new(width, height, boundary).map_err(Error)?;
        if occupation.len() != lattice.len() || !j.is_finite() || !mu.is_finite() {
            return Err(Error("invalid lattice gas"));
        }
        let mut out = Self {
            lattice,
            occupation,
            j,
            mu,
            energy: 0.0,
        };
        out.energy = out.recompute_energy();
        if !out.energy.is_finite() {
            return Err(Error("energy overflow"));
        }
        Ok(out)
    }
    pub fn particles(&self) -> usize {
        self.occupation.iter().filter(|n| **n).count()
    }
    pub fn occupation(&self) -> &[bool] {
        &self.occupation
    }
    pub fn recompute_energy(&self) -> f64 {
        -self.j
            * self
                .lattice
                .bonds()
                .filter(|(i, j)| self.occupation[*i] && self.occupation[*j])
                .count() as f64
            - self.mu * self.particles() as f64
    }
    pub fn physical_energy(&self) -> f64 {
        -self.j
            * self
                .lattice
                .bonds()
                .filter(|(i, j)| self.occupation[*i] && self.occupation[*j])
                .count() as f64
    }
}
impl State for LatticeGas {
    type Params = ();
    type Change = (usize, f64, f64);
    fn valid_state(&self, _: &()) -> bool {
        self.lattice.valid()
            && self.occupation.len() == self.lattice.len()
            && self.j.is_finite()
            && self.mu.is_finite()
            && self.energy.is_finite()
            && (self.energy - self.recompute_energy()).abs() <= 1e-9 * (1.0 + self.energy.abs())
    }
    fn energy(&self, _: &mut ()) -> f64 {
        self.energy
    }
    fn propose_change(&self, r: &mut impl RngExt) -> Self::Change {
        let i = r.random_range(0..self.occupation.len());
        let dn = if self.occupation[i] { -1.0 } else { 1.0 };
        let neighbors = self
            .lattice
            .neighbors(i)
            .expect("valid site")
            .iter()
            .filter(|j| self.occupation[**j])
            .count();
        (i, self.energy, -dn * (self.j * neighbors as f64 + self.mu))
    }
    fn apply_change(&mut self, (i, old, d): Self::Change) {
        self.occupation[i] = !self.occupation[i];
        self.energy = old + d;
    }
    fn revert_change(&mut self, (i, old, _): Self::Change) {
        self.occupation[i] = !self.occupation[i];
        self.energy = old;
    }
    fn delta_energy(&self, c: &Self::Change, _: &()) -> Option<f64> {
        Some(c.2)
    }
}
