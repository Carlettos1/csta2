//! Sparse joint DOS coordinates (field-free energy E0, total magnetization M).
use crate::*;

#[derive(Clone, Debug)]
pub struct JointGrid {
    cells: Vec<(f64, f64)>,
}
impl JointGrid {
    /// Cells must be lexicographically sorted, unique, and finite. Sparse cells
    /// are an explicit support declaration, never inferred from absent visits.
    pub fn new(cells: Vec<(f64, f64)>, max_cells: usize) -> Result<Self> {
        if cells.is_empty()
            || cells.len() > max_cells
            || cells.iter().any(|(e, m)| !e.is_finite() || !m.is_finite())
            || cells.windows(2).any(|c| c[0] >= c[1])
        {
            return Err(Error::Invalid("invalid joint support or memory budget"));
        }
        Ok(Self { cells })
    }
    pub fn cells(&self) -> &[(f64, f64)] {
        &self.cells
    }
    fn index(&self, cell: (f64, f64)) -> Option<usize> {
        self.cells
            .binary_search_by(|x| x.partial_cmp(&cell).unwrap())
            .ok()
    }
}
#[derive(Clone, Debug)]
pub struct JointDos {
    grid: JointGrid,
    log_dos: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct JointThermodynamics {
    pub log_partition: f64,
    pub energy: f64,
    pub magnetization: f64,
    pub probabilities: Vec<f64>,
}
impl JointDos {
    pub fn new(grid: JointGrid, log_dos: Vec<f64>) -> Result<Self> {
        if log_dos.len() != grid.cells.len() || log_dos.iter().any(|x| !x.is_finite()) {
            return Err(Error::Invalid(
                "joint DOS needs finite mass on declared support",
            ));
        }
        Ok(Self { grid, log_dos })
    }
    pub fn normalize_log_count(&mut self, log_count: f64) -> Result<()> {
        if !log_count.is_finite() {
            return Err(Error::Invalid("invalid state count"));
        }
        let z = log_sum_exp(&self.log_dos)?;
        let next: Vec<_> = self.log_dos.iter().map(|x| x - z + log_count).collect();
        if next.iter().any(|x| !x.is_finite()) {
            return Err(Error::Numerical("joint normalization overflow"));
        }
        self.log_dos = next;
        Ok(())
    }
    pub fn evaluate(&self, beta: f64, field: f64) -> Result<JointThermodynamics> {
        if !beta.is_finite() || !field.is_finite() {
            return Err(Error::Invalid("invalid beta or field"));
        }
        let energies: Vec<_> = self.grid.cells.iter().map(|(e, m)| e - field * m).collect();
        if energies.iter().any(|x| !x.is_finite()) {
            return Err(Error::Numerical("joint energy overflow"));
        }
        let reference = energies
            .iter()
            .copied()
            .reduce(|a, b| if beta >= 0.0 { a.min(b) } else { a.max(b) })
            .unwrap();
        let offset = self
            .log_dos
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let logs: Vec<_> = energies
            .iter()
            .zip(&self.log_dos)
            .map(|(e, g)| {
                csta_metropolis::boltzmann_log_ratio(beta, reference, *e)
                    .map(|r| (g - offset) + r)
                    .map_err(|_| Error::Numerical("joint weights invalid"))
            })
            .collect::<Result<_>>()?;
        let shifted_z = log_sum_exp(&logs)?;
        if !shifted_z.is_finite() {
            return Err(Error::Numerical("joint weights unrepresentable"));
        }
        let z = (offset - beta * reference) + shifted_z;
        if !z.is_finite() {
            return Err(Error::Numerical("joint partition overflow"));
        }
        let probabilities: Vec<_> = logs.iter().map(|x| (x - shifted_z).exp()).collect();
        let energy = probabilities
            .iter()
            .zip(&energies)
            .map(|(p, e)| p * e)
            .sum::<f64>();
        let magnetization = probabilities
            .iter()
            .zip(&self.grid.cells)
            .map(|(p, (_, m))| p * m)
            .sum::<f64>();
        if !energy.is_finite() || !magnetization.is_finite() {
            return Err(Error::Numerical("joint moment overflow"));
        }
        Ok(JointThermodynamics {
            log_partition: z,
            energy,
            magnetization,
            probabilities,
        })
    }
    pub fn energy_marginal(&self) -> Result<WLData> {
        let mut energies = Vec::new();
        let mut logs: Vec<f64> = Vec::new();
        for ((e, _), g) in self.grid.cells.iter().zip(&self.log_dos) {
            if energies.last() == Some(e) {
                let last = logs.last_mut().unwrap();
                *last = log_sum_exp(&[*last, *g])?;
            } else {
                energies.push(*e);
                logs.push(*g);
            }
        }
        let n = energies.len();
        WLData::new(EnergyGrid::discrete(energies)?, logs, vec![0; n])
    }
}
struct JointState<S, F> {
    state: S,
    grid: JointGrid,
    magnetization: F,
}
impl<S: State, F: Fn(&S) -> f64> State for JointState<S, F> {
    type Params = S::Params;
    type Change = S::Change;
    fn energy(&self, p: &mut Self::Params) -> f64 {
        let e = self.state.energy(p);
        let m = (self.magnetization)(&self.state);
        if !e.is_finite() || !m.is_finite() {
            return f64::NAN;
        }
        self.grid.index((e, m)).map_or(-1.0, |i| i as f64)
    }
    fn propose_change(&self, r: &mut impl RngExt) -> Self::Change {
        self.state.propose_change(r)
    }
    fn apply_change(&mut self, c: Self::Change) {
        self.state.apply_change(c)
    }
    fn revert_change(&mut self, c: Self::Change) {
        self.state.revert_change(c)
    }
    fn log_proposal_ratio(&self, c: &Self::Change) -> f64 {
        self.state.log_proposal_ratio(c)
    }
}
pub struct JointRun<S: State> {
    pub state: S,
    pub params: S::Params,
    pub data: JointDos,
    pub diagnostics: Diagnostics,
}
impl<S: State> JointRun<S> {
    pub fn is_complete(&self) -> bool {
        self.diagnostics.stop_reason == Some(StopReason::TargetReached)
    }
}
/// Uses existing WL/SAMC budgets and cancellation. State::energy must be E0,
/// without an external magnetic field term. Magnetization closure is read-only.
pub fn run_joint<S: State, R: RngExt>(
    state: S,
    params: S::Params,
    rng: R,
    grid: JointGrid,
    magnetization: impl Fn(&S) -> f64,
    config: Config,
    cancel: Option<&AtomicBool>,
) -> Result<JointRun<S>>
where
    S::Params: Clone,
{
    let index = EnergyGrid::discrete((0..grid.cells.len()).map(|i| i as f64).collect())?;
    let wrapped = JointState {
        state,
        grid: grid.clone(),
        magnetization,
    };
    let result = run(
        wrapped,
        params,
        rng,
        RawWangLandauData::on_grid(index)?,
        config,
        cancel,
    )?;
    Ok(JointRun {
        state: result.state.state,
        params: result.params,
        data: JointDos::new(grid, result.data.dos)?,
        diagnostics: result.diagnostics,
    })
}
