use csta::{State, wl};

// This fixture deliberately implements csta::State rather than wl::State. The
// test catches accidental duplicate registry/path trait identities at integration.
struct Fixed;
impl State for Fixed {
    type Params = ();
    type Change = ();
    fn energy(&self, _: &mut ()) -> f64 {
        0.0
    }
    fn propose_change(&self, _: &mut impl rand::Rng) {}
    fn apply_change(&mut self, _: ()) {}
    fn revert_change(&mut self, _: ()) {}
}
#[test]
fn umbrella_state_runs_through_wl() {
    let config = wl::Config {
        preliminary_stages: 0,
        sampling_steps: 3,
        ..wl::Config::default()
    };
    let result = wl::run(
        Fixed,
        (),
        rand::rng(),
        wl::RawWangLandauData::on_grid(wl::EnergyGrid::discrete(vec![0.0]).unwrap()).unwrap(),
        config,
        None,
    )
    .unwrap();
    assert!(result.is_complete());
    assert_eq!(result.data.total_visits(), 3);
}
