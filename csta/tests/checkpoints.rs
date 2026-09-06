#![cfg(feature = "checkpoint")]
use csta::{CheckpointRng, Metropolis, State, models::Ising2D, wl};
use rand::{RngExt, SeedableRng};
fn path(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("csta-{}-{name}.json", std::process::id()))
}
#[test]
fn canonical_exact_resume_and_validation() {
    let file = path("canonical");
    let model = Ising2D::aligned(3).unwrap();
    let mut original = Metropolis::with_all(model, (), 0.4, 0, CheckpointRng::seed_from_u64(20));
    for _ in 0..35 {
        original.try_step().unwrap();
    }
    original.save(&file, "ising-v1").unwrap();
    let mut resumed = Metropolis::<Ising2D, CheckpointRng>::load(&file, "ising-v1").unwrap();
    assert!(Metropolis::<Ising2D, CheckpointRng>::load(&file, "ising-v2").is_err());
    for _ in 0..100 {
        assert_eq!(original.try_step().unwrap(), resumed.try_step().unwrap());
        assert_eq!(original.state.spins(), resumed.state.spins());
    }
    assert_eq!(original.attempted_moves, resumed.attempted_moves);
    assert_eq!(original.rng.random::<u64>(), resumed.rng.random::<u64>());
    // Failed save validates before replacing a good snapshot.
    original.beta = f64::NAN;
    assert!(original.save(&file, "ising-v1").is_err());
    assert!(Metropolis::<Ising2D, CheckpointRng>::load(&file, "ising-v1").is_ok());
    std::fs::write(&file, b"{broken").unwrap();
    assert!(Metropolis::<Ising2D, CheckpointRng>::load(&file, "ising-v1").is_err());
    std::fs::remove_file(file).unwrap();
}
fn config() -> wl::Config {
    wl::Config {
        preliminary_stages: 2,
        min_stage_steps: 5,
        visits_per_bin: 2,
        sampling_steps: 1000,
        max_steps: 20_000,
        ..Default::default()
    }
}
// A 3x3 periodic square lattice has integer energies; enumerate its exact support.
fn grid() -> wl::EnergyGrid {
    let mut energies = Vec::new();
    for bits in 0..512 {
        let spins = (0..9)
            .map(|i| if bits & (1 << i) == 0 { -1 } else { 1 })
            .collect();
        let model = Ising2D::new(3, 3, csta::models::Boundary::Periodic, 1.0, 0.0, spins).unwrap();
        energies.push(model.energy(&mut ()));
    }
    energies.sort_by(f64::total_cmp);
    energies.dedup();
    wl::EnergyGrid::discrete(energies).unwrap()
}
#[test]
fn serial_wl_exact_resume_during_warmup_and_production() {
    for steps in [1, 200, 500] {
        let file = path(&format!("serial-{steps}"));
        let mut a = wl::Session::new(
            Ising2D::aligned(3).unwrap(),
            (),
            CheckpointRng::seed_from_u64(5),
            wl::RawWangLandauData::on_grid(grid()).unwrap(),
            config(),
        )
        .unwrap();
        a.advance(steps, None).unwrap();
        a.save(&file, "ising-v1").unwrap();
        let mut b = wl::Session::<Ising2D, CheckpointRng>::load(&file, "ising-v1").unwrap();
        a.advance(20_000, None).unwrap();
        b.advance(20_000, None).unwrap();
        assert_eq!(a.diagnostics(), b.diagnostics());
        assert_eq!(a.data().dos(), b.data().dos());
        assert_eq!(a.data().lifetime_bins(), b.data().lifetime_bins());
        assert_eq!(a.state().spins(), b.state().spins());
        assert!(a.finish().is_complete());
        assert!(b.finish().is_complete());
        std::fs::remove_file(file).unwrap();
    }
}
#[test]
fn parallel_wl_exact_barrier_resume() {
    let file = path("parallel");
    let g = grid();
    let n = g.len();
    let config = wl::ParallelConfig {
        windows: std::iter::once(0..n).collect(),
        walkers_per_window: 2,
        exchange_every: 13,
        seed: 2,
        run: config(),
    };
    let mut a = wl::ParallelSession::<_, CheckpointRng>::new(
        g,
        vec![Ising2D::aligned(3).unwrap(); 2],
        (),
        config,
    )
    .unwrap();
    a.advance(4, None).unwrap();
    a.save(&file, "ising-v1").unwrap();
    let mut b = wl::ParallelSession::<Ising2D, CheckpointRng>::load(&file, "ising-v1").unwrap();
    a.advance(10_000, None).unwrap();
    b.advance(10_000, None).unwrap();
    let a = a.finish().unwrap();
    let b = b.finish().unwrap();
    assert_eq!(a.stop_reason, wl::StopReason::TargetReached);
    assert_eq!(a.exchange_attempts, b.exchange_attempts);
    assert_eq!(a.merged.unwrap().dos(), b.merged.unwrap().dos());
    for (a, b) in a.walkers.iter().zip(b.walkers) {
        assert_eq!(a.diagnostics, b.diagnostics);
        assert_eq!(a.state.spins(), b.state.spins());
    }
    std::fs::remove_file(file).unwrap();
}

#[test]
fn snapshot_corruption_and_cancel_resume() {
    let file = path("corrupt");
    let mut session = wl::Session::new(
        wl::models::Ising::<4>::new([1; 4]).unwrap(),
        (),
        CheckpointRng::seed_from_u64(42),
        wl::RawWangLandauData::on_grid(wl::models::Ising::<4>::grid().unwrap()).unwrap(),
        config(),
    )
    .unwrap();
    session.advance(30, None).unwrap();
    session.save(&file, "chain-v1").unwrap();
    let valid: serde_json::Value = serde_json::from_slice(&std::fs::read(&file).unwrap()).unwrap();
    for pointer in [
        "/format",
        "/value/walker/data/visits",
        "/value/walker/bin",
        "/value/walker/diagnostics/accepted",
        "/value/walker/state/energy",
    ] {
        let mut bad = valid.clone();
        *bad.pointer_mut(pointer).unwrap() = serde_json::json!(100_000);
        std::fs::write(&file, serde_json::to_vec(&bad).unwrap()).unwrap();
        assert!(
            wl::Session::<wl::models::Ising<4>, CheckpointRng>::load(&file, "chain-v1").is_err(),
            "{pointer}"
        );
    }
    let flag = std::sync::atomic::AtomicBool::new(true);
    session.advance(1, Some(&flag)).unwrap();
    assert_eq!(
        session.diagnostics().stop_reason,
        Some(wl::StopReason::Cancelled)
    );
    session.save(&file, "chain-v1").unwrap();
    let mut resumed =
        wl::Session::<wl::models::Ising<4>, CheckpointRng>::load(&file, "chain-v1").unwrap();
    session.advance(20_000, None).unwrap();
    resumed.advance(20_000, None).unwrap();
    assert_eq!(session.data().dos(), resumed.data().dos());
    resumed.save(&file, "chain-v1").unwrap();
    let mut done =
        wl::Session::<wl::models::Ising<4>, CheckpointRng>::load(&file, "chain-v1").unwrap();
    let before = done.data().dos().to_vec();
    done.advance(10, None).unwrap();
    assert_eq!(before, done.data().dos());
    std::fs::remove_file(file).unwrap();
}
#[test]
fn masked_dos_parallel_checkpoint() {
    let file = path("masked");
    let g = wl::models::Ising::<6>::grid().unwrap();
    let cold = wl::models::Ising::<6>::new([1; 6]).unwrap();
    let hot = wl::models::Ising::<6>::new([1, -1, 1, -1, 1, -1]).unwrap();
    let config = wl::ParallelConfig {
        windows: vec![0..3, 1..4],
        walkers_per_window: 1,
        exchange_every: 7,
        seed: 42,
        run: wl::Config {
            preliminary_stages: 0,
            sampling_steps: 1000,
            max_steps: 2000,
            ..Default::default()
        },
    };
    let mut a =
        wl::ParallelSession::<_, CheckpointRng>::new(g, vec![cold, hot], (), config).unwrap();
    a.advance(4, None).unwrap();
    a.save(&file, "masked-v1").unwrap();
    let mut b =
        wl::ParallelSession::<wl::models::Ising<6>, CheckpointRng>::load(&file, "masked-v1")
            .unwrap();
    a.advance(2000, None).unwrap();
    b.advance(2000, None).unwrap();
    let a = a.finish().unwrap();
    let b = b.finish().unwrap();
    assert_eq!(a.exchange_attempts, b.exchange_attempts);
    assert_eq!(a.accepted_exchanges, b.accepted_exchanges);
    assert_eq!(a.merged.unwrap().dos(), b.merged.unwrap().dos());
    std::fs::remove_file(file).unwrap();
}
#[test]
fn temperature_exchange_checkpoint_and_cancelled_chunk() {
    use csta::replica::ReplicaExchange;
    let file = path("temperature");
    let states = vec![Ising2D::aligned(3).unwrap(); 3];
    let mut a =
        ReplicaExchange::<_, CheckpointRng>::seeded(states, (), vec![0.1, 0.3, 0.5], 42).unwrap();
    a.run(5, 11, None).unwrap();
    let flag = std::sync::atomic::AtomicBool::new(true);
    a.run(1, 11, Some(&flag)).unwrap();
    a.save(&file, "temperature-v1").unwrap();
    let mut b = ReplicaExchange::<Ising2D, CheckpointRng>::load(&file, "temperature-v1").unwrap();
    assert!(b.run(1, 12, None).is_err());
    a.run(40, 11, None).unwrap();
    b.run(40, 11, None).unwrap();
    assert_eq!(a.swap_accepts, b.swap_accepts);
    assert_eq!(a.round_trips, b.round_trips);
    for (a, b) in a.replicas().iter().zip(b.replicas()) {
        assert_eq!(a.state.spins(), b.state.spins());
        assert_eq!(a.attempted_moves, b.attempted_moves);
    }
    std::fs::remove_file(file).unwrap();
}
