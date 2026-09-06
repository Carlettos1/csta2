#[cfg(feature = "checkpoint")]
fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    use csta::{CheckpointRng, Metropolis, models::Ising2D};
    use rand::SeedableRng;
    let path = std::env::temp_dir().join(format!("csta-example-{}.json", std::process::id()));
    let mut run = Metropolis::with_all(
        Ising2D::aligned(4)?,
        (),
        0.4,
        100,
        CheckpointRng::seed_from_u64(42),
    );
    run.run(None)?;
    run.save(&path, "ising-j1-h0-v1")?;
    let mut resumed = Metropolis::<Ising2D, CheckpointRng>::load(&path, "ising-j1-h0-v1")?;
    resumed.run(None)?;
    println!("Resumed cumulative attempts: {}", resumed.attempted_moves);
    std::fs::remove_file(path)?;
    Ok(())
}
#[cfg(not(feature = "checkpoint"))]
fn main() {
    eprintln!("Run with --features checkpoint to enable persistence.");
}
