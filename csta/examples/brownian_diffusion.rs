//! Free diffusion in three dimensions: ensemble MSD = 6 D t.
//! Run `cargo run -p csta --release --example brownian_diffusion`.
use csta::{Vec3f64, gaussian, statistics::Moments};
use rand::{SeedableRng, rngs::StdRng};

const DIFFUSION: f64 = 0.7;
const DT: f64 = 0.02;
const PATHS: usize = 4096;
const TIMES: [usize; 5] = [25, 50, 100, 200, 400];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut msd = [Moments::default(); TIMES.len()];
    let scale = (2.0 * DIFFUSION * DT).sqrt();
    for path in 0..PATHS {
        let mut rng = StdRng::seed_from_u64(1000 + path as u64);
        let mut position = Vec3f64::default();
        let mut observation = 0;
        for step in 1..=TIMES[TIMES.len() - 1] {
            position += Vec3f64::new(
                gaussian(&mut rng, 0.0, scale)?,
                gaussian(&mut rng, 0.0, scale)?,
                gaussian(&mut rng, 0.0, scale)?,
            );
            if step == TIMES[observation] {
                msd[observation].push(position.len_squared())?;
                observation += 1;
            }
        }
    }
    println!("time,msd,msd_se,msd_exact,diffusion_estimate,diffusion_se,diffusion_exact");
    for (steps, moments) in TIMES.into_iter().zip(msd) {
        let time = steps as f64 * DT;
        let mean = moments.mean().ok_or("missing trajectories")?;
        let se = (moments.variance().ok_or("need multiple trajectories")? / PATHS as f64).sqrt();
        println!(
            "{time:.4},{mean:.8},{se:.8},{:.8},{:.8},{:.8},{DIFFUSION}",
            6.0 * DIFFUSION * time,
            mean / (6.0 * time),
            se / (6.0 * time)
        );
    }
    eprintln!(
        "Each SE uses {PATHS} separate trajectories at fixed time. Rows share trajectories and are correlated; they are not independent measurements for a naive line-fit error. Positions are unwrapped."
    );
    Ok(())
}
