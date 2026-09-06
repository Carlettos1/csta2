//! Versioned JSON snapshots. Write a sibling temporary file, sync, then rename.
//! Model tags are caller-controlled schema/Hamiltonian versions. Load only local
//! trusted snapshots; structural model validation belongs to the model's loader.
use serde::{Deserialize, Serialize, de::DeserializeOwned};
use std::{
    fs::OpenOptions,
    io::{Read, Write},
    path::Path,
};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;
#[derive(Serialize, Deserialize)]
struct Envelope<T> {
    format: u32,
    implementation: String,
    model: String,
    rust_type: String,
    value: T,
}
pub fn save<T: Serialize>(path: impl AsRef<Path>, model: &str, value: &T) -> Result<()> {
    let path = path.as_ref();
    if model.is_empty() {
        return Err("checkpoint needs a model version tag".into());
    }
    let bytes = serde_json::to_vec(&Envelope {
        format: 1,
        implementation: "csta-3/checkpoint-1/rand-0.10.2".into(),
        model: model.into(),
        rust_type: std::any::type_name::<T>().into(),
        value,
    })?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut temporary = None;
    for i in 0..1000 {
        let candidate = parent.join(format!(".csta-checkpoint-{}-{i}.tmp", std::process::id()));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(file) => {
                temporary = Some((candidate, file));
                break;
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e.into()),
        }
    }
    let (temporary, mut file) = temporary.ok_or("cannot create checkpoint temporary file")?;
    let result = (|| -> Result<()> {
        file.write_all(&bytes)?;
        file.sync_all()?;
        drop(file);
        std::fs::rename(&temporary, path)?;
        std::fs::File::open(parent)?.sync_all()?;
        Ok(())
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(temporary);
    }
    result
}
pub fn load<T: DeserializeOwned>(path: impl AsRef<Path>, model: &str) -> Result<T> {
    // Bounded read rejects accidental huge files before JSON allocation.
    const LIMIT: u64 = 256 * 1024 * 1024;
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take(LIMIT + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 > LIMIT {
        return Err("checkpoint exceeds 256 MiB limit".into());
    }
    let e: Envelope<T> = serde_json::from_slice(&bytes)?;
    if e.format != 1
        || e.implementation != "csta-3/checkpoint-1/rand-0.10.2"
        || e.model != model
        || e.rust_type != std::any::type_name::<T>()
    {
        return Err("incompatible checkpoint schema, model or RNG type".into());
    }
    Ok(e.value)
}
/// Preserve finite values and -infinity exactly; JSON otherwise encodes infinity
/// as null. DOS bit patterns are validated by the sampler after decoding.
pub mod float_bits {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    pub fn serialize<S: Serializer>(values: &[f64], s: S) -> Result<S::Ok, S::Error> {
        values
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<_>>()
            .serialize(s)
    }
    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<Vec<f64>, D::Error> {
        Ok(Vec::<u64>::deserialize(d)?
            .into_iter()
            .map(f64::from_bits)
            .collect())
    }
}
