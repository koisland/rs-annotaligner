use std::str::FromStr;

use eyre::bail;

pub mod align;
pub mod bed4;
pub mod bedpe;
pub mod paf;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputType {
    PAF,
    BEDPE,
}

impl FromStr for OutputType {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_ascii_lowercase().as_str() {
            "paf" => OutputType::PAF,
            "bedpe" => OutputType::BEDPE,
            _ => bail!("Invalid output type ({s}"),
        })
    }
}
