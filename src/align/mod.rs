use std::str::FromStr;

use eyre::bail;

pub mod global;
pub mod local;

#[derive(Debug, Clone, Copy)]
pub enum Mode {
    /// Global alignment via Needleman-Wunsch with affine gap penalties
    Global,
    /// Local alignment via Smith-Waterman with affine gap penalties
    Local,
}

impl FromStr for Mode {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "local" => Ok(Mode::Local),
            "global" => Ok(Mode::Global),
            _ => bail!("Invalid mode: {s}"),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum TraceOp {
    /// Match or mismatch at (i, j)
    M,
    /// Gap in query (insertion in target) at (i, j)
    X,
    /// Gap in target (deletion in target) at (i, j)
    Y,
}

#[derive(Debug, Clone, Copy)]
pub struct Trace {
    /// Trace operation
    pub op: TraceOp,
    /// Index in target
    pub i: usize,
    /// Index in query
    pub j: usize,
}

impl Trace {
    pub fn new(op: TraceOp, i: usize, j: usize) -> Self {
        Trace { op, i, j }
    }
}
