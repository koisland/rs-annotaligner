use std::str::FromStr;

use eyre::bail;

pub mod global;
pub mod local;

#[derive(Debug, Clone, Copy)]
pub enum TraceOp {
    /// Match or mismatch at (i, j)
    M,
    /// Insertion in target (deletion in query) at (i, j)
    X,
    /// Gap in target (insertion in query) at (i, j)
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CigarOp {
    Match,
    Mismatch,
    Insertion,
    Deletion,
}

impl From<CigarOp> for char {
    fn from(op: CigarOp) -> Self {
        match op {
            CigarOp::Match => '=',
            CigarOp::Mismatch => 'X',
            CigarOp::Insertion => 'I',
            CigarOp::Deletion => 'D',
        }
    }
}

impl FromStr for CigarOp {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "match" | "Match" => CigarOp::Match,
            "mismatch" | "Mismatch" => CigarOp::Mismatch,
            "insertion" | "Insertion" => CigarOp::Insertion,
            "deletion" | "Deletion" => CigarOp::Deletion,
            _ => bail!("Invalid operation ({s})"),
        })
    }
}

impl Trace {
    pub fn new(op: TraceOp, i: usize, j: usize) -> Self {
        Trace { op, i, j }
    }
}
