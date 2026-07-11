use clap::Parser;

use std::{num::NonZeroUsize, path::PathBuf};

use crate::{
    align::Mode,
    io::{OutputType, bed4::DEF_NAME_COL},
};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Cli {
    /// First BED file.
    #[arg(short = 't', long)]
    pub infile_target: PathBuf,
    /// Second BED file.
    #[arg(short = 'q', long)]
    pub infile_query: PathBuf,
    /// 1-based index for labels
    #[arg(short = 'c', long, default_value_t = DEF_NAME_COL)]
    pub label_col: NonZeroUsize,
    /// Output file. If omitted, defaults to stdout.
    #[arg(short = 'o', long)]
    pub outfile: Option<PathBuf>,
    /// Output file type.
    /// Either BEDPE or PAF.
    #[arg(short = 'y', long, default_value = "BEDPE")]
    pub outfile_type: OutputType,
    /// Match score
    #[arg(short = 'm', long, default_value_t = 2.0)]
    pub score_match: f32,
    /// Mismatch score
    #[arg(short = 'x', long, default_value_t = -1.0)]
    pub score_mismatch: f32,
    /// Gap-open penalty
    #[arg(short = 'p', long, default_value_t = -4.0)]
    pub score_gap_open: f32,
    /// Gap-extension penalty
    #[arg(short = 'e', long, default_value_t = -1.0)]
    pub score_gap_ext: f32,
    /// Alignment mode. Either global or local.
    /// Local only returns the highest scoring alignment.
    #[arg(short = 'a', long, default_value = "global")]
    pub mode: Mode,
    /// Minimum alignment score. Only valid with local alignment.
    #[arg(short = 's', long, default_value_t = 5.0)]
    pub minimum_aln_score: f32,
}
