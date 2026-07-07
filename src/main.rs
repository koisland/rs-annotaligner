use clap::Parser;

mod align;
mod cli;
mod io;

use crate::{
    align::needleman_wuncsh_affine,
    io::{read_bed, write_bedpe},
};

fn main() -> eyre::Result<()> {
    let args = cli::Cli::parse();

    let bed_target = read_bed(&args.infile_target, Some(args.label_col))?;
    let bed_query = read_bed(&args.infile_query, Some(args.label_col))?;

    let alns = needleman_wuncsh_affine(
        &bed_target,
        &bed_query,
        args.score_match,
        args.score_mismatch,
        args.score_gap_open,
        args.score_gap_ext,
    );
    write_bedpe(&alns, args.outfile)?;

    Ok(())
}
