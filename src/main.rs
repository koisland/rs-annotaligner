mod align;
mod cli;
mod io;

use crate::{
    align::{global::needleman_wuncsh_affine, local::smith_waterman_affine},
    cli::Commands,
    io::{align::write_alignments, bed4::read_bed4},
};
use clap::Parser;

fn main() -> eyre::Result<()> {
    let args = cli::Cli::parse();

    match args.command {
        Commands::Global { args } => {
            let bed_target = read_bed4(&args.infile_target, Some(args.label_col))?;
            let bed_query = read_bed4(&args.infile_query, Some(args.label_col))?;
            let aln = needleman_wuncsh_affine(
                &bed_target,
                &bed_query,
                args.score_match,
                args.score_mismatch,
                args.score_gap_open,
                args.score_gap_ext,
            )?;
            // Only generates best for now?
            write_alignments(vec![aln], args.outfile_type, args.outfile.as_deref())?;
        }
        Commands::Local {
            args,
            minimum_aln_score,
        } => {
            let bed_target = read_bed4(&args.infile_target, Some(args.label_col))?;
            let bed_query = read_bed4(&args.infile_query, Some(args.label_col))?;

            let alns = smith_waterman_affine(
                &bed_target,
                &bed_query,
                args.score_match,
                args.score_mismatch,
                args.score_gap_open,
                args.score_gap_ext,
                minimum_aln_score,
            )?;
            write_alignments(alns, args.outfile_type, args.outfile.as_deref())?;
        }
    };

    Ok(())
}
