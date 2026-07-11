use eyre::bail;
use itertools::Itertools;

use crate::io::align::Alignment;

/// A [`PAF`](https://github.com/lh3/miniasm/blob/master/PAF.md) record for aligned [`BEDPE`] records.
///
/// # Reference
/// * https://lh3.github.io/minimap2/minimap2.html
/// * https://github.com/lh3/miniasm/blob/master/PAF.md
#[derive(Debug, Clone)]
pub struct PAF(Vec<Alignment>);

impl PAF {
    /// Create new PAF and calculate MAPQ.
    ///
    /// Calculate as follows:
    /// `30 * log_10(100.0 * 0.5 * (identity + (dp_score / top_dp_score)))`
    /// * Equal weight given to DP score and gap-compressed identity.
    /// * Limited to 0 to 60.
    pub(crate) fn new(mut alignments: Vec<Alignment>) -> eyre::Result<Self> {
        let best_aln_score = {
            let Some(best_aln) = alignments.iter().max() else {
                bail!("No possible alignments")
            };
            best_aln.score
        };

        // MAPQ from bwa-mem based on base quality and the best and second-best hit
        // * https://genome.cshlp.org/content/18/11/1851
        // Minimap2's is based on the first and second best chain and the number of anchors
        // * https://academic.oup.com/bioinformatics/article/34/18/3094/4994778?login=false
        // This score should penalize short alignments via the max DP score and sequence identity.
        // * We cap the score at 60.
        for aln in &mut alignments {
            let score = 100.0 * *((aln.gap_cmp_divergence + (aln.score / best_aln_score)) / 2.0);
            let mapq = (30.0 * score.log10()) as u8;
            aln.mapq = Some(mapq)
        }
        Ok(PAF(alignments))
    }

    /// Get rows as iterator.
    pub fn rows(&self) -> impl Iterator<Item = &Alignment> {
        self.0.iter()
    }

    /// Generate alignment as PAF row strings sorted by MAPQ.
    pub fn as_str(&self) -> impl Iterator<Item = String> {
        self.rows()
            .sorted_by(|a, b| a.mapq.cmp(&b.mapq))
            .rev()
            .map(|aln| aln.as_paf_row())
    }
}

#[cfg(test)]
mod test {
    use crate::{
        align::local::smith_waterman_affine,
        io::{bed4::read_bed4, paf::PAF},
    };
    use std::{
        fs::File,
        io::{BufRead, BufReader},
        path::PathBuf,
    };

    const INPUT_DIR: &str = "test/data/input";
    const EXP_DIR: &str = "test/data/output";

    #[test]
    fn test_local_paf() {
        let t = PathBuf::from(INPUT_DIR).join("target_local.bed");
        let q = PathBuf::from(INPUT_DIR).join("query_local.bed");
        let exp = PathBuf::from(EXP_DIR).join("basic_example_local.paf");

        let rec_t = read_bed4(&t, None).unwrap();
        let rec_q = read_bed4(&q, None).unwrap();

        // Filter with min score of 6.0
        let res = smith_waterman_affine(&rec_t, &rec_q, 2.0, -1.0, -4.0, -1.0, 6.0).unwrap();
        let paf = PAF::new(res).unwrap();

        let exp_reader = BufReader::new(File::open(exp).unwrap());
        for (line, exp_line) in paf.as_str().zip(exp_reader.lines().map_while(Result::ok)) {
            assert_eq!(line, exp_line)
        }
    }
}
