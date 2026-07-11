use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write, stdout},
    path::Path,
};

use digits_iterator::DigitsExtension;
use eyre::{ContextCompat, bail};
use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::{
    align::CigarOp,
    io::{OutputType, bedpe::BEDPE, paf::PAF},
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// Records in alignment
    /// * `nr` tag (`i`) - Number of BED4 records
    pub records: Vec<BEDPE>,
    /// DP alignment score
    /// * `AS` tag (`i`)
    pub score: OrderedFloat<f32>,
    /// CIGAR string
    /// * `cg` tag (`Z`)
    pub cigar: String,
    /// Number of matching bases (Based on intervals)
    pub num_matches: usize,
    /// Number of bases (Based on intervals)
    pub num_bases: usize,
    /// Number of mismatches and alignment gaps
    /// * `NM` tag (`i`)
    pub num_mismatches_gaps: usize,
    /// Gap compressed divergence
    /// * `de` tag (f)
    pub gap_cmp_divergence: OrderedFloat<f32>,
    /// MAPQ (From 0 to 60)
    /// * Based on gap-compressed sequence divergence and max DP score.
    pub mapq: Option<u8>,
    /// Target interval
    pub titv: (String, u64, u64),
    /// Query interval
    pub qitv: (String, u64, u64),
    /// Target annotation gap percentage
    /// * `tg` tag (`f`)
    pub tgap_perc: OrderedFloat<f32>,
    /// Query annotation gap percentage
    /// * `qg` tag (`f`)
    pub qgap_perc: OrderedFloat<f32>,
}

impl Alignment {
    pub fn new(records: Vec<BEDPE>, score: OrderedFloat<f32>) -> eyre::Result<Self> {
        if records.is_empty() {
            bail!("No valid records")
        }
        let mut cigar = String::new();
        let mut op_counts: HashMap<CigarOp, (u64, u64)> = HashMap::new();
        // NOTE: This does not check if multiple alignments are present.
        let (mut min_tst, mut max_tend, mut min_qst, mut max_qend) =
            (u64::MAX, 0u64, u64::MAX, 0u64);
        // Keep track of the gap size.
        let (mut tgap_len, mut qgap_len) = (0, 0);
        let (mut prev_tst, mut prev_tend, mut prev_qst, mut prev_qend) = (None, None, None, None);
        let (mut tchrom, mut qchrom) = (None, None);

        // Group by target chrom, query chrom, and operation.
        for ((chrom_1, chrom_2, op), recs) in &records
            .iter()
            .chunk_by(|rec| (&rec.chrom_1, &rec.chrom_2, rec.op))
        {
            if chrom_1.as_ref().is_some_and(|_| tchrom.is_none()) {
                tchrom = chrom_1.as_ref().map(|c| c.to_owned())
            }
            if chrom_2.as_ref().is_some_and(|_| qchrom.is_none()) {
                qchrom = chrom_2.as_ref().map(|c| c.to_owned())
            }
            let (tst, tend, qst, qend) =
                recs.into_iter()
                    .fold((u64::MAX, 0u64, u64::MAX, 0u64), |acc, x| {
                        (
                            std::cmp::min(acc.0, x.chrom_1_st.unwrap_or(acc.0)),
                            std::cmp::max(acc.1, x.chrom_1_end.unwrap_or(acc.1)),
                            std::cmp::min(acc.2, x.chrom_2_st.unwrap_or(acc.2)),
                            std::cmp::max(acc.3, x.chrom_2_end.unwrap_or(acc.3)),
                        )
                    });
            // Track and update min and max
            min_tst = std::cmp::min(tst, min_tst);
            max_tend = std::cmp::max(tend, max_tend);
            min_qst = std::cmp::min(qst, min_qst);
            max_qend = std::cmp::max(qend, max_qend);

            // Record gap size, if any.
            // If has coordinate, set to prev.
            if let (Some(prev_tst), Some(prev_tend)) = (prev_tst, prev_tend) {
                // No overlap
                let valid_titvs = tst != u64::MAX && prev_tst != u64::MAX;
                if !(tst < prev_tend && tend > prev_tst) && valid_titvs {
                    let tdiff = std::cmp::max(tst, prev_tst) - std::cmp::min(tend, prev_tend);
                    // eprintln!("{prev_tst},{prev_tend}|{tst},{tend}|t{tdiff}");
                    tgap_len += tdiff;
                }
            }
            if let (Some(prev_qst), Some(prev_qend)) = (prev_qst, prev_qend) {
                // No overlap
                let valid_qitvs = qst != u64::MAX && prev_qst != u64::MAX;
                if !(qst < prev_qend && tend > prev_qst) && valid_qitvs {
                    let qdiff = std::cmp::max(qst, prev_qst) - std::cmp::min(qend, prev_qend);
                    // eprintln!("{qst},{qend}|{prev_qst},{prev_qend}|q{qdiff}");
                    qgap_len += qdiff;
                }
            }
            prev_tst = Some(tst);
            prev_tend = Some(tend);
            prev_qst = Some(qst);
            prev_qend = Some(qend);

            let length = match op {
                CigarOp::Match | CigarOp::Mismatch | CigarOp::Deletion => tend - tst,
                CigarOp::Insertion => qend - qst,
            };
            // Count for estimated interval sequence divergence.
            op_counts
                .entry(op)
                .and_modify(|(op_length, count)| {
                    *op_length += length;
                    *count += 1;
                })
                .or_insert((length, 1));

            let char_op: char = op.into();
            // Iterate thru digits and add one char at time.
            for digit in length.digits() {
                cigar.push(
                    char::from_digit(digit as u32, 10).with_context(|| {
                        let chrom_1 = chrom_1.as_deref().unwrap_or(".");
                        let chrom_2 = chrom_2.as_deref().unwrap_or(".");
                        format!("Invalid digit {digit} in {length} for records between {chrom_1}{tst}-{tend} and {chrom_2}{qst}-{qend}")
                    })?,
                );
            }
            cigar.push(char_op);
        }

        // If not valid alignment
        if qchrom.is_none() || tchrom.is_none() {
            bail!("{records:?}")
        }

        // From https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
        let matches = op_counts
            .get(&CigarOp::Match)
            .map_or(0.0, |(l, _)| *l as f32);
        let mismatches = op_counts
            .get(&CigarOp::Mismatch)
            .map_or(0.0, |(l, _)| *l as f32);
        // Only consider counts.
        let (deletions, n_deletions) = op_counts
            .get(&CigarOp::Deletion)
            .map_or((0, 0.0), |(l, c)| (*l as usize, *c as f32));
        let (insertions, n_insertions) = op_counts
            .get(&CigarOp::Insertion)
            .map_or((0, 0.0), |(l, c)| (*l as usize, *c as f32));
        let n_cols = op_counts.values().map(|v| v.0 as f32).sum::<f32>();
        // Gap-compressed identity
        let de = OrderedFloat(matches / (n_cols - n_deletions - n_insertions));

        // Track length of alignment and the gaps in between
        let titv = (tchrom.expect("No tchrom"), min_tst, max_tend);
        let tlen = titv.2 - titv.1;
        let qitv = (qchrom.expect("No qchrom"), min_qst, max_qend);
        let qlen = qitv.2 - qitv.1;

        Ok(Self {
            records,
            score,
            cigar,
            num_bases: n_cols as usize,
            num_matches: matches as usize,
            num_mismatches_gaps: mismatches as usize + deletions + insertions,
            gap_cmp_divergence: de,
            titv,
            qitv,
            tgap_perc: OrderedFloat(tgap_len as f32 / tlen as f32),
            qgap_perc: OrderedFloat(qgap_len as f32 / qlen as f32),
            mapq: None,
        })
    }

    pub(crate) fn as_paf_row(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}\tnr:i:{}\tAS:i:{}\tNM:i:{}\tde:f:{}\ttg:f:{}\tqg:f:{}",
            // query chrom
            self.qitv.0,
            // query length
            // We don't know this.
            self.qitv.2 - self.qitv.1,
            // query start
            self.qitv.1,
            // query end
            self.qitv.2,
            // target chrom
            self.titv.0,
            // target length
            // We don't know this.
            self.titv.2 - self.titv.1,
            // target start
            self.titv.1,
            // target end
            self.titv.2,
            // number of matching bases (based on interval length)
            self.num_matches,
            // number of base including gaps.
            self.num_bases,
            // MAPQ
            self.mapq.unwrap_or(255),
            // cg:Z: - cigar
            self.cigar,
            // nr:i: - number of records
            self.records.len(),
            // AS:i: - DP score
            self.score,
            // NM:i: - number of mismatches/gaps
            self.num_mismatches_gaps,
            // de:f: - gap-compressed sequence divergence
            self.gap_cmp_divergence,
            // tg:f: - target sequence annotation gap percent
            self.tgap_perc,
            // qg:f: - query sequence annotation gap percent
            self.qgap_perc,
        )
    }
}

impl PartialOrd for Alignment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Alignment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score.cmp(&other.score)
    }
}

pub fn write_alignments(
    alignments: Vec<Alignment>,
    output_type: OutputType,
    path: Option<&Path>,
) -> eyre::Result<()> {
    let mut writer = if let Some(outfile) = path {
        Box::new(BufWriter::new(File::create(outfile)?)) as Box<dyn Write>
    } else {
        Box::new(BufWriter::new(stdout()))
    };
    match output_type {
        OutputType::PAF => {
            // Create PAF and calculate MAPQ
            let paf = PAF::new(alignments)?;
            for line in paf.as_str() {
                writeln!(&mut writer, "{line}")?;
            }
        }
        OutputType::BEDPE => {
            writeln!(&mut writer, "{}", BEDPE::header())?;
            for (i, aln) in alignments.into_iter().enumerate() {
                for record in aln.records {
                    writeln!(&mut writer, "{}", record.as_str(*aln.score, i + 1))?;
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::{
        align::global::needleman_wuncsh_affine,
        io::bed4::{BED4, read_bed4},
    };

    const INPUT_DIR: &str = "test/data/input";

    // Using sequence since easier.
    // From https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
    const REF: &str = "CCAGTGTGGCCGATaCCCcagGTtgGCACGCATCGTTGCCTTGGTAAGC";
    const QRY: &str = "CCAGTGTGGCCGATgCCCGTGCtACGCATCGTTGCCTTGGTAAGC";

    fn sequence_to_bed(seq: &str) -> Vec<BED4> {
        seq.chars()
            .enumerate()
            .map(|(i, c)| BED4 {
                chrom: "chr1".to_owned(),
                st: i as u64,
                end: (i + 1) as u64,
                name: String::from(c),
            })
            .collect()
    }

    #[test]
    fn test_alignment_metrics() {
        let target = sequence_to_bed(REF);
        let query = sequence_to_bed(QRY);

        let res = needleman_wuncsh_affine(&target, &query, 1.0, -2.0, -2.0, -1.0).unwrap();
        // de
        assert_eq!(res.gap_cmp_divergence, 0.9148936);
        // cigar
        assert_eq!(res.cigar, "14=1X3=3D2=2D2=1I22=");
        // No gaps in interval because made from sequence
        assert_eq!(res.titv, (String::from("chr1"), 0, 49));
        assert_eq!(res.qitv, (String::from("chr1"), 0, 45));
        assert_eq!(res.tgap_perc, 0.0);
        assert_eq!(res.qgap_perc, 0.0);
    }

    #[test]
    fn test_gap_annot_alignment() {
        let t = PathBuf::from(INPUT_DIR).join("target_gap.bed");
        let q = PathBuf::from(INPUT_DIR).join("query_gap.bed");

        let rec_t = read_bed4(&t, None).unwrap();
        let rec_q = read_bed4(&q, None).unwrap();

        let res = needleman_wuncsh_affine(&rec_t, &rec_q, 2.0, -1.0, -4.0, -1.0).unwrap();

        assert_eq!(res.titv, (String::from("chr1"), 1, 20));
        assert_eq!(res.qitv, (String::from("chr1"), 1, 9));
        // Gap in target BED annotation aligned over but recorded.
        // chr1	5	6	Anger
        // chr1	10	11	Creek
        assert_eq!(res.tgap_perc, 0.21052632);
        assert_eq!(res.qgap_perc, 0.0);
    }
}
