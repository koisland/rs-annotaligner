use flate2::read::MultiGzDecoder;
use ordered_float::OrderedFloat;
use std::{
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    str::FromStr,
};

use eyre::bail;
use itertools::Itertools;

use crate::{align::CigarOp, io::align::Alignment};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BEDPE {
    pub chrom_1: Option<String>,
    pub chrom_1_st: Option<u64>,
    pub chrom_1_end: Option<u64>,
    pub chrom_1_name: Option<String>,
    pub chrom_2: Option<String>,
    pub chrom_2_st: Option<u64>,
    pub chrom_2_end: Option<u64>,
    pub chrom_2_name: Option<String>,
    pub op: CigarOp,
}

impl BEDPE {
    pub fn header() -> &'static str {
        "#tchrom\ttst\ttend\tqchrom\tqst\tqend\tname\tscore\top\tn_aln"
    }

    pub fn as_str(&self, score: f32, num: usize) -> String {
        match (&self.chrom_1, &self.chrom_2) {
            (Some(chrom_1), Some(chrom_2)) => {
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}~{}\t{}\t{:?}\t{}",
                    chrom_1,
                    self.chrom_1_st.as_ref().unwrap(),
                    self.chrom_1_end.as_ref().unwrap(),
                    chrom_2,
                    self.chrom_2_st.as_ref().unwrap(),
                    self.chrom_2_end.as_ref().unwrap(),
                    self.chrom_1_name.as_ref().unwrap(),
                    self.chrom_2_name.as_ref().unwrap(),
                    score,
                    self.op,
                    num
                )
            }
            (Some(chrom_1), None) => {
                format!(
                    "{}\t{}\t{}\t.\t.\t.\t{}~.\t{}\t{:?}\t{}",
                    chrom_1,
                    self.chrom_1_st.as_ref().unwrap(),
                    self.chrom_1_end.as_ref().unwrap(),
                    self.chrom_1_name.as_ref().unwrap(),
                    score,
                    self.op,
                    num
                )
            }
            (None, Some(chrom_2)) => {
                format!(
                    ".\t.\t.\t{}\t{}\t{}\t.~{}\t{}\t{:?}\t{}",
                    chrom_2,
                    self.chrom_2_st.as_ref().unwrap(),
                    self.chrom_2_end.as_ref().unwrap(),
                    self.chrom_2_name.as_ref().unwrap(),
                    score,
                    self.op,
                    num
                )
            }
            (None, None) => todo!(),
        }
    }
}

#[allow(dead_code)]
/// Read [BEDPE] file for one alignment.
pub fn read_bedpe(path: &Path) -> eyre::Result<Alignment> {
    let fh = File::open(path)?;
    let reader = if path.extension().is_some_and(|ext| ext == OsStr::new("gz")) {
        Box::new(BufReader::new(MultiGzDecoder::new(fh))) as Box<dyn BufRead>
    } else {
        Box::new(BufReader::new(fh))
    };
    let mut records = vec![];
    let mut aln_score: Option<f32> = None;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let Some((
            chrom_1,
            chrom_1_st,
            chrom_1_end,
            chrom_2,
            chrom_2_st,
            chrom_2_end,
            name,
            score,
            op,
            _num,
        )) = line.splitn(10, '\t').collect_tuple()
        else {
            bail!("Invalid line {line}")
        };
        let Some((chrom_1_name, chrom_2_name)) = name.split('~').collect_tuple() else {
            bail!("{line} has name column that contains more than one ~")
        };
        let (chrom_1, chrom_1_st, chrom_1_end, chrom_1_name) = if chrom_1 == "." {
            (None, None, None, None)
        } else {
            (
                Some(chrom_1.to_owned()),
                Some(chrom_1_st.parse()?),
                Some(chrom_1_end.parse()?),
                Some(chrom_1_name.parse()?),
            )
        };
        let (chrom_2, chrom_2_st, chrom_2_end, chrom_2_name) = if chrom_2 == "." {
            (None, None, None, None)
        } else {
            (
                Some(chrom_2.to_owned()),
                Some(chrom_2_st.parse()?),
                Some(chrom_2_end.parse()?),
                Some(chrom_2_name.parse()?),
            )
        };
        let score = score.parse::<f32>()?;
        if aln_score.is_some_and(|s| s != score) {
            bail!("Different score found in records. BEDPE should represent only one alignment.")
        } else if aln_score.is_none() {
            aln_score = Some(score)
        }
        records.push(BEDPE {
            chrom_1,
            chrom_1_st,
            chrom_1_end,
            chrom_1_name,
            chrom_2,
            chrom_2_st,
            chrom_2_end,
            chrom_2_name,
            op: CigarOp::from_str(op)?,
        });
    }

    Alignment::new(records, OrderedFloat(aln_score.unwrap()))
}
