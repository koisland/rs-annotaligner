use flate2::read::MultiGzDecoder;
use std::{
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write, stdout},
    path::Path,
};

use eyre::bail;
use itertools::Itertools;

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
}

impl BEDPE {
    pub fn as_row(&self) -> String {
        match (&self.chrom_1, &self.chrom_2) {
            (Some(chrom_1), Some(chrom_2)) => {
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}~{}\t{}",
                    chrom_1,
                    self.chrom_1_st.as_ref().unwrap(),
                    self.chrom_1_end.as_ref().unwrap(),
                    chrom_2,
                    self.chrom_2_st.as_ref().unwrap(),
                    self.chrom_2_end.as_ref().unwrap(),
                    self.chrom_1_name.as_ref().unwrap(),
                    self.chrom_2_name.as_ref().unwrap(),
                    if self.chrom_1_name == self.chrom_2_name {
                        "Match"
                    } else {
                        "Mismatch"
                    },
                )
            }
            (Some(chrom_1), None) => {
                format!(
                    "{}\t{}\t{}\t.\t.\t.\t{}~.\tDeletion",
                    chrom_1,
                    self.chrom_1_st.as_ref().unwrap(),
                    self.chrom_1_end.as_ref().unwrap(),
                    self.chrom_1_name.as_ref().unwrap(),
                )
            }
            (None, Some(chrom_2)) => {
                format!(
                    ".\t.\t.\t{}\t{}\t{}\t.~{}\tInsertion",
                    chrom_2,
                    self.chrom_2_st.as_ref().unwrap(),
                    self.chrom_2_end.as_ref().unwrap(),
                    self.chrom_2_name.as_ref().unwrap(),
                )
            }
            (None, None) => todo!(),
        }
    }
}

#[allow(dead_code)]
/// Read [BEDPE] file.
pub fn read_bedpe(path: &Path) -> eyre::Result<Vec<BEDPE>> {
    let fh = File::open(path)?;
    let reader = if path.extension().is_some_and(|ext| ext == OsStr::new("gz")) {
        Box::new(BufReader::new(MultiGzDecoder::new(fh))) as Box<dyn BufRead>
    } else {
        Box::new(BufReader::new(fh))
    };

    let mut records = vec![];
    for line in reader.lines() {
        let line = line?;
        let Some((
            chrom_1,
            chrom_1_st,
            chrom_1_end,
            chrom_2,
            chrom_2_st,
            chrom_2_end,
            name,
            _score,
        )) = line.splitn(8, '\t').collect_tuple()
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
        records.push(BEDPE {
            chrom_1,
            chrom_1_st,
            chrom_1_end,
            chrom_1_name,
            chrom_2,
            chrom_2_st,
            chrom_2_end,
            chrom_2_name,
        });
    }

    Ok(records)
}

/// Write [BEDPE] alignments to a file or stdout.
pub fn write_bedpe(alns: &[BEDPE], path: Option<&Path>) -> eyre::Result<()> {
    let mut writer = if let Some(outfile) = path {
        Box::new(BufWriter::new(File::create(outfile)?)) as Box<dyn Write>
    } else {
        Box::new(BufWriter::new(stdout()))
    };
    for record in alns {
        writeln!(&mut writer, "{}", record.as_row())?;
    }
    Ok(())
}
