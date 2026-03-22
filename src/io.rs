use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write, stdout},
    num::NonZeroUsize, path::PathBuf,
};
use flate2::read::MultiGzDecoder;

use eyre::bail;
use itertools::Itertools;

use crate::align::BEDPE;

pub const DEF_NAME_COL: NonZeroUsize = NonZeroUsize::new(4).unwrap();

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BED4 {
    pub chrom: String,
    pub st: u64,
    pub end: u64,
    pub name: String,
}

pub fn read_bed(path: &PathBuf, col_index: Option<NonZeroUsize>) -> eyre::Result<Vec<BED4>> {
    let fh = File::open(path)?;
    let reader = if path.ends_with("gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(fh))) as Box<dyn BufRead>
    } else {
        Box::new(BufReader::new(fh))
    };
    let mut rows = Vec::new();

    // Three for chrom, st, and end.
    let skip_cols = col_index.unwrap_or(DEF_NAME_COL).get() - 3;
    for line in reader.lines() {
        let line = line?;
        let Some((chrom, st, end, other)) = line.splitn(4, '\t').collect_tuple() else {
            bail!("Invalid line {line}")
        };
        // Skip cols and take value as name column
        let Some(name) = other.split('\t').take(skip_cols).next() else {
            bail!("Invalid column for name column: {}", skip_cols + 3)
        };
        let st = st.parse::<u64>()?;
        let end = end.parse::<u64>()?;

        rows.push(BED4 {
            chrom: chrom.to_owned(),
            st,
            end,
            name: name.to_owned(),
        });
    }

    Ok(rows)
}

pub fn write_bedpe(alns: &[BEDPE], path: Option<PathBuf>) -> eyre::Result<()> {
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

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::io::read_bed;

    #[test]
    fn test_read_bed() {
        let res = read_bed(
            &PathBuf::from("test/data/HG008-N_v6.3_chr7_hap2_57312660-64850688.bed.gz"),
            None,
        )
        .unwrap();

        for line in res {
            println!("{line:?}")
        }
    }
}
