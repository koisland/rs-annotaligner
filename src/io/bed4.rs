use flate2::read::MultiGzDecoder;
use std::{
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader},
    num::NonZeroUsize,
    path::Path,
};

use eyre::bail;
use itertools::Itertools;

pub const DEF_NAME_COL: NonZeroUsize = NonZeroUsize::new(4).unwrap();

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BED4 {
    pub chrom: String,
    pub st: u64,
    pub end: u64,
    pub name: String,
}

/// Read a [BED4] file and optionally choose a column index as its name label.
pub fn read_bed4(path: &Path, col_index: Option<NonZeroUsize>) -> eyre::Result<Vec<BED4>> {
    let fh = File::open(path)?;
    let reader = if path.extension().is_some_and(|ext| ext == OsStr::new("gz")) {
        Box::new(BufReader::new(MultiGzDecoder::new(fh))) as Box<dyn BufRead>
    } else {
        Box::new(BufReader::new(fh))
    };
    let mut rows = Vec::new();

    // Four for chrom, st, end, and name.
    // Nothing is skipped by default.
    let skip_cols = col_index.unwrap_or(DEF_NAME_COL).get().saturating_sub(4);
    for line in reader.lines() {
        let line = line?;
        let Some((chrom, st, end, other)) = line.splitn(4, '\t').collect_tuple() else {
            bail!("Invalid line {line}")
        };
        // Skip cols and take value as name column
        let Some(name) = other.split('\t').nth(skip_cols) else {
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

#[cfg(test)]
mod test {
    use std::path::Path;

    use super::{BED4, read_bed4};

    const INFILE: &str = "test/data/input/alternate_col.bed";

    fn expected() -> Vec<BED4> {
        vec![
            BED4 {
                chrom: "chr1".to_owned(),
                st: 1,
                end: 2,
                name: "And".to_owned(),
            },
            BED4 {
                chrom: "chr1".to_owned(),
                st: 2,
                end: 3,
                name: "Bow".to_owned(),
            },
            BED4 {
                chrom: "chr1".to_owned(),
                st: 3,
                end: 4,
                name: "Cow".to_owned(),
            },
        ]
    }

    #[test]
    fn test_read_bed4() {
        let infile = Path::new(INFILE);
        let res = read_bed4(infile, None).unwrap();
        assert_eq!(res, expected())
    }

    #[test]
    fn test_read_bed4_change_col() {
        let infile = Path::new(INFILE);
        let res = read_bed4(infile, Some(5.try_into().unwrap())).unwrap();

        // 4th column is capitalized while 5th is not.
        let expected: Vec<BED4> = expected()
            .into_iter()
            .map(|mut bed| {
                bed.name = bed.name.to_ascii_lowercase();
                bed
            })
            .collect();
        assert_eq!(res, expected)
    }
}
