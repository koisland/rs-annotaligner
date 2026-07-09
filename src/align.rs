use std::str::FromStr;

use eyre::bail;

use crate::io::BED4;

#[derive(Debug, Clone, Copy)]
pub enum Mode {
    /// Global alignment via Needleman-Wunsch with affine gap penalties
    Global,
    /// Local alignment via Smith-Waterman with affine gap penalties
    Local,
}

impl FromStr for Mode {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "local" => Ok(Mode::Local),
            "global" => Ok(Mode::Global),
            _ => bail!("Invalid mode: {s}"),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum TraceOp {
    /// Match or mismatch at (i, j)
    M,
    /// Gap in query (insertion in target) at (i, j)
    X,
    /// Gap in target (deletion in target) at (i, j)
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

impl Trace {
    pub fn new(op: TraceOp, i: usize, j: usize) -> Self {
        Trace { op, i, j }
    }
}

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

pub fn needleman_wuncsh_affine(
    bed_target: &[BED4],
    bed_query: &[BED4],
    score_match: f32,
    score_mismatch: f32,
    score_gap_open: f32,
    score_gap_ext: f32,
) -> Vec<BEDPE> {
    /*
    Global alignment with affine gap penalties.
    Three DP matrices:
      M[i][j] — best score ending with a match/mismatch at (i,j)
      X[i][j] — best score ending with a gap in seq2 (insert in seq1) at (i,j)
      Y[i][j] — best score ending with a gap in seq1 (delete from seq1) at (i,j)
    */
    let n = bed_target.len();
    let m = bed_query.len();

    let mut M = vec![vec![f32::NEG_INFINITY; m + 1]; n + 1];
    let mut X = vec![vec![f32::NEG_INFINITY; m + 1]; n + 1];
    let mut Y = vec![vec![f32::NEG_INFINITY; m + 1]; n + 1];

    let mut trace_m: Vec<Vec<Option<Trace>>> = vec![vec![None; m + 1]; n + 1];
    let mut trace_x: Vec<Vec<Option<Trace>>> = vec![vec![None; m + 1]; n + 1];
    let mut trace_y: Vec<Vec<Option<Trace>>> = vec![vec![None; m + 1]; n + 1];

    // Cast bool to usize and index rather than if statement.
    let match_mismatch_score = [score_mismatch, score_match];

    // Init top-left corner of mtx
    //    t0 t1 t2
    // q0 0
    // q1
    // q2
    M[0][0] = 0.0;
    for i in 1..n + 1 {
        if i == 1 {
            X[i][0] = score_gap_open + score_gap_ext;
            trace_x[i][0] = Some(Trace {
                op: TraceOp::M,
                i: 1,
                j: 0,
            })
        } else {
            X[i][0] = X[i - 1][0] + score_gap_ext;
            trace_x[i][0] = Some(Trace {
                op: TraceOp::X,
                i: 1,
                j: 0,
            });
        }
    }
    for j in 1..m + 1 {
        if j == 1 {
            Y[0][j] = score_gap_open + score_gap_ext;
            trace_x[0][j] = Some(Trace {
                op: TraceOp::M,
                i: 0,
                j: 1,
            })
        } else {
            Y[0][j] = Y[0][j - 1] + score_gap_ext;
            trace_y[0][j] = Some(Trace {
                op: TraceOp::Y,
                i: 0,
                j: 1,
            });
        }
    }

    // https://cseweb.ucsd.edu/classes/wi12/cse282-a/Lecture03_Ch06_Alignment.pdf
    for i in tqdm::tqdm(1..n + 1) {
        let target_i = &bed_target[i - 1];
        for j in 1..m + 1 {
            let query_j = &bed_query[j - 1];

            let sub_score = match_mismatch_score[usize::from(target_i.name == query_j.name)];
            // Diagonal (* = current, (#) = evaluating)
            //    t0 t1 t2
            // q0 (#)
            // q1    *
            // q2
            let cand_m = [
                (Y[i - 1][j - 1] + sub_score, Trace::new(TraceOp::Y, 1, 1)),
                (X[i - 1][j - 1] + sub_score, Trace::new(TraceOp::X, 1, 1)),
                (M[i - 1][j - 1] + sub_score, Trace::new(TraceOp::M, 1, 1)),
            ];
            let max_cand_m = cand_m
                .into_iter()
                .max_by(|a, b| a.0.total_cmp(&b.0))
                .unwrap();
            M[i][j] = max_cand_m.0;
            trace_m[i][j] = Some(max_cand_m.1);

            // Left (* = current, (#) = evaluating)
            //    t0  t1 t2
            // q0
            // q1 (#) *
            // q2
            let cand_x = [
                // Extend
                (X[i - 1][j] + score_gap_ext, Trace::new(TraceOp::X, 1, 0)),
                // Open
                (
                    M[i - 1][j] + score_gap_open + score_gap_ext,
                    Trace::new(TraceOp::M, 1, 0),
                ),
            ];
            let max_cand_x = cand_x
                .into_iter()
                .max_by(|a, b| a.0.total_cmp(&b.0))
                .unwrap();
            X[i][j] = max_cand_x.0;
            trace_x[i][j] = Some(max_cand_x.1);

            // Top (* = current, (#) = evaluating)
            //    t0  t1  t2
            // q0     (#)
            // q1     *
            // q2
            let cand_y = [
                // Extend
                (Y[i][j - 1] + score_gap_ext, Trace::new(TraceOp::Y, 0, 1)),
                // Open
                (
                    M[i][j - 1] + score_gap_open + score_gap_ext,
                    Trace::new(TraceOp::M, 0, 1),
                ),
            ];
            // Max by returns last not first.
            let max_cand_y = cand_y
                .into_iter()
                .max_by(|a, b| a.0.total_cmp(&b.0))
                .unwrap();
            Y[i][j] = max_cand_y.0;
            trace_y[i][j] = Some(max_cand_y.1);

            // eprintln!(
            //     "({i}, {j})\n\t{:?}\n\t\t{:?}\n\t{:?}\n\t\t{:?}\n\t{:?}\n\t\t{:?}",
            //     max_cand_m, cand_m, max_cand_x, cand_x, max_cand_y, cand_y
            // );
        }
    }

    // End state
    let end_scores = [
        (Y[n][m], TraceOp::Y),
        (X[n][m], TraceOp::X),
        (M[n][m], TraceOp::M),
    ];
    let (_, mut state) = end_scores
        .into_iter()
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap();
    // eprintln!("({n}, {m})\n\t{state:?}\n\t\t{end_scores:?}");

    // Traceback
    let (mut i, mut j) = (n, m);
    let mut alns: Vec<BEDPE> = Vec::with_capacity(std::cmp::max(n, m));
    // for (i, row) in trace_m.iter().enumerate() {
    //     eprintln!("M{i}: {row:?}");
    // }
    // for (i, row) in trace_x.iter().enumerate() {
    //     eprintln!("X{i}: {row:?}");
    // }
    // for (i, row) in trace_y.iter().enumerate() {
    //     eprintln!("Y{i}: {row:?}");
    // }
    while i > 0 || j > 0 {
        let (aln_target, aln_query, trace) = match state {
            TraceOp::M => (
                Some(&bed_target[i - 1]),
                Some(&bed_query[j - 1]),
                trace_m[i][j].unwrap(),
            ),
            TraceOp::X => (Some(&bed_target[i - 1]), None, trace_x[i][j].unwrap()),
            TraceOp::Y => (None, Some(&bed_query[j - 1]), trace_y[i][j].unwrap()),
        };
        // eprintln!("({i}, {j}) {trace:?}");
        let bedpe = match (aln_target, aln_query) {
            (Some(target), Some(query)) => BEDPE {
                chrom_1: Some(target.chrom.to_owned()),
                chrom_1_st: Some(target.st),
                chrom_1_end: Some(target.end),
                chrom_1_name: Some(target.name.to_owned()),
                chrom_2: Some(query.chrom.to_owned()),
                chrom_2_st: Some(query.st),
                chrom_2_end: Some(query.end),
                chrom_2_name: Some(query.name.to_owned()),
            },
            (Some(target), None) => BEDPE {
                chrom_1: Some(target.chrom.to_owned()),
                chrom_1_st: Some(target.st),
                chrom_1_end: Some(target.end),
                chrom_1_name: Some(target.name.to_owned()),
                chrom_2: None,
                chrom_2_st: None,
                chrom_2_end: None,
                chrom_2_name: None,
            },
            (None, Some(query)) => BEDPE {
                chrom_1: None,
                chrom_1_st: None,
                chrom_1_end: None,
                chrom_1_name: None,
                chrom_2: Some(query.chrom.to_owned()),
                chrom_2_st: Some(query.st),
                chrom_2_end: Some(query.end),
                chrom_2_name: Some(query.name.to_owned()),
            },
            (None, None) => unreachable!(),
        };
        alns.push(bedpe);
        i -= trace.i;
        j -= trace.j;
        state = trace.op
    }

    alns.reverse();
    alns
}

pub fn smith_waterman_affine(
    bed_target: &[BED4],
    bed_query: &[BED4],
    score_match: f32,
    score_mismatch: f32,
    score_gap_open: f32,
    score_gap_ext: f32,
) -> Vec<BEDPE> {
    // https://informatika.stei.itb.ac.id/~rinaldi.munir/Stmik/2024-2025/Makalah2025/Makalah-IF2211-Strategi-Algoritma-2025%20(98).pdf
    // Adapted from https://github.com/varel183/Makalah_STIMA_13523008/blob/main/src/main.py
    let n = bed_target.len();
    let m = bed_query.len();

    let mut M = vec![vec![0.0; m + 1]; n + 1];
    let mut X = vec![vec![0.0; m + 1]; n + 1];
    let mut Y = vec![vec![0.0; m + 1]; n + 1];

    let mut max_pos: ((usize, usize), f32) = ((0, 0), 0.0);
    let match_mismatch_score = [score_mismatch, score_match];

    // TODO: Max heap of size n to store positions for traceback.

    for i in tqdm::tqdm(1..n + 1) {
        let target_i = &bed_target[i - 1];
        for j in 1..m + 1 {
            let query_j = &bed_query[j - 1];
            let sub_score = match_mismatch_score[usize::from(target_i.name == query_j.name)];

            // Recurrence relation
            let cand_m = [
                // 0 is provided if all scores are negative
                // to avoid extending bad alignments
                0.0,
                M[i - 1][j - 1] + sub_score,
                X[i - 1][j - 1] + sub_score,
                Y[i - 1][j - 1] + sub_score,
            ];
            let max_cand_m = cand_m.into_iter().max_by(|a, b| a.total_cmp(b)).unwrap();
            M[i][j] = max_cand_m;

            let cand_x = [
                // Transition to gap state
                M[i - 1][j] + score_gap_open,
                // Stay within gap
                X[i - 1][j] + score_gap_ext,
            ];
            let max_cand_x = cand_x.into_iter().max_by(|a, b| a.total_cmp(b)).unwrap();
            X[i][j] = max_cand_x;

            let cand_y = [M[i][j - 1] + score_gap_open, X[i][j - 1] + score_gap_ext];
            let max_cand_y = cand_y.into_iter().max_by(|a, b| a.total_cmp(b)).unwrap();
            Y[i][j] = max_cand_y;

            let current_score = [M[i][j], X[i][j], Y[i][j]]
                .into_iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap();

            max_pos =
                std::cmp::max_by(max_pos, ((i, j), current_score), |a, b| a.1.total_cmp(&b.1));
        }
    }

    let ((mut i, mut j), _score) = max_pos;
    let score_m = M[i][j];
    let score_x = X[i][j];
    let score_y = Y[i][j];

    let mut current_op = if score_m >= score_x && score_m >= score_y {
        TraceOp::M
    } else if score_x >= score_y {
        TraceOp::X
    } else {
        TraceOp::Y
    };

    let mut alns: Vec<BEDPE> = Vec::with_capacity(std::cmp::max(n, m));
    while i > 0 && j > 0 {
        let score_m = M[i][j];
        let score_x = X[i][j];
        let score_y = Y[i][j];
        let max_score = [score_m, score_x, score_y]
            .into_iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();

        if max_score <= 0.0 {
            break;
        }

        let target_i = &bed_target[i - 1];
        let query_j = &bed_query[j - 1];
        let sub_score = match_mismatch_score[usize::from(target_i.name == query_j.name)];

        match current_op {
            TraceOp::M => {
                alns.push(BEDPE {
                    chrom_1: Some(target_i.chrom.to_owned()),
                    chrom_1_st: Some(target_i.st),
                    chrom_1_end: Some(target_i.end),
                    chrom_1_name: Some(target_i.name.to_owned()),
                    chrom_2: Some(query_j.chrom.to_owned()),
                    chrom_2_st: Some(query_j.st),
                    chrom_2_end: Some(query_j.end),
                    chrom_2_name: Some(query_j.name.to_owned()),
                });
                let next_score_m = M[i - 1][j - 1];
                let next_score_x = X[i - 1][j - 1];
                if score_m == next_score_m + sub_score {
                    current_op = TraceOp::M
                } else if score_m == next_score_x + sub_score {
                    current_op = TraceOp::X
                } else {
                    current_op = TraceOp::Y
                }
                i -= 1;
                j -= 1;
            }
            TraceOp::X => {
                alns.push(BEDPE {
                    chrom_1: Some(target_i.chrom.to_owned()),
                    chrom_1_st: Some(target_i.st),
                    chrom_1_end: Some(target_i.end),
                    chrom_1_name: Some(target_i.name.to_owned()),
                    chrom_2: None,
                    chrom_2_st: None,
                    chrom_2_end: None,
                    chrom_2_name: None,
                });
                let next_score_x = X[i - 1][j];
                if score_x == next_score_x + score_gap_open {
                    current_op = TraceOp::M
                } else {
                    current_op = TraceOp::X
                }
                i -= 1;
            }
            TraceOp::Y => {
                alns.push(BEDPE {
                    chrom_1: None,
                    chrom_1_st: None,
                    chrom_1_end: None,
                    chrom_1_name: None,
                    chrom_2: Some(query_j.chrom.to_owned()),
                    chrom_2_st: Some(query_j.st),
                    chrom_2_end: Some(query_j.end),
                    chrom_2_name: Some(query_j.name.to_owned()),
                });
                let next_score_y = Y[i][j - 1];
                if score_y == next_score_y + score_gap_open {
                    current_op = TraceOp::M
                } else {
                    current_op = TraceOp::Y
                }
                j -= 1;
            }
        }
    }
    alns.reverse();
    alns
}

#[cfg(test)]
mod test {
    use crate::{
        align::{needleman_wuncsh_affine, smith_waterman_affine},
        io::{read_bed4, read_bedpe},
    };
    use std::path::PathBuf;

    const INPUT_DIR: &str = "test/data/input";
    const EXP_DIR: &str = "test/data/output";

    #[test]
    fn test_global() {
        let t =
            PathBuf::from(INPUT_DIR).join("HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz");
        let q = PathBuf::from(INPUT_DIR)
            .join("HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz");
        let exp = PathBuf::from(EXP_DIR).join("HG008-TN_chr6_chr7_fusion.bed.gz");

        let rec_t = read_bed4(&t, None).unwrap();
        let rec_q = read_bed4(&q, None).unwrap();

        let res = needleman_wuncsh_affine(&rec_t, &rec_q, 2.0, -1.0, -4.0, -1.0);
        let exp_res = read_bedpe(&exp).unwrap();
        assert_eq!(res, exp_res)
    }

    #[test]
    fn test_global_small() {
        let t = PathBuf::from(INPUT_DIR).join("target.bed");
        let q = PathBuf::from(INPUT_DIR).join("query.bed");
        let exp = PathBuf::from(EXP_DIR).join("basic_example.bedpe");

        let rec_t = read_bed4(&t, None).unwrap();
        let rec_q = read_bed4(&q, None).unwrap();

        let res = needleman_wuncsh_affine(&rec_t, &rec_q, 2.0, -1.0, -4.0, -1.0);
        let exp_res = read_bedpe(&exp).unwrap();
        assert_eq!(res, exp_res)
    }

    #[test]
    fn test_local_small() {
        let t = PathBuf::from(INPUT_DIR).join("target_local.bed");
        let q = PathBuf::from(INPUT_DIR).join("query_local.bed");
        let exp = PathBuf::from(EXP_DIR).join("basic_example_local.bedpe");

        let rec_t = read_bed4(&t, None).unwrap();
        let rec_q = read_bed4(&q, None).unwrap();

        let res = smith_waterman_affine(&rec_t, &rec_q, 2.0, -1.0, -4.0, -1.0);
        let exp_res = read_bedpe(&exp).unwrap();
        assert_eq!(res, exp_res)
    }
}
