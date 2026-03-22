use crate::io::BED4;

#[derive(Debug, Clone, Copy)]
pub enum TraceOp {
    /// Match or mismatch at (i, j)
    MatchMismatch,
    /// Gap in query (insertion in target) at (i, j)
    Insertion,
    /// Gap in target (deletion in target) at (i, j)
    Deletion,
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

pub struct BEDPE {
    chrom_1: Option<String>,
    chrom_1_st: Option<u64>,
    chrom_1_end: Option<u64>,
    chrom_1_name: Option<String>,
    chrom_2: Option<String>,
    chrom_2_st: Option<u64>,
    chrom_2_end: Option<u64>,
    chrom_2_name: Option<String>,
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
                    if self.chrom_1_name == self.chrom_2_name { "Match" } else { "Mismatch" },
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
            },
            (None, Some(chrom_2)) => {
                format!(
                    ".\t.\t.\t{}\t{}\t{}\t.~{}\tInsertion",
                    chrom_2,
                    self.chrom_2_st.as_ref().unwrap(),
                    self.chrom_2_end.as_ref().unwrap(),
                    self.chrom_2_name.as_ref().unwrap(),
                )
            },
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

    let mut trace_M: Vec<Vec<Option<Trace>>> = vec![vec![None; m + 1]; n + 1];
    let mut trace_X: Vec<Vec<Option<Trace>>> = vec![vec![None; m + 1]; n + 1];
    let mut trace_Y: Vec<Vec<Option<Trace>>> = vec![vec![None; m + 1]; n + 1];

    // Cast bool to usize and index rather than if statement.
    let match_mismatch_score = [score_mismatch, score_match];

    // Init top-left corner of mtx
    //    t0 t1 t2
    // q0 0
    // q1
    // q2
    M[0][0] = 0.0;
    //    t0 t1 t2
    // q0 0  sg
    // q1
    // q2
    X[1][0] = score_gap_open + score_gap_ext;
    trace_X[1][0] = Some(Trace {
        op: TraceOp::MatchMismatch,
        i: 1,
        j: 0,
    });
    //    t0 t1 t2
    // q0 0
    // q1 sg
    // q2
    Y[0][1] = score_gap_open + score_gap_ext;
    trace_Y[0][1] = Some(Trace {
        op: TraceOp::MatchMismatch,
        i: 0,
        j: 1,
    });

    for i in 2..n + 1 {
        //    t0 t1 t2
        // q0 0  sg sg+sg_pre
        // q1
        // q2
        X[i][0] = X[i - 1][0] + score_gap_ext;
        trace_X[i][0] = Some(Trace {
            op: TraceOp::Insertion,
            i: 1,
            j: 0,
        });
    }
    for j in 2..m + 1 {
        //    t0 t1 t2
        // q0 0
        // q1 sg
        // q2 sg+sg_pre
        Y[0][j] = Y[0][j - 1] + score_gap_ext;
        trace_Y[0][j] = Some(Trace {
            op: TraceOp::Deletion,
            i: 0,
            j: 1,
        });
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
            let cand_M = [
                (
                    M[i - 1][j - 1] + sub_score,
                    Trace::new(TraceOp::MatchMismatch, 1, 1),
                ),
                (
                    X[i - 1][j - 1] + sub_score,
                    Trace::new(TraceOp::Insertion, 1, 1),
                ),
                (
                    Y[i - 1][j - 1] + sub_score,
                    Trace::new(TraceOp::Deletion, 1, 1),
                ),
            ];
            let cand_M = cand_M
                .into_iter()
                .max_by(|a, b| a.0.total_cmp(&b.0))
                .unwrap();
            M[i][j] = cand_M.0;
            trace_M[i][j] = Some(cand_M.1);

            // Left (* = current, (#) = evaluating)
            //    t0  t1 t2
            // q0
            // q1 (#) *
            // q2
            let cand_X = std::cmp::max_by(
                // Open
                (
                    M[i - 1][j] + score_gap_open + score_gap_ext,
                    Trace::new(TraceOp::MatchMismatch, 1, 0),
                ),
                // Extend
                (
                    X[i - 1][j] + score_gap_ext,
                    Trace::new(TraceOp::Insertion, 1, 0),
                ),
                |a, b| a.0.total_cmp(&b.0),
            );
            X[i][j] = cand_X.0;
            trace_X[i][j] = Some(cand_X.1);

            // Top (* = current, (#) = evaluating)
            //    t0  t1  t2
            // q0     (#)
            // q1     *
            // q2
            let cand_Y = std::cmp::max_by(
                // Open
                (
                    M[i][j - 1] + score_gap_open + score_gap_ext,
                    Trace::new(TraceOp::MatchMismatch, 0, 1),
                ),
                // Extend
                (
                    Y[i][j - 1] + score_gap_ext,
                    Trace::new(TraceOp::Deletion, 0, 1),
                ),
                |a, b| a.0.total_cmp(&b.0),
            );
            Y[i][j] = cand_Y.0;
            trace_Y[i][j] = Some(cand_Y.1);
        }
    }

    // End state
    let end_scores = [
        (M[n][m], TraceOp::MatchMismatch),
        (X[n][m], TraceOp::Deletion),
        (Y[n][m], TraceOp::Insertion),
    ];
    let (_, mut state) = end_scores
        .into_iter()
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap();

    // Traceback
    let (mut i, mut j) = (n, m);
    let mut alns: Vec<BEDPE> = Vec::with_capacity(std::cmp::max(n, m));
    while i > 0 || j > 0 {
        let (aln_target, aln_query, trace) = match state {
            TraceOp::MatchMismatch => (
                Some(&bed_target[i - 1]),
                Some(&bed_query[j - 1]),
                trace_M[i][j].unwrap(),
            ),
            TraceOp::Insertion => (Some(&bed_target[i - 1]), None, trace_X[i][j].unwrap()),
            TraceOp::Deletion => (None, Some(&bed_query[j - 1]), trace_Y[i][j].unwrap()),
        };
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
) {
    todo!()
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::{align::needleman_wuncsh_affine, io::read_bed};

    #[test]
    fn test_global() {
        let t = PathBuf::from("test/data/input/HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz");
        let q = PathBuf::from("test/data/input/HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz");
        let rec_t = read_bed(&t, None).unwrap();
        let rec_q = read_bed(&q, None).unwrap();
        needleman_wuncsh_affine(&rec_t, &rec_q, 2.0, -1.0, -4.0, -1.0);
    }
}
