use crate::{
    align::TraceOp,
    io::{bed4::BED4, bedpe::BEDPE},
};

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
    use super::smith_waterman_affine;
    use crate::io::{bed4::read_bed4, bedpe::read_bedpe};
    use std::path::PathBuf;

    const INPUT_DIR: &str = "test/data/input";
    const EXP_DIR: &str = "test/data/output";

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
