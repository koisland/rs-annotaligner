"""
Global alignment (Needleman–Wunsch) of annotations from two BED files.

Inputs:
  - Two BED files. Label is in column 4.
Outputs (defaults; override with CLI):
  - output_alignment.tsv  : aligned labels (two columns, '-' for gaps)
  - output_gaps.tsv       : runs of deletions/insertions/mismatches with coordinates

Scoring defaults (override with CLI):
  --match 2  --mismatch -1  --gap-open -4  --gap-ext -1

Author: Fedor Ryabov fedorrik1@gmail.com
"""

from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import sys
import pandas as pd
from tqdm import tqdm


# ----------------------------
# I/O: read BED into arrays
# ----------------------------
def read_bed(path: str | Path, col_index: int = 3):
    """
    Read BED (tab-delimited) and return three parallel lists:
      labels: list[str]    (column `col_index`, default 3 -> 4th column)
      starts: list[int]    (column 2; 0-based start as in BED)
      ends:   list[int]    (column 3; 1-based end as in BED)
    Lines starting with '#' ignored.
    """
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, engine="python")
    if df.shape[1] <= col_index or df.shape[1] < 4:
        raise ValueError(f"{path} has fewer than 4 columns.")
    # convert start/end to int, keep as-is if missing -> raise
    starts = df[1].astype(int).tolist()
    ends   = df[2].astype(int).tolist()
    labels = df[col_index].astype(str).tolist()
    return labels, starts, ends


# ----------------------------
# Alignment core
# ----------------------------
def sub_score(a: str, b: str, match: float, mismatch: float) -> float:
    return match if a == b else mismatch


def needleman_wunsch_affine(
    seq1: List[str],
    seq2: List[str],
    match: float = 2.0,
    mismatch: float = -1.0,
    gap_open: float = -4.0,
    gap_ext: float = -1.0,
) -> Tuple[List[str], List[str]]:
    """
    Global alignment with affine gap penalties.
    Three DP matrices:
      M[i][j] — best score ending with a match/mismatch at (i,j)
      X[i][j] — best score ending with a gap in seq2 (insert in seq1) at (i,j)
      Y[i][j] — best score ending with a gap in seq1 (delete from seq1) at (i,j)
    """
    n, m = len(seq1), len(seq2)
    NEG_INF = -math.inf

    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    trace_M = [[None] * (m + 1) for _ in range(n + 1)]
    trace_X = [[None] * (m + 1) for _ in range(n + 1)]
    trace_Y = [[None] * (m + 1) for _ in range(n + 1)]

    # Init
    M[0][0] = 0.0
    for i in range(1, n + 1):
        if i == 1:
            X[i][0] = gap_open + gap_ext
            trace_X[i][0] = ('M', 1, 0)
        else:
            X[i][0] = X[i - 1][0] + gap_ext
            trace_X[i][0] = ('X', 1, 0)
    for j in range(1, m + 1):
        if j == 1:
            Y[0][j] = gap_open + gap_ext
            trace_Y[0][j] = ('M', 0, 1)
        else:
            Y[0][j] = Y[0][j - 1] + gap_ext
            trace_Y[0][j] = ('Y', 0, 1)

    # Fill with progress
    for i in tqdm(range(1, n + 1), desc="Filling DP matrix", unit="row"):
        ai = seq1[i - 1]
        for j in range(1, m + 1):
            bj = seq2[j - 1]

            cand_M = [
                (M[i - 1][j - 1] + sub_score(ai, bj, match, mismatch), ('M', 1, 1)),
                (X[i - 1][j - 1] + sub_score(ai, bj, match, mismatch), ('X', 1, 1)),
                (Y[i - 1][j - 1] + sub_score(ai, bj, match, mismatch), ('Y', 1, 1)),
            ]
            max_cand_M = max(cand_M, key=lambda x: x[0])
            M[i][j], trace_M[i][j] = max_cand_M

            cand_X = [
                (M[i - 1][j] + gap_open + gap_ext, ('M', 1, 0)),  # open
                (X[i - 1][j] + gap_ext,           ('X', 1, 0)),  # extend
            ]
            max_cand_X = max(cand_X, key=lambda x: x[0])
            X[i][j], trace_X[i][j] = max_cand_X

            cand_Y = [
                (M[i][j - 1] + gap_open + gap_ext, ('M', 0, 1)),  # open
                (Y[i][j - 1] + gap_ext,           ('Y', 0, 1)),  # extend
            ]
            max_cand_Y = max(cand_Y, key=lambda x: x[0])
            Y[i][j], trace_Y[i][j] = max_cand_Y

            # print(f"({i}, {j})\n\t{max_cand_M}\n\t\t{cand_M}\n\t{max_cand_X}\n\t\t{cand_X}\n\t{max_cand_Y}\n\t\t{cand_Y}", file=sys.stderr)

    # End state
    end_scores = [(M[n][m], 'M'), (X[n][m], 'X'), (Y[n][m], 'Y')]
    _, state = max(end_scores, key=lambda x: x[0])
    # print(f"({n}, {m})\n\t{state}\n\t\t{end_scores}", file=sys.stderr)

    # Traceback
    i, j = n, m
    aln1, aln2 = [], []
    # for i, r in enumerate(trace_M):
    #     print(f"M{i}: {r}", file=sys.stderr)
    # for i, r in enumerate(trace_X):
    #     print(f"X{i}: {r}", file=sys.stderr)
    # for i, r in enumerate(trace_Y):
    #     print(f"Y{i}: {r}", file=sys.stderr)
    while i > 0 or j > 0:
        if state == 'M':
            prev_state, di, dj = trace_M[i][j]
            aln1.append(seq1[i - 1]); aln2.append(seq2[j - 1])
        elif state == 'X':
            prev_state, di, dj = trace_X[i][j]
            aln1.append(seq1[i - 1]); aln2.append('-')
        elif state == 'Y':
            prev_state, di, dj = trace_Y[i][j]
            aln1.append('-');          aln2.append(seq2[j - 1])
        else:
            raise RuntimeError("Invalid state in traceback.")
        # print(f"({i}, {j}) ({prev_state}, {di}, {dj})", file=sys.stderr)
        i -= di; j -= dj; state = prev_state

    aln1.reverse(); aln2.reverse()
    return aln1, aln2


# ----------------------------
# Post-processing: gap & mismatch runs + coordinates
# ----------------------------
def _span_from_indices(starts: List[int], ends: List[int],
                       start_idx_1based: int, end_idx_1based: int) -> Tuple[Optional[int], Optional[int]]:
    """
    Compute coordinate span (min start, max end) for indices in [start_idx_1based, end_idx_1based] inclusive.
    If the range is empty (end < start), return (None, None).
    """
    if end_idx_1based < start_idx_1based:
        return None, None
    i0 = start_idx_1based - 1
    i1 = end_idx_1based - 1
    seg_starts = starts[i0:i1+1]
    seg_ends   = ends[i0:i1+1]
    return (min(seg_starts), max(seg_ends)) if seg_starts and seg_ends else (None, None)


def extract_gap_runs(
    alnA: List[str],
    alnB: List[str],
    name_a: str,
    name_b: str,
    starts_a: List[int],
    ends_a: List[int],
    starts_b: List[int],
    ends_b: List[int],
) -> pd.DataFrame:
    """
    Unified table of contiguous events:
      - deletion_in_2 : gaps in B (tokens come from A)
      - insertion_in_2: gaps in A (tokens come from B)
      - mismatch      : both present but unequal (tokens from both)

    Adds coordinate spans:
      a_start/a_end : min..max over covered A-rows
      b_start/b_end : min..max over covered B-rows
    """
    assert len(alnA) == len(alnB), "Aligned sequences must have same length"
    L = len(alnA)

    posA = 0  # original (ungapped) counters, 1-based when incremented
    posB = 0

    rows: List[Dict] = []
    k = 0
    gap_id = 1

    while k < L:
        a = alnA[k]
        b = alnB[k]

        # ---- GAP IN B (deletion_in_2) ----
        if a != '-' and b == '-':
            aln_start = k + 1
            posA_left, posB_left = posA, posB
            tokA = []
            run_len = 0
            while k < L and alnA[k] != '-' and alnB[k] == '-':
                tokA.append(alnA[k]); posA += 1; run_len += 1; k += 1
            aln_end = aln_start + run_len - 1
            # covered A indices: (posA_left+1 .. posA_right)
            a_l = posA_left + 1
            a_r = posA
            a_start, a_end = _span_from_indices(starts_a, ends_a, a_l, a_r)
            # no B bases consumed during the run
            b_start, b_end = None, None

            rows.append({
                "gap_id": gap_id,
                "file_a": name_a,
                "file_b": name_b,
                "type": "deletion_in_2",
                "aln_start": aln_start,
                "aln_end": aln_end,
                "length": run_len,
                "seq_from_file": name_a,
                "tokens": " ".join(tokA),
                "tokens_a": " ".join(tokA),
                "tokens_b": "",
                "posA_left": posA_left,
                "posB_left": posB_left,
                "posA_right": posA,
                "posB_right": posB,
                "a_start": a_start, "a_end": a_end,
                "b_start": b_start, "b_end": b_end,
            })
            gap_id += 1
            continue

        # ---- GAP IN A (insertion_in_2) ----
        if a == '-' and b != '-':
            aln_start = k + 1
            posA_left, posB_left = posA, posB
            tokB = []
            run_len = 0
            while k < L and alnA[k] == '-' and alnB[k] != '-':
                tokB.append(alnB[k]); posB += 1; run_len += 1; k += 1
            aln_end = aln_start + run_len - 1
            # covered B indices: (posB_left+1 .. posB_right)
            b_l = posB_left + 1
            b_r = posB
            b_start, b_end = _span_from_indices(starts_b, ends_b, b_l, b_r)
            a_start, a_end = None, None

            rows.append({
                "gap_id": gap_id,
                "file_a": name_a,
                "file_b": name_b,
                "type": "insertion_in_2",
                "aln_start": aln_start,
                "aln_end": aln_end,
                "length": run_len,
                "seq_from_file": name_b,
                "tokens": " ".join(tokB),
                "tokens_a": "",
                "tokens_b": " ".join(tokB),
                "posA_left": posA_left,
                "posB_left": posB_left,
                "posA_right": posA,
                "posB_right": posB,
                "a_start": a_start, "a_end": a_end,
                "b_start": b_start, "b_end": b_end,
            })
            gap_id += 1
            continue

        # ---- MISMATCH RUN ----
        if a != '-' and b != '-' and a != b:
            aln_start = k + 1
            posA_left, posB_left = posA, posB
            tokA, tokB = [], []
            run_len = 0
            while k < L and alnA[k] != '-' and alnB[k] != '-' and alnA[k] != alnB[k]:
                tokA.append(alnA[k]); tokB.append(alnB[k])
                posA += 1; posB += 1; run_len += 1; k += 1
            aln_end = aln_start + run_len - 1
            # both sides consumed
            a_l = posA_left + 1; a_r = posA
            b_l = posB_left + 1; b_r = posB
            a_start, a_end = _span_from_indices(starts_a, ends_a, a_l, a_r)
            b_start, b_end = _span_from_indices(starts_b, ends_b, b_l, b_r)

            rows.append({
                "gap_id": gap_id,
                "file_a": name_a,
                "file_b": name_b,
                "type": "mismatch",
                "aln_start": aln_start,
                "aln_end": aln_end,
                "length": run_len,
                "seq_from_file": "both",
                "tokens": "",
                "tokens_a": " ".join(tokA),
                "tokens_b": " ".join(tokB),
                "posA_left": posA_left,
                "posB_left": posB_left,
                "posA_right": posA,
                "posB_right": posB,
                "a_start": a_start, "a_end": a_end,
                "b_start": b_start, "b_end": b_end,
            })
            gap_id += 1
            continue

        # ---- MATCH ----
        if a != '-' and b != '-':  # equal
            posA += 1; posB += 1; k += 1
            continue

        # safety
        k += 1

    cols = [
        "gap_id", "file_a", "file_b", "type",
        "aln_start", "aln_end", "length",
        "seq_from_file", "tokens", "tokens_a", "tokens_b",
        "posA_left", "posB_left", "posA_right", "posB_right",
        "a_start", "a_end", "b_start", "b_end",
    ]
    return pd.DataFrame(rows, columns=cols)


# ----------------------------
# Helpers
# ----------------------------
def alignment_to_dataframe(
    aln1: List[str],
    aln2: List[str],
    col_a: str,
    col_b: str,
) -> pd.DataFrame:
    return pd.DataFrame({col_a: aln1, col_b: aln2})


# ----------------------------
# CLI
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Global alignment of annotations from two BED files."
    )
    p.add_argument("bed_a", help="First BED file.")
    p.add_argument("bed_b", help="Second BED file.")
    p.add_argument("--col", type=int, default=4,
                   help="1-based column index with labels (default: 4).")

    # outputs with defaults
    p.add_argument("--out-align", default="output_alignment.tsv",
                   help="Aligned labels TSV (default: output_alignment.tsv).")
    p.add_argument("--out-gaps", default="output_gaps.tsv",
                   help="Gap/mismatch runs TSV (default: output_gaps.tsv).")

    # scoring
    p.add_argument("--match", type=float, default=2.0, help="Match score (default: 2).")
    p.add_argument("--mismatch", type=float, default=-1.0, help="Mismatch score (default: -1).")
    p.add_argument("--gap-open", type=float, default=-4.0, help="Gap-open penalty (default: -4).")
    p.add_argument("--gap-ext", type=float, default=-1.0, help="Gap-extension penalty (default: -1).")
    return p.parse_args()


def main():
    args = parse_args()
    col_index0 = args.col - 1
    if col_index0 < 0:
        raise SystemExit("--col must be >= 1")

    # Read inputs
    labels_a, starts_a, ends_a = read_bed(args.bed_a, col_index=col_index0)
    labels_b, starts_b, ends_b = read_bed(args.bed_b, col_index=col_index0)

    # Align
    alnA, alnB = needleman_wunsch_affine(
        labels_a, labels_b,
        match=args.match, mismatch=args.mismatch,
        gap_open=args.gap_open, gap_ext=args.gap_ext,
    )

    # Write alignment
    name_a = Path(args.bed_a).name
    name_b = Path(args.bed_b).name
    align_df = alignment_to_dataframe(alnA, alnB, name_a, name_b)
    align_df.to_csv(args.out_align, sep="\t", index=True)

    # Extract runs (with coordinate spans) and write
    gaps_df = extract_gap_runs(
        alnA, alnB, name_a, name_b,
        starts_a, ends_a, starts_b, ends_b,
    )
    gaps_df.to_csv(args.out_gaps, sep="\t", index=False)


if __name__ == "__main__":
    main()

