# rs-annotaligner
[![CI](https://github.com/koisland/rs-annotaligner/actions/workflows/ci.yaml/badge.svg)](https://github.com/koisland/rs-annotaligner/actions/workflows/ci.yaml)

Rust port of [`annotaligner`](https://github.com/fedorrik/annotaligner) with additional features.

## Why?
* Learning exercise. Never implemented basic alignment algorithms myself.
* Outputs of `annotaligner` could be better (no coordinates). Ideally should be [`BEDPE`](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) or [`PAF`](https://github.com/lh3/miniasm/blob/master/PAF.md).
* Rust for performance and correctness.

## Usage
Clone repo.
```bash
# --recursive is optional for comparison against annotaligner submodule
git clone https://github.com/koisland/rs-annotaligner --recursive
```

Compile with Rust.
```bash
cargo build --release
```

### Alignment
Global or local alignment with affine gap penalties.
* Global alignment only outputs best alignment.
* Local alignment supports multiple alignments.

#### BEDPE
```bash
./target/release/rs-annotaligner global \
    -t test/data/input/target.bed \
    -q test/data/input/query.bed \
    -y bedpe
# ./target/release/rs-annotaligner local
```

<table>
    <tr>
        <td>Target</td>
        <td>Query</td>
        <td>Output (Global)</td>
        <td>Output (Local)</td>
    </tr>
    <tr>
<td>

```
chr1	1	2	L1
chr1	2	3	G1
chr1	3	4	P1
chr1	4	5	S1
chr1	5	6	S1
chr1	7	8	K1
chr1	8	9	Q1
chr1	9	10	T1
chr1	10	11	G1
chr1	11	12	K1
chr1	12	13	G1
chr1	13	14	S1

```

</td>
<td>

```
chr1	1	2	L1
chr1	2	3	N1
chr1	4	5	I1
chr1	5	6	T1
chr1	7	8	K1
chr1	8	9	S1
chr1	9	10	A1
chr1	10	11	G1
chr1	11	12	K1
chr1	12	13	G1
chr1	13	14	A1
```

</td>
<td>

```
#tchrom	tst	tend	qchrom	qst	qend	name	score	op	n_aln
chr1	1	2	chr1	1	2	L1~L1	-1	Match	1
chr1	2	3	.	.	.	G1~.	-1	Deletion	1
chr1	3	4	chr1	2	3	P1~N1	-1	Mismatch	1
chr1	4	5	chr1	4	5	S1~I1	-1	Mismatch	1
chr1	5	6	chr1	5	6	S1~T1	-1	Mismatch	1
chr1	7	8	chr1	7	8	K1~K1	-1	Match	1
chr1	8	9	chr1	8	9	Q1~S1	-1	Mismatch	1
chr1	9	10	chr1	9	10	T1~A1	-1	Mismatch	1
chr1	10	11	chr1	10	11	G1~G1	-1	Match	1
chr1	11	12	chr1	11	12	K1~K1	-1	Match	1
chr1	12	13	chr1	12	13	G1~G1	-1	Match	1
chr1	13	14	chr1	13	14	S1~A1	-1	Mismatch	1
```

</td>

<td>

```
#tchrom	tst	tend	qchrom	qst	qend	name	score	op	n_aln
chr1	10	11	chr1	10	11	G1~G1	6	Match	1
chr1	11	12	chr1	11	12	K1~K1	6	Match	1
chr1	12	13	chr1	12	13	G1~G1	6	Match	1
...
```

</td>

</table>

#### PAF
[Pairwise mApping Format](https://github.com/lh3/miniasm/blob/master/PAF.md#paf-a-pairwise-mapping-format) adds more alignment metadata.

The following tags are added:
* `cg` - CIGAR string
* `nr` - Number of aligned annotations (custom)
* `AS` - DP score
* `NM` - Number of mismatches/gaps
* `de` - Gap-compressed identity based on interval lengths.
* `tg` - Target annotation gap percentage. If aligned over annotation gaps, will be >0.0% (custom)
* `qg` - Query annotation gap percentage. If aligned over annotation gaps, will be >0.0% (custom)

MAPQ is calculated based on averaging `de` and `AS` relative to the best, max-scoring alignment.
* Minimum of `0` and maximum of `60`
* Formula: `30 * log_10(100.0 * (0.5 * (identity + (dp_score / top_dp_score))))`

Column 2 and 7 are just the aligned length and the same as `nr`.

```bash
./target/release/rs-annotaligner global \
    -t test/data/input/target.bed \
    -q test/data/input/query.bed \
    -y paf
# ./target/release/rs-annotaligner global
```

<table>
    <tr>
        <td>Target</td>
        <td>Query</td>
        <td>Output (Global)</td>
        <td>Output (Local)</td>
    </tr>
    <tr>
<td>

```
chr1	1	2	And
chr1	2	3	Bow
chr1	3	4	Cow
chr1	4	5	Anger
chr1	5	6	Anger
chr1	10	11	Creek
chr1	11	13	Bag
chr1	13	15	Banana
chr1	15	20	Aww
```

</td>
<td>

```
chr1	1	2	And
chr1	2	3	Bow
chr1	3	4	Cow
chr1	4	5	Zed
chr1	5	6	Zed
chr1	6	7	Creek
chr1	7	8	Bag
chr1	8	9	Bag
```

</td>
<td>

```
chr1    8       1       9       +       chr1    19      1       20      6       15      55      cg:Z:3=2X3=2D5X nr:i:9  AS:i:2  NM:i:9  de:f:0.42857143 tg:f:0.21052632 qg:f:0
```

</td>

<td>

```
chr1    7       1       8       +       chr1    12      1       13      6       8       58      cg:Z:3=2X3=     nr:i:7  AS:i:8  NM:i:2  de:f:0.75       tg:f:0.33333334 qg:f:0
chr1    3       1       4       +       chr1    3       1       4       3       3       58      cg:Z:3= nr:i:3  AS:i:6  NM:i:0  de:f:1  tg:f:0  qg:f:0
chr1    8       1       9       +       chr1    14      1       15      6       10      56      cg:Z:3=2X3=2X   nr:i:8  AS:i:7  NM:i:4  de:f:0.6        tg:f:0.2857143  qg:f:0
chr1    6       1       7       +       chr1    10      1       11      4       6       55      cg:Z:3=2X1=     nr:i:6  AS:i:6  NM:i:2  de:f:0.6666667  tg:f:0.4        qg:f:0
chr1    4       1       5       +       chr1    4       1       5       3       4       55      cg:Z:3=1X       nr:i:4  AS:i:5  NM:i:1  de:f:0.75       tg:f:0  qg:f:0
...
```

</td>

</table>

## Test

### Suite
Test suite and all examples.
```bash
cargo test --release
```

### Against annotaligner with HG008 example
Chromosomal fusion of chr6 and chr7 in HG008-T.
```bash
cargo run --release -- \
    global \
    -t test/data/input/HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz \
    -q test/data/input/HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz
```

```bash
# pip install -r test/requirements.txt
python test/annotaligner/annotaligner.py \
    test/data/input/HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz \
    test/data/input/HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz
```

## TO-DO
* [x] Gzip input
* [x] Smith-Waterman
* [x] Output as PAF
* [x] Output as BEDPE
* [ ] Benchmarks
