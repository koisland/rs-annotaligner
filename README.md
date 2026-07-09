# rs-annotaligner
[![CI](https://github.com/koisland/rs-annotaligner/actions/workflows/ci.yaml/badge.svg)](https://github.com/koisland/rs-annotaligner/actions/workflows/ci.yaml)

Rust port of [`annotaligner`](https://github.com/fedorrik/annotaligner) with additional features.

## Why?
* Learning exercise. Never implemented basic alignment algorithms myself.
* Outputs of `annotaligner` could be better (no coordinates). Ideally shoule be [`BEDPE`](https://bedtools.readthedocs.io/en/latest/content/general-usage.html) or [`paf`](https://github.com/lh3/miniasm/blob/master/PAF.md).
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

### BEDPE
Global or local alignment with affine gap penalties.
```bash
./target/release/rs-annotaligner \
    -t test/data/input/target.bed \
    -q test/data/input/query.bed \
    -a global # Or local
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
chr1	1	2	chr1	1	2	L1~L1	Match
chr1	2	3	.	.	.	G1~.	Deletion
chr1	3	4	chr1	2	3	P1~N1	Mismatch
chr1	4	5	chr1	4	5	S1~I1	Mismatch
chr1	5	6	chr1	5	6	S1~T1	Mismatch
chr1	7	8	chr1	7	8	K1~K1	Match
chr1	8	9	chr1	8	9	Q1~S1	Mismatch
chr1	9	10	chr1	9	10	T1~A1	Mismatch
chr1	10	11	chr1	10	11	G1~G1	Match
chr1	11	12	chr1	11	12	K1~K1	Match
chr1	12	13	chr1	12	13	G1~G1	Match
chr1	13	14	chr1	13	14	S1~A1	Mismatch
```

</td>

<td>

```
chr1	10	11	chr1	10	11	G1~G1	Match
chr1	11	12	chr1	11	12	K1~K1	Match
chr1	12	13	chr1	12	13	G1~G1	Match
```

</td>

</table>

### PAF
```bash
TODO
```

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
    -t test/data/input/HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz \
    -q test/data/input/HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz
```

```bash
# pip install -r test/requirements.txt
python test/annotaligner/annotaligner.py \
    test/data/input/HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz \
    test/data/input/HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz
```

## TODO
* [x] Gzip input
* [x] Smith-Waterman
* [ ] Output as PAF
* [x] Output as BEDPE
* [ ] Benchmarks
