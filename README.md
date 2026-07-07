# rs-annotaligner
Rust port of [`annotaligner`](https://github.com/fedorrik/annotaligner) with additional features.

## Why?
* Learning exercise. Never implemented basic alignment algorithms myself.
* Outputs of `annotaligner` could be better (no coordinates). BEDPE or paf.
* Rust for performance and correctness.

## Usage
Clone repo.
```bash
# --recursive is optional for comparison against annotaligner submodule
git clone https://github.com/koisland/rs-annotaligner --recursive
```

Compile.
```bash
cargo build --release
```

### BEDPE
```bash
./target/release/rs-annotaligner \
    -t test/data/input/target.bed \
    -q test/data/input/query.bed
```

<table>
    <tr>
        <td>Target</td>
        <td>Query</td>
        <td>Output</td>
    </tr>
    <tr>
<td>

```
chr1	1	2	And
chr1	2	3	Bow
chr1	3	4	Cow
chr1	4	5	Anger
chr1	5	6	Anger
chr1	6	7	Creek
chr1	7	8	Bag
chr1	8	9	Banana
chr1	9	10	Aww
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
chr1	1	2	chr1	1	2	And~And	Match
chr1	2	3	chr1	2	3	Bow~Bow	Match
chr1	3	4	chr1	3	4	Cow~Cow	Match
chr1	4	5	chr1	4	5	Anger~Zed	Mismatch
chr1	5	6	chr1	5	6	Anger~Zed	Mismatch
chr1	6	7	chr1	6	7	Creek~Creek	Match
chr1	7	8	chr1	7	8	Bag~Bag	Match
chr1	8	9	.	.	.	Banana~.	Deletion
chr1	9	10	chr1	8	9	Aww~Bag	Mismatch
```

</td>
</table>

### PAF
```bash
TODO
```

## Test
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
* [ ] Smith-Waterman
* [ ] Output as PAF
* [x] Output as BEDPE
* [ ] Benchmarks
