# rs-annotaligner
Rust port of [`annotaligner`](https://github.com/fedorrik/annotaligner) with additional features.

## Why?
* Learning exercise. Never implemented an aligner myself.
* Outputs of `annotaligner` could be better (no coordinates). BEDPE or paf.
* Rust for performance and correctness.

## Test
```bash
cargo run --release -- -t test/data/HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz -q test/data/HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz
```

```bash
python test/annotaligner.py test/data/HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz test/data/HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz
```

## TODO
* [x] Gzip input
* [ ] Smith-Waterman
* [ ] Output as PAF
* [x] Output as BEDPE
* [ ] Benchmarks
