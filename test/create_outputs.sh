#!/bin/bash

set -euo pipefail

cargo run --release -- local \
    -t test/data/input/target_local.bed \
    -q test/data/input/query_local.bed \
    -s 6  > test/data/output/basic_example_local.bedpe
cargo run --release -- local \
    -t test/data/input/target_local.bed \
    -q test/data/input/query_local.bed \
    -y paf \
    -s 6  > test/data/output/basic_example_local.paf
cargo run --release -- global \
    -t test/data/input/target_local.bed \
    -q test/data/input/query_local.bed > test/data/output/basic_example_global.bedpe
cargo run --release -- global \
    -t test/data/input/HG008-N_v6.3_chr7_hap2_57312660-64850688_stv.bed.gz \
    -q test/data/input/HG008-T_v3.2_chr6_chr7_chr11_hap2_60228206-67527215_stv.bed.gz \
| gzip > test/data/output/HG008-TN_chr6_chr7_fusion.bed.gz
cargo run --release -- global \
    -t test/data/input/target.bed \
    -q test/data/input/query.bed > test/data/output/basic_example.bedpe
