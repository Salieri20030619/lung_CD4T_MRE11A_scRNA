#!/bin/bash
set -eu
O=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/independent_screen_20261016_v1
mkdir -p "$O"/time "$O"/provenance
trap 'r=$?; echo "$r" > "$O/exit_status"' EXIT
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python "$O/code/screen.py" "$O"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla "$O/code/time.R" "$O"
touch "$O/COMPLETED"
