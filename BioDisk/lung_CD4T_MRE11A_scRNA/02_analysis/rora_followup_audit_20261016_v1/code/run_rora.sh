#!/bin/bash
set -eu
O=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/02_analysis/rora_followup_audit_20261016_v1
trap 'r=$?; echo "$r" > "$O/exit_status"' EXIT
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python "$O/code/rora_check.py" "$O"
