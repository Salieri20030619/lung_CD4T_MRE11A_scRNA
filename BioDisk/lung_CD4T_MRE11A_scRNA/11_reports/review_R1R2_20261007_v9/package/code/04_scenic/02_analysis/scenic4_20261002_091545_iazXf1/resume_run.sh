#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
OUT="$ROOT/02_analysis/scenic4_20261002_091545_iazXf1"
PREV1="$ROOT/02_analysis/mechanism_step1_20260927_152218_s4V5eF"
SCENIC1="$ROOT/02_analysis/scenic1_20260929_103513_xRoEnr"
PY="$ROOT/software/pyscenic_0.12.1_py310/bin/python"
AUDIT="$ROOT/tmp/scenic4_repair_20261002_1205"
export ROOT OUT PREV1 SCENIC1
export WORKERS=4 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR="$AUDIT/matplotlib"
exec 9>"$ROOT/tmp/lung_scenic4.lock"
flock -n 9 || { echo SCENIC4_ALREADY_RUNNING; exit 1; }
trap 'rc=$?; printf "RESUME_EXIT_CODE=%s TIME=%s\n" "$rc" "$(date -Is)"; if ((rc != 0)); then printf "LUNG_SCENIC4_RESUME_FAIL OUTPUT=%s\n" "$OUT"; fi' EXIT
printf 'RESUME_STARTED=%s OUTPUT=%s\n' "$(date -Is)" "$OUT"
sha256sum -c "$AUDIT/original_sha256.txt"
"$PY" -B -u - <<'PREFLIGHT'
import os,runpy
from pathlib import Path
runpy.run_path(str(Path(os.environ['OUT'])/'validate.py'),run_name='dependency_preflight')
print('FULL_RUNTIME_IMPORT_PREFLIGHT=PASS',flush=True)
PREFLIGHT
printf 'REUSING_FROZEN_INPUTS=YES ORIGINAL_VALIDATE_UNCHANGED=YES\n'
"$PY" -B -u "$OUT/validate.py"
