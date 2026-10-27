#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
OUT="$ROOT/02_analysis/scenic1_20260929_103513_xRoEnr"
PREV1="$ROOT/02_analysis/mechanism_step1_20260927_152218_s4V5eF"
PREV2="$ROOT/02_analysis/mechanism_step2_20260927_154919_EsN4gk"
PREV7="$ROOT/02_analysis/mechanism_step7_20260929_100321_6hm4dZ"
CACHE="$ROOT/00_raw/SCENIC_resources_v10"
PY="$ROOT/software/pyscenic_0.12.1_py310/bin/python"
export ROOT OUT PREV1 PREV2 PREV7 CACHE
export WORKERS=8 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
if [[ -n "${SCENIC_PROXY:-}" ]]; then
 export http_proxy="$SCENIC_PROXY" https_proxy="$SCENIC_PROXY" HTTP_PROXY="$SCENIC_PROXY" HTTPS_PROXY="$SCENIC_PROXY"
 export no_proxy=localhost,127.0.0.1 NO_PROXY=localhost,127.0.0.1
fi
bash -n "$0"
[[ -d "$OUT" && -x "$PY" && -f "$OUT/smoke.py" && -f "$OUT/pipeline.py" ]]
exec 9>"$ROOT/tmp/lung_scenic1.lock"
flock -n 9 || { echo SCENIC1_ALREADY_RUNNING; exit 1; }
trap 'rc=$?; if ((rc != 0)); then printf "LUNG_SCENIC1_RESUME_FAIL EXIT_CODE=%s OUTPUT=%s\n" "$rc" "$OUT"; fi' EXIT
[[ ! -e "$OUT/status.json" && ! -d "$OUT/Mouse" && ! -d "$OUT/Human" ]] || { echo STOP_EXISTING_ANALYSIS_OUTPUT_REQUIRES_REVIEW; exit 1; }
"$PY" -m py_compile "$OUT/smoke.py" "$OUT/pipeline.py"
echo SYNTAX_CHECK=PASS
cp -- "$0" "$OUT/resume_environment.sh"
printf 'setuptools==80.9.0\n' > "$OUT/environment_repair_constraints.txt"
"$PY" -m pip install --no-deps --timeout 120 --retries 8 'setuptools==80.9.0'
"$PY" -m pip check
"$PY" -c 'import pkg_resources, ctxcore; from pyscenic.cli.utils import load_signatures; from importlib.metadata import version; print("SETUPTOOLS="+version("setuptools")); print("PKG_RESOURCES_IMPORT_PASS")'
"$PY" -m pip freeze > "$OUT/environment_freeze.txt"
SMOKE=$(mktemp -d "$OUT/environment_smoke_resume_XXXXXX")
"$PY" -u "$OUT/smoke.py" "$SMOKE"
"$PY" -u "$OUT/pipeline.py"
