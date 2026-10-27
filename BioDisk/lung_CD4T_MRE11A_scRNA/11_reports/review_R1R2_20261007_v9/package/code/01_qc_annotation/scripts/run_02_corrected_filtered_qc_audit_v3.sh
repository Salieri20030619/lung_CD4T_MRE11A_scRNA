#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
LOG="${ROOT}/logs/02_corrected_filtered_qc_audit_v3.log"

exec > "${LOG}" 2>&1

echo "===== V3 AUDIT START ====="
date
hostname

Rscript \
    "${ROOT}/scripts/02_corrected_filtered_qc_audit_v3.R"

echo "===== COMPLETION ====="
echo "status=PASS_WITH_COVERAGE_WARNING"
date
