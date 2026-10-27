#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
LOG_FILE="${PROJECT_ROOT}/logs/10_CD4_DEG_robustness_audit.log"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate R_4.5.2

Rscript \
  "${PROJECT_ROOT}/scripts/10_CD4_DEG_robustness_audit.R" \
  > "${LOG_FILE}" 2>&1
