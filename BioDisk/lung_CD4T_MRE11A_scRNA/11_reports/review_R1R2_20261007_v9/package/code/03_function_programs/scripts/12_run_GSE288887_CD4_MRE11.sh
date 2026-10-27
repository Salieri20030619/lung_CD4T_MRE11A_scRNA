#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
LOG_FILE="${PROJECT_ROOT}/logs/12_GSE288887_CD4_MRE11.log"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate R_4.5.2

Rscript \
  "${PROJECT_ROOT}/scripts/12_GSE288887_CD4_MRE11_analysis.R" \
  > "${LOG_FILE}" 2>&1
