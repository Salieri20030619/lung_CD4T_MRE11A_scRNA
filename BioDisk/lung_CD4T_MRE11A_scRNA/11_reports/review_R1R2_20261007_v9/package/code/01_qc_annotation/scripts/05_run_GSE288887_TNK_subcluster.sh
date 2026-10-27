#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
LOG_FILE="${PROJECT_ROOT}/logs/05_GSE288887_TNK_subcluster.log"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate R_4.5.2

Rscript "${PROJECT_ROOT}/scripts/05_GSE288887_TNK_subcluster.R" \
  > "${LOG_FILE}" 2>&1
