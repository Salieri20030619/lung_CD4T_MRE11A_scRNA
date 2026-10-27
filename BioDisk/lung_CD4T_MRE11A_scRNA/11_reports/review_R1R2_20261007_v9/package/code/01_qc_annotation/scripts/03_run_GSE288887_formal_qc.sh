#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
LOG="${ROOT}/logs/03_GSE288887_formal_qc_doublet.log"

exec > "${LOG}" 2>&1

echo "===== GSE288887 FORMAL PREPROCESSING START ====="
date
hostname
echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNKNOWN}"

Rscript \
    "${ROOT}/scripts/03_GSE288887_formal_qc_doublet.R"

cat > "${ROOT}/00_project/PROJECT_STATUS.md" <<'STATUS'
# Project status

- GSE288887 selected as the primary dataset.
- GSE174725 retained as a candidate validation dataset only.
- GSE288887 filtered count matrices validated.
- Formal per-sample QC completed.
- Per-sample doublet detection completed using scDblFinder.
- Clean singlet cells merged into one GSE-level Seurat object.
- Clustering, annotation, CD4+ T-cell extraction and DEG analysis have not started.
STATUS

echo "===== COMPLETION ====="
echo "status=PASS"
date
