#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
LOG="${ROOT}/logs/04_GSE288887_broad_clustering.log"

exec > "${LOG}" 2>&1

echo "===== BROAD CLUSTERING START ====="
date
hostname
echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-UNKNOWN}"

Rscript \
    "${ROOT}/scripts/04a_GSE288887_broad_clustering.R"

echo "===== CLUSTERING CHECKPOINT PASS ====="

Rscript \
    "${ROOT}/scripts/04b_GSE288887_marker_lineage.R"

cat > "${ROOT}/00_project/PROJECT_STATUS.md" <<'STATUS'
# Project status

- GSE288887 selected as the primary dataset.
- Formal QC and per-sample doublet removal completed.
- Final clean singlet cells: 192,575.
- Broad normalization, PCA and Harmony correction completed.
- Multi-resolution broad clustering completed.
- Broad cluster markers and preliminary lineage scores completed.
- Preliminary lineage labels require marker-based review before CD4+ T-cell extraction.
- GSE174725 remains a candidate validation dataset.
STATUS

echo "===== COMPLETION ====="
echo "status=PASS"
date
