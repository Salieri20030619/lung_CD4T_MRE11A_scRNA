#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
RAW288="${ROOT}/01_raw_data/GSE288887"
EXT288="${RAW288}/extracted"
REUSED288="${RAW288}/reused_filtered"
ARCHIVE="${RAW288}/GSE288887_RAW.tar"
REPAIR_LOG="${ROOT}/logs/02_h5_repair_all.log"
AUDIT_LOG="${ROOT}/logs/02_corrected_filtered_qc_audit.log"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288887/suppl"
STAMP="$(date +%Y%m%d_%H%M%S)"

mkdir -p "${EXT288}" "${REUSED288}"

exec > >(tee "${REPAIR_LOG}") 2>&1

NEW_FILES=(
"GSM8778048_4859_2_filtered_feature_bc_matrix.h5"
"GSM8778049_6281_1_filtered_feature_bc_matrix.h5"
"GSM8778050_6281_2_filtered_feature_bc_matrix.h5"
"GSM8778051_7539_filtered_feature_bc_matrix.h5"
"GSM8778052_8981_1_filtered_feature_bc_matrix.h5"
"GSM8778053_8981_2_filtered_feature_bc_matrix.h5"
"GSM8778054_9367_filtered_feature_bc_matrix.h5"
"GSM8778055_9564_filtered_feature_bc_matrix.h5"
"GSM8808454_5465_1_filtered_feature_bc_matrix.h5"
"GSM8808455_5465_2_filtered_feature_bc_matrix.h5"
)

REUSED_FILES=(
"GSE288887_4546_1_filtered_feature_bc_matrix.h5"
"GSE288887_4546_2_filtered_feature_bc_matrix.h5"
"GSE288887_4736_1_filtered_feature_bc_matrix.h5"
"GSE288887_4736_2_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD092_1_LU_Whole_C2_X5SCR_F02526_HLKH5DSXY_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD094_1_LU_Whole_C2_X5SCR_F02524_HLKH5DSXY_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD098_1_LU_Whole_C2_X5SCR_F02528_HLKH5DSXY_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD103_1_LU_Whole_C1_X5SCR_F02617_HLKH5DSXY_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD105_1_LU_Whole_C1_X5SCR_F02609_HLKH5DSXY_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD106_1_LU_Whole_C1_X5SCR_F02611_HLKH5DSXY_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD108_1_LU_Whole_C1_X5SC2_F03548_H3HLJDSX2-H5W7FDSX2_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD111_1_LU_Whole_C1_X5SC2_F03552_H3HLJDSX2-H5W7FDSX2_filtered_feature_bc_matrix.h5"
"GSE288887_IPF_VUHD117_1_LU_MNC_C1_X5SC2_F03708_H3HLJDSX2-H5W7FDSX2_filtered_feature_bc_matrix.h5"
"GSE288887_vuhd068_filtered_feature_bc_matrix.h5"
"GSE288887_vuhd070_filtered_feature_bc_matrix.h5"
)

validate_h5() {
    local path="$1"

    [[ -s "${path}" ]] || return 1

    Rscript - "${path}" <<'RSCRIPT' >/dev/null 2>&1
suppressPackageStartupMessages(library(hdf5r))

args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]

handle <- H5File$new(path, mode = "r")
on.exit(handle$close_all(), add = TRUE)

if (!"matrix" %in% names(handle)) {
    stop("Missing /matrix group")
}

matrix_group <- handle[["matrix"]]

required <- c(
    "data",
    "indices",
    "indptr",
    "shape",
    "barcodes",
    "features"
)

missing <- setdiff(required, names(matrix_group))

if (length(missing) > 0) {
    stop(
        paste(
            "Missing matrix objects:",
            paste(missing, collapse = ",")
        )
    )
}

shape <- matrix_group[["shape"]][]

if (length(shape) != 2 || any(shape <= 0)) {
    stop("Invalid matrix shape")
}
RSCRIPT
}

repair_new_file() {
    local filename="$1"
    local target="${EXT288}/${filename}"
    local temporary="${target}.downloading_${STAMP}"
    local archive_entry

    if validate_h5 "${target}"; then
        echo -e "PASS\tNEW\t${filename}\t$(stat -c '%s' "${target}")"
        return 0
    fi

    if [[ -e "${target}" ]]; then
        mv "${target}" "${target}.corrupt_${STAMP}"
        echo "BACKUP_NEW=${target}.corrupt_${STAMP}"
    else
        echo "MISSING_NEW=${target}"
    fi

    archive_entry="$(
        tar -tf "${ARCHIVE}" |
        awk -v wanted="${filename}" '
        {
            n = split($0, item, "/")
            if (item[n] == wanted) {
                print $0
                exit
            }
        }'
    )"

    if [[ -z "${archive_entry}" ]]; then
        echo "ERROR: archive entry not found for ${filename}"
        return 1
    fi

    tar -xOf "${ARCHIVE}" "${archive_entry}" > "${temporary}"

    if ! validate_h5 "${temporary}"; then
        echo "ERROR: extracted H5 remains invalid: ${filename}"
        return 1
    fi

    mv "${temporary}" "${target}"

    echo -e "REPAIRED\tNEW\t${filename}\t$(stat -c '%s' "${target}")"
}

repair_reused_file() {
    local filename="$1"
    local target="${REUSED288}/${filename}"
    local temporary="${target}.downloading_${STAMP}"
    local url="${BASE_URL}/${filename}"

    if validate_h5 "${target}"; then
        echo -e "PASS\tREUSED\t${filename}\t$(stat -c '%s' "${target}")"
        return 0
    fi

    if [[ -e "${target}" ]]; then
        mv "${target}" "${target}.corrupt_${STAMP}"
        echo "BACKUP_REUSED=${target}.corrupt_${STAMP}"
    else
        echo "MISSING_REUSED=${target}"
    fi

    if command -v wget >/dev/null 2>&1; then
        wget \
            --tries=30 \
            --timeout=60 \
            --retry-connrefused \
            --waitretry=5 \
            --output-document="${temporary}" \
            "${url}"
    else
        curl \
            --location \
            --fail \
            --retry 30 \
            --retry-delay 5 \
            --connect-timeout 60 \
            --output "${temporary}" \
            "${url}"
    fi

    if ! validate_h5 "${temporary}"; then
        echo "ERROR: downloaded H5 remains invalid: ${filename}"
        return 1
    fi

    mv "${temporary}" "${target}"

    echo -e "REPAIRED\tREUSED\t${filename}\t$(stat -c '%s' "${target}")"
}

echo "===== H5 REPAIR START ====="
date
hostname

echo "===== CHECK 10 NEW FILTERED MATRICES ====="

for filename in "${NEW_FILES[@]}"; do
    repair_new_file "${filename}"
done

echo "===== CHECK 15 REUSED FILTERED MATRICES ====="

for filename in "${REUSED_FILES[@]}"; do
    repair_reused_file "${filename}"
done

echo "===== FINAL FULL VALIDATION ====="

PASS_N=0
FAIL_N=0

for filename in "${NEW_FILES[@]}"; do
    target="${EXT288}/${filename}"

    if validate_h5 "${target}"; then
        echo -e "FINAL_PASS\t${target}\t$(stat -c '%s' "${target}")"
        PASS_N=$((PASS_N + 1))
    else
        echo -e "FINAL_FAIL\t${target}"
        FAIL_N=$((FAIL_N + 1))
    fi
done

for filename in "${REUSED_FILES[@]}"; do
    target="${REUSED288}/${filename}"

    if validate_h5 "${target}"; then
        echo -e "FINAL_PASS\t${target}\t$(stat -c '%s' "${target}")"
        PASS_N=$((PASS_N + 1))
    else
        echo -e "FINAL_FAIL\t${target}"
        FAIL_N=$((FAIL_N + 1))
    fi
done

echo "validated_pass=${PASS_N}"
echo "validated_fail=${FAIL_N}"

if [[ "${PASS_N}" -ne 25 || "${FAIL_N}" -ne 0 ]]; then
    echo "STATUS=H5_VALIDATION_FAILED"
    exit 1
fi

echo "===== RESTART CORRECTED AUDIT ====="

nohup bash "${ROOT}/scripts/02_corrected_filtered_qc_audit.sh" \
    >/dev/null 2>&1 &

AUDIT_PID=$!

echo "${AUDIT_PID}" \
    > "${ROOT}/logs/02_corrected_filtered_qc_audit.pid"

echo "AUDIT_PID=${AUDIT_PID}"
echo "AUDIT_LOG=${AUDIT_LOG}"
echo "STATUS=REPAIRED_AND_RESTARTED"
date
