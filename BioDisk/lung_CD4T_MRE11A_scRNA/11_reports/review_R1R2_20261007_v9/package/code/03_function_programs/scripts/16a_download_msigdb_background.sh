#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
CACHE_DIR=/mnt/BioDisk/wangmingyuan/.cache/R/msigdbr

URL='https://zenodo.org/records/18968178/files/msigdb.2026.1.zip?download=1'
EXPECTED_MD5='512ba99c6827141a9d471972b812d4ac'
EXPECTED_SIZE=44335790

ZIP_PATH="${CACHE_DIR}/msigdb.2026.1.zip"
PART_PATH="${ZIP_PATH}.part"
OLD_PART="${ZIP_PATH}.curltmp"
SUMMARY_PATH="${CACHE_DIR}/msigdb.2026.1.summary.rds"

mkdir -p "${CACHE_DIR}"

echo "===== MSIGDB BACKGROUND DOWNLOAD ====="
date
echo "URL=${URL}"
echo "ZIP_PATH=${ZIP_PATH}"
echo "EXPECTED_SIZE=${EXPECTED_SIZE}"
echo "EXPECTED_MD5=${EXPECTED_MD5}"

# 接管此前 R curl 留下的临时文件，以便断点续传
if [[ -s "${OLD_PART}" && ! -e "${PART_PATH}" ]]; then
    mv "${OLD_PART}" "${PART_PATH}"
    echo "[RECOVERED] ${OLD_PART} -> ${PART_PATH}"
fi

# 如果最终 ZIP 已存在，先检查是否已经完整
download_required=1

if [[ -s "${ZIP_PATH}" ]]; then
    current_md5=$(md5sum "${ZIP_PATH}" | awk '{print $1}')

    if [[ "${current_md5}" == "${EXPECTED_MD5}" ]]; then
        echo "[EXISTS_AND_VALID] ${ZIP_PATH}"
        download_required=0
    else
        if [[ ! -e "${PART_PATH}" ]]; then
            mv "${ZIP_PATH}" "${PART_PATH}"
            echo "[RESUME_EXISTING_ZIP] ${ZIP_PATH} -> ${PART_PATH}"
        else
            backup="${ZIP_PATH}.invalid.$(date +%Y%m%d_%H%M%S)"
            mv "${ZIP_PATH}" "${backup}"
            echo "[BACKUP_INVALID_ZIP] ${backup}"
        fi
    fi
fi

if [[ "${download_required}" -eq 1 ]]; then
    completed=0

    for attempt in $(seq 1 100); do
        current_size=0
        [[ -f "${PART_PATH}" ]] && current_size=$(stat -c '%s' "${PART_PATH}")

        echo "===== DOWNLOAD ATTEMPT ${attempt}/100 ====="
        echo "[RESUME_FROM] ${current_size}/${EXPECTED_SIZE} bytes"
        date

        curl \
            --fail \
            --location \
            --continue-at - \
            --connect-timeout 30 \
            --max-time 1800 \
            --speed-limit 1024 \
            --speed-time 300 \
            --retry 5 \
            --retry-delay 15 \
            --retry-all-errors \
            --silent \
            --show-error \
            --output "${PART_PATH}" \
            "${URL}" &

        curl_pid=$!

        # 每30秒向日志打印进度，避免长时间看不到输出
        while kill -0 "${curl_pid}" 2>/dev/null; do
            sleep 30
            current_size=0
            [[ -f "${PART_PATH}" ]] && current_size=$(stat -c '%s' "${PART_PATH}")
            percent=$(awk -v a="${current_size}" -v b="${EXPECTED_SIZE}" \
                'BEGIN {if (b > 0) printf "%.2f", 100*a/b; else print "0.00"}')
            echo "[PROGRESS] ${current_size}/${EXPECTED_SIZE} bytes (${percent}%)"
        done

        if wait "${curl_pid}"; then
            current_size=$(stat -c '%s' "${PART_PATH}")
            echo "[TRANSFER_FINISHED] size=${current_size}"

            current_md5=$(md5sum "${PART_PATH}" | awk '{print $1}')
            echo "[MD5] observed=${current_md5}"

            if [[ "${current_md5}" == "${EXPECTED_MD5}" ]]; then
                mv "${PART_PATH}" "${ZIP_PATH}"
                completed=1
                echo "[DOWNLOAD_PASS] ${ZIP_PATH}"
                break
            fi

            echo "[MD5_MISMATCH] 文件保留在 ${PART_PATH}"
            exit 1
        else
            rc=$?
            current_size=0
            [[ -f "${PART_PATH}" ]] && current_size=$(stat -c '%s' "${PART_PATH}")
            echo "[INTERRUPTED] rc=${rc}; retained_bytes=${current_size}"
            echo "[RETRY_IN_20_SECONDS]"
            sleep 20
        fi
    done

    if [[ "${completed}" -ne 1 ]]; then
        echo "[FAILED] Download did not complete after 100 attempts"
        exit 1
    fi
fi

echo "===== ZIP VALIDATION ====="
unzip -tq "${ZIP_PATH}"

echo "===== EXTRACTION ====="
unzip -oq "${ZIP_PATH}" -d "${CACHE_DIR}"

if [[ ! -s "${SUMMARY_PATH}" ]]; then
    echo "[ERROR] Expected file not found: ${SUMMARY_PATH}"
    exit 1
fi

echo "[CACHE_PASS] ${SUMMARY_PATH}"

echo "===== MSIGDB CACHE TEST ====="
conda run -n R_4.5.2 Rscript -e '
suppressPackageStartupMessages(library(msigdbr))
x <- msigdbr(db_species = "HS", collection = "H")
cat("Hallmark rows:", nrow(x), "\n")
cat("Hallmark gene sets:", length(unique(x$gs_name)), "\n")
stopifnot(nrow(x) > 0)
'

echo "===== START PATHWAY ANALYSIS ====="
conda run -n R_4.5.2 Rscript \
    "${PROJECT_ROOT}/scripts/16_GSE288887_CD4_pathway_analysis.R"

echo "===== FINAL STATUS ====="
echo "STATUS=PASS"
date
