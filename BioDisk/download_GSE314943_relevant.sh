#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
DATASET_ROOT="${PROJECT_ROOT}/00_raw/GSE314943"
META_DIR="${DATASET_ROOT}/metadata"
SCRNA_DIR="${DATASET_ROOT}/scRNA"
SCATAC_DIR="${DATASET_ROOT}/scATAC"
LOG_DIR="${PROJECT_ROOT}/logs"
MINIML_TGZ="${META_DIR}/GSE314943_family.xml.tgz"
MINIML_XML="${META_DIR}/GSE314943_family.xml"
MANIFEST="${META_DIR}/GSE314943_Saline_Silica_manifest.tsv"
JOBS_NUL="${META_DIR}/GSE314943_download_jobs.nul"
AUDIT_TSV="${META_DIR}/GSE314943_download_audit.tsv"
MD5_TSV="${META_DIR}/GSE314943_download_md5.tsv"
LOG_FILE="${LOG_DIR}/download_GSE314943_relevant.log"
MIN_FREE_KB=$((15 * 1024 * 1024))

mkdir -p "${META_DIR}" "${SCRNA_DIR}" "${SCATAC_DIR}" "${LOG_DIR}"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "STATUS=START"
echo "DATASET=GSE314943"
echo "SCOPE=Saline_and_Silica_scRNA_scATAC"
echo "START_TIME=$(date '+%F %T %z')"

AVAILABLE_KB=$(df -Pk "${DATASET_ROOT}" | tail -n 1 | tr -s ' ' | cut -d ' ' -f 4)
echo "AVAILABLE_KB=${AVAILABLE_KB}"
if [[ ! "${AVAILABLE_KB}" =~ ^[0-9]+$ ]] || (( AVAILABLE_KB < MIN_FREE_KB )); then
    echo "ERROR=INSUFFICIENT_DISK_SPACE_REQUIRE_AT_LEAST_15GB_FREE"
    exit 1
fi

curl -L --fail --show-error --silent \
    --retry 10 --retry-delay 10 --retry-all-errors \
    --continue-at - \
    --output "${MINIML_TGZ}.part" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE314nnn/GSE314943/miniml/GSE314943_family.xml.tgz"
mv "${MINIML_TGZ}.part" "${MINIML_TGZ}"
gzip -t "${MINIML_TGZ}"
tar --no-same-owner -xzf "${MINIML_TGZ}" -C "${META_DIR}"
test -s "${MINIML_XML}"

conda run --no-capture-output -n batch_bench python - "${MINIML_XML}" "${MANIFEST}" "${JOBS_NUL}" <<'PY'
import csv
import os
import sys
import xml.etree.ElementTree as ET

xml_path, manifest_path, jobs_path = sys.argv[1:4]
namespace = {"m": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
root = ET.parse(xml_path).getroot()
rows = []

for sample in root.findall("m:Sample", namespace):
    title = sample.findtext("m:Title", default="", namespaces=namespace).strip()
    if not ("Saline" in title or "Silica" in title):
        continue
    if "scRNA" in title:
        modality = "scRNA"
    elif "scATAC" in title:
        modality = "scATAC"
    else:
        continue
    group = "Saline" if "Saline" in title else "Silica"
    replicate = title.split("_")[-1].split(",")[0].strip()
    accession = sample.findtext(
        "m:Accession", default="", namespaces=namespace
    ).strip()
    for element in sample.findall("m:Supplementary-Data", namespace):
        if not element.text:
            continue
        url = element.text.strip().replace("ftp://", "https://", 1)
        filename = os.path.basename(url)
        rows.append(
            [accession, title, modality, group, replicate, url, filename]
        )

rows.sort(key=lambda row: ({"scRNA": 0, "scATAC": 1}[row[2]], row[3], int(row[4]), row[6]))
expected_libraries = {
    (modality, group, replicate)
    for modality in ("scRNA", "scATAC")
    for group in ("Saline", "Silica")
    for replicate in ("1", "2", "3")
}
observed_libraries = {(row[2], row[3], row[4]) for row in rows}

if observed_libraries != expected_libraries:
    missing = sorted(expected_libraries - observed_libraries)
    extra = sorted(observed_libraries - expected_libraries)
    raise SystemExit(f"Library design mismatch; missing={missing}; extra={extra}")
if len(rows) != 32:
    raise SystemExit(f"Expected 32 supplementary files, observed {len(rows)}")

with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ["accession", "title", "modality", "group", "replicate", "url", "filename"]
    )
    writer.writerows(rows)

with open(jobs_path, "wb") as handle:
    for row in rows:
        for value in (row[2], row[5], row[6]):
            handle.write(value.encode("utf-8") + b"\0")

print(f"MANIFEST_FILES={len(rows)}")
print(f"MANIFEST_LIBRARIES={len(observed_libraries)}")
print("MANIFEST_CHECK=PASS")
PY

download_one() {
    local modality="$1"
    local url="$2"
    local filename="$3"
    local target_dir="${DATASET_ROOT}/${modality}"
    local target="${target_dir}/${filename}"
    local partial="${target}.part"

    mkdir -p "${target_dir}"
    if [[ -s "${target}" ]]; then
        echo "DOWNLOAD_SKIP=${modality}/${filename}"
        return 0
    fi

    echo "DOWNLOAD_START=${modality}/${filename}"
    curl -L --fail --show-error --silent \
        --retry 100 --retry-delay 10 --retry-all-errors --connect-timeout 30 --speed-time 120 --speed-limit 1024 --http1.1 \
        --continue-at - \
        --output "${partial}" \
        "${url}"
    test -s "${partial}"
    mv "${partial}" "${target}"
    echo "DOWNLOAD_DONE=${modality}/${filename}"
}

export -f download_one
export DATASET_ROOT

if ! xargs -0 -n 3 -P 1 bash -c 'download_one "$1" "$2" "$3"' _ < "${JOBS_NUL}"; then
    echo "ERROR=ONE_OR_MORE_DOWNLOADS_FAILED"
    exit 1
fi

conda run --no-capture-output -n batch_bench python - "${MANIFEST}" "${DATASET_ROOT}" "${AUDIT_TSV}" <<'PY'
import csv
import gzip
import os
import sys

manifest_path, dataset_root, audit_path = sys.argv[1:4]
audits = []
failures = []

with open(manifest_path, encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

for row in rows:
    path = os.path.join(dataset_root, row["modality"], row["filename"])
    status = "PASS"
    detail = "nonempty_file"
    try:
        size = os.path.getsize(path)
        if size <= 0:
            raise ValueError("empty file")
        if path.endswith(".h5"):
            with open(path, "rb") as handle:
                if handle.read(8) != b"\x89HDF\r\n\x1a\n":
                    raise ValueError("invalid HDF5 signature")
            detail = "valid_HDF5_signature"
        elif path.endswith(".mtx.gz"):
            with gzip.open(path, "rt") as handle:
                if not handle.readline().startswith("%%MatrixMarket"):
                    raise ValueError("invalid MatrixMarket header")
            detail = "valid_MatrixMarket_header"
        elif path.endswith(".gz"):
            with gzip.open(path, "rb") as handle:
                handle.read(1)
            detail = "readable_gzip_stream"
    except Exception as exc:
        size = os.path.getsize(path) if os.path.exists(path) else 0
        status = "FAIL"
        detail = str(exc)
        failures.append(path)
    audits.append([row["accession"], row["modality"], row["group"], row["filename"], size, status, detail])

with open(audit_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["accession", "modality", "group", "filename", "bytes", "status", "detail"])
    writer.writerows(audits)

print(f"AUDIT_FILES={len(audits)}")
print(f"AUDIT_FAILURES={len(failures)}")
if failures:
    for path in failures:
        print(f"AUDIT_FAILED_FILE={path}")
    raise SystemExit(1)
print("STRUCTURE_CHECK=PASS")
PY

while IFS=$'\t' read -r accession title modality group replicate url filename; do
    if [[ "${accession}" == "accession" ]]; then
        continue
    fi
    if [[ "${filename}" == *.gz ]]; then
        gzip -t "${DATASET_ROOT}/${modality}/${filename}"
    fi
done < "${MANIFEST}"
echo "GZIP_CHECK=PASS"

find "${SCRNA_DIR}" "${SCATAC_DIR}" -maxdepth 1 -type f ! -name '*.part' -print0 \
    | sort -z \
    | xargs -0 md5sum > "${MD5_TSV}"

EXPECTED_FILES=$(($(wc -l < "${MANIFEST}") - 1))
AUDITED_PASS=$(awk -F '\t' 'NR > 1 && $6 == "PASS" {count++} END {print count+0}' "${AUDIT_TSV}")
PARTIAL_FILES=$(find "${SCRNA_DIR}" "${SCATAC_DIR}" -maxdepth 1 -type f -name '*.part' | wc -l)

echo "EXPECTED_FILES=${EXPECTED_FILES}"
echo "AUDITED_PASS=${AUDITED_PASS}"
echo "PARTIAL_FILES=${PARTIAL_FILES}"

if [[ "${EXPECTED_FILES}" -ne 32 || "${AUDITED_PASS}" -ne 32 || "${PARTIAL_FILES}" -ne 0 ]]; then
    echo "ERROR=FINAL_COUNT_CHECK_FAILED"
    exit 1
fi

du -sh "${DATASET_ROOT}"
echo "END_TIME=$(date '+%F %T %z')"
echo "GSE314943_DOWNLOAD_PASS"
