#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
RAW288="${ROOT}/01_raw_data/GSE288887"
EXT288="${RAW288}/extracted"
REUSED288="${RAW288}/reused_filtered"
RAW174="${ROOT}/01_raw_data/GSE174725"
EXT174="${RAW174}/extracted"
AUDIT="${ROOT}/03_data_audit/corrected_v2"
LOG="${ROOT}/logs/02_corrected_filtered_qc_audit.log"

mkdir -p \
    "${REUSED288}" \
    "${AUDIT}"

exec > "${LOG}" 2>&1

echo "===== CORRECTED AUDIT START ====="
date
hostname
echo "ROOT=${ROOT}"

download_file() {
    local url="$1"
    local output="$2"

    if [[ -s "${output}" ]]; then
        echo "[EXISTS] $(basename "${output}")"
        return 0
    fi

    echo "[DOWNLOAD] $(basename "${output}")"

    if command -v wget >/dev/null 2>&1; then
        wget \
            --continue \
            --tries=30 \
            --timeout=60 \
            --retry-connrefused \
            --waitretry=5 \
            --output-document="${output}" \
            "${url}"
    else
        curl \
            --location \
            --fail \
            --retry 30 \
            --retry-delay 5 \
            --connect-timeout 60 \
            --continue-at - \
            --output "${output}" \
            "${url}"
    fi
}

BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288887/suppl"

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

echo "===== DOWNLOAD REUSED FILTERED MATRICES ====="

for filename in "${REUSED_FILES[@]}"; do
    download_file \
        "${BASE_URL}/${filename}" \
        "${REUSED288}/${filename}"
done

echo "===== FILTERED MATRIX INVENTORY ====="

find "${EXT288}" "${REUSED288}" \
    -type f \
    -name '*filtered_feature_bc_matrix.h5' \
    -printf '%p\t%s\n' |
    sort |
    tee "${AUDIT}/GSE288887_filtered_file_inventory.tsv"

FILTERED_N=$(
    find "${EXT288}" "${REUSED288}" \
        -type f \
        -name '*filtered_feature_bc_matrix.h5' |
    wc -l
)

RAW_N=$(
    find "${EXT288}" \
        -type f \
        -name '*raw_feature_bc_matrix.h5' |
    wc -l
)

echo "filtered_h5_count=${FILTERED_N}"
echo "raw_h5_excluded=${RAW_N}"

if [[ "${FILTERED_N}" -ne 25 ]]; then
    echo "ERROR: expected 25 currently downloadable filtered H5 files"
    exit 1
fi

cat > "${ROOT}/scripts/02_corrected_filtered_qc_audit.R" <<'RSCRIPT'
suppressPackageStartupMessages({
    library(Matrix)
    library(Seurat)
})

root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

ext288 <- file.path(
    root,
    "01_raw_data",
    "GSE288887",
    "extracted"
)

reused288 <- file.path(
    root,
    "01_raw_data",
    "GSE288887",
    "reused_filtered"
)

raw174 <- file.path(
    root,
    "01_raw_data",
    "GSE174725"
)

ext174 <- file.path(
    raw174,
    "extracted"
)

audit <- file.path(
    root,
    "03_data_audit",
    "corrected_v2"
)

dir.create(audit, recursive = TRUE, showWarnings = FALSE)

qvalue <- function(x, probability) {
    if (length(x) == 0 || all(is.na(x))) {
        return(NA_real_)
    }

    as.numeric(
        quantile(
            x,
            probability,
            na.rm = TRUE,
            names = FALSE
        )
    )
}

safe_median <- function(x) {
    if (length(x) == 0 || all(is.na(x))) {
        return(NA_real_)
    }

    median(x, na.rm = TRUE)
}

extract_expression_matrix <- function(x) {
    if (inherits(x, "Matrix") || is.matrix(x)) {
        return(x)
    }

    if (is.list(x)) {
        if ("Gene Expression" %in% names(x)) {
            return(x[["Gene Expression"]])
        }

        return(x[[1]])
    }

    stop("Unsupported Read10X object")
}

calculate_metrics <- function(mat) {
    total_counts <- Matrix::colSums(mat)
    detected_features <- Matrix::colSums(mat > 0)

    mitochondrial <- grepl(
        "^MT-",
        rownames(mat),
        ignore.case = FALSE
    )

    if (any(mitochondrial)) {
        mitochondrial_counts <- Matrix::colSums(
            mat[mitochondrial, , drop = FALSE]
        )

        pct_mito <- ifelse(
            total_counts > 0,
            100 * mitochondrial_counts / total_counts,
            NA_real_
        )
    } else {
        pct_mito <- rep(NA_real_, ncol(mat))
    }

    mre11_symbols <- intersect(
        c("MRE11", "MRE11A"),
        rownames(mat)
    )

    if (length(mre11_symbols) > 0) {
        mre11_counts <- Matrix::colSums(
            mat[mre11_symbols, , drop = FALSE]
        )

        mre11_detected <- mre11_counts > 0
    } else {
        mre11_detected <- rep(NA, ncol(mat))
    }

    list(
        n_genes = nrow(mat),
        total_counts = total_counts,
        detected_features = detected_features,
        pct_mito = pct_mito,
        mre11_detected = mre11_detected,
        mre11_symbols = paste(
            mre11_symbols,
            collapse = ";"
        )
    )
}

summarize_metrics <- function(
    dataset,
    sample_id,
    subject_id,
    group,
    metrics,
    source_type
) {
    common_pass <- (
        metrics$detected_features >= 200 &
        metrics$pct_mito <= 15
    )

    strict_pass <- (
        metrics$detected_features >= 200 &
        metrics$pct_mito < 10
    )

    if (dataset == "GSE288887") {
        author_basic_pass <- (
            metrics$total_counts > 100 &
            metrics$pct_mito < 10
        )

        author_rule <- "UMI>100;MT<10%;Scrublet_not_available"
    } else {
        author_basic_pass <- common_pass
        author_rule <- "nFeature>=200;MT<=15%"
    }

    data.frame(
        dataset = dataset,
        sample_id = sample_id,
        subject_id = subject_id,
        group = group,
        source_type = source_type,
        n_genes = metrics$n_genes,
        n_cells = length(metrics$total_counts),
        total_umi = sum(metrics$total_counts),
        median_umi = safe_median(metrics$total_counts),
        q05_umi = qvalue(metrics$total_counts, 0.05),
        q95_umi = qvalue(metrics$total_counts, 0.95),
        median_features = safe_median(
            metrics$detected_features
        ),
        q05_features = qvalue(
            metrics$detected_features,
            0.05
        ),
        q95_features = qvalue(
            metrics$detected_features,
            0.95
        ),
        median_pct_mito = safe_median(
            metrics$pct_mito
        ),
        q95_pct_mito = qvalue(
            metrics$pct_mito,
            0.95
        ),
        common_pass_n = sum(
            common_pass,
            na.rm = TRUE
        ),
        common_pass_pct = 100 * mean(
            common_pass,
            na.rm = TRUE
        ),
        strict_pass_n = sum(
            strict_pass,
            na.rm = TRUE
        ),
        strict_pass_pct = 100 * mean(
            strict_pass,
            na.rm = TRUE
        ),
        author_basic_pass_n = sum(
            author_basic_pass,
            na.rm = TRUE
        ),
        author_basic_pass_pct = 100 * mean(
            author_basic_pass,
            na.rm = TRUE
        ),
        author_rule = author_rule,
        mre11_symbol_found = metrics$mre11_symbols,
        mre11_detected_n = sum(
            metrics$mre11_detected,
            na.rm = TRUE
        ),
        mre11_detected_pct = 100 * mean(
            metrics$mre11_detected,
            na.rm = TRUE
        ),
        stringsAsFactors = FALSE
    )
}

metadata288 <- data.frame(
    file_name = c(
        "GSM8778048_4859_2_filtered_feature_bc_matrix.h5",
        "GSM8778049_6281_1_filtered_feature_bc_matrix.h5",
        "GSM8778050_6281_2_filtered_feature_bc_matrix.h5",
        "GSM8778051_7539_filtered_feature_bc_matrix.h5",
        "GSM8778052_8981_1_filtered_feature_bc_matrix.h5",
        "GSM8778053_8981_2_filtered_feature_bc_matrix.h5",
        "GSM8778054_9367_filtered_feature_bc_matrix.h5",
        "GSM8778055_9564_filtered_feature_bc_matrix.h5",
        "GSM8808454_5465_1_filtered_feature_bc_matrix.h5",
        "GSM8808455_5465_2_filtered_feature_bc_matrix.h5",
        "GSE288887_4546_1_filtered_feature_bc_matrix.h5",
        "GSE288887_4546_2_filtered_feature_bc_matrix.h5",
        "GSE288887_4736_1_filtered_feature_bc_matrix.h5",
        "GSE288887_4736_2_filtered_feature_bc_matrix.h5",
        "GSE288887_vuhd068_filtered_feature_bc_matrix.h5",
        "GSE288887_vuhd070_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD092_1_LU_Whole_C2_X5SCR_F02526_HLKH5DSXY_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD094_1_LU_Whole_C2_X5SCR_F02524_HLKH5DSXY_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD098_1_LU_Whole_C2_X5SCR_F02528_HLKH5DSXY_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD103_1_LU_Whole_C1_X5SCR_F02617_HLKH5DSXY_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD105_1_LU_Whole_C1_X5SCR_F02609_HLKH5DSXY_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD106_1_LU_Whole_C1_X5SCR_F02611_HLKH5DSXY_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD108_1_LU_Whole_C1_X5SC2_F03548_H3HLJDSX2-H5W7FDSX2_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD111_1_LU_Whole_C1_X5SC2_F03552_H3HLJDSX2-H5W7FDSX2_filtered_feature_bc_matrix.h5",
        "GSE288887_IPF_VUHD117_1_LU_MNC_C1_X5SC2_F03708_H3HLJDSX2-H5W7FDSX2_filtered_feature_bc_matrix.h5"
    ),
    sample_id = c(
        "CWP89",
        "CWP109_MF",
        "CWP109_LF",
        "CWP127",
        "CWP144_MF",
        "CWP144_LF",
        "CWP149",
        "CWP152",
        "CWP94_MF",
        "CWP94_LF",
        "CWP85_LF",
        "CWP85_MF",
        "CWP86_MF",
        "CWP86_LF",
        "CTL068",
        "CTL070",
        "CTL092",
        "CTL094",
        "CTL098",
        "CTL103",
        "CTL105",
        "CTL106",
        "CTL108",
        "CTL111",
        "CTL117"
    ),
    subject_id = c(
        "CWP89",
        "CWP109",
        "CWP109",
        "CWP127",
        "CWP144",
        "CWP144",
        "CWP149",
        "CWP152",
        "CWP94",
        "CWP94",
        "CWP85",
        "CWP85",
        "CWP86",
        "CWP86",
        "CTL068",
        "CTL070",
        "CTL092",
        "CTL094",
        "CTL098",
        "CTL103",
        "CTL105",
        "CTL106",
        "CTL108",
        "CTL111",
        "CTL117"
    ),
    group = c(
        rep("CWP", 14),
        rep("Control", 11)
    ),
    region = c(
        "unspecified",
        "more_fibrotic",
        "less_fibrotic",
        "unspecified",
        "more_fibrotic",
        "less_fibrotic",
        "unspecified",
        "unspecified",
        "more_fibrotic",
        "less_fibrotic",
        "less_fibrotic",
        "more_fibrotic",
        "more_fibrotic",
        "less_fibrotic",
        rep("nonfibrotic_control", 11)
    ),
    source_type = c(
        rep("GSE288887_new", 10),
        rep("reused_GSE227136", 4),
        rep("reused_GSE135893", 2),
        rep("reused_GSE227136", 9)
    ),
    stringsAsFactors = FALSE
)

write.table(
    metadata288,
    file.path(audit, "GSE288887_sample_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

filtered_files288 <- c(
    list.files(
        ext288,
        pattern = "filtered_feature_bc_matrix\\.h5$",
        full.names = TRUE,
        recursive = TRUE
    ),
    list.files(
        reused288,
        pattern = "filtered_feature_bc_matrix\\.h5$",
        full.names = TRUE,
        recursive = TRUE
    )
)

filtered_files288 <- sort(unique(filtered_files288))

if (length(filtered_files288) != 25) {
    stop(
        paste0(
            "Expected 25 filtered GSE288887 files, found ",
            length(filtered_files288)
        )
    )
}

unmapped_files <- setdiff(
    basename(filtered_files288),
    metadata288$file_name
)

missing_files <- setdiff(
    metadata288$file_name,
    basename(filtered_files288)
)

if (length(unmapped_files) > 0) {
    stop(
        paste(
            "Unmapped GSE288887 files:",
            paste(unmapped_files, collapse = ";")
        )
    )
}

if (length(missing_files) > 0) {
    stop(
        paste(
            "Missing GSE288887 files:",
            paste(missing_files, collapse = ";")
        )
    )
}

cat("===== GSE288887 FILTERED-ONLY AUDIT =====\n")

qc288_list <- vector(
    "list",
    length(filtered_files288)
)

all288_counts <- numeric()
all288_features <- numeric()
all288_mito <- numeric()
all288_mre11 <- logical()
all288_genes <- integer()
all288_symbols <- character()
all288_groups <- character()
all288_subjects <- character()
all288_samples <- character()

for (i in seq_along(filtered_files288)) {
    path <- filtered_files288[[i]]
    file_name <- basename(path)

    meta_index <- match(
        file_name,
        metadata288$file_name
    )

    meta <- metadata288[meta_index, , drop = FALSE]

    cat(
        sprintf(
            "[%d/%d] %s | %s | %s\n",
            i,
            length(filtered_files288),
            meta$sample_id,
            meta$group,
            file_name
        )
    )

    matrix_object <- Read10X_h5(
        filename = path,
        use.names = TRUE,
        unique.features = TRUE
    )

    mat <- extract_expression_matrix(
        matrix_object
    )

    metrics <- calculate_metrics(mat)

    qc288_list[[i]] <- summarize_metrics(
        dataset = "GSE288887",
        sample_id = meta$sample_id,
        subject_id = meta$subject_id,
        group = meta$group,
        metrics = metrics,
        source_type = meta$source_type
    )

    all288_counts <- c(
        all288_counts,
        metrics$total_counts
    )

    all288_features <- c(
        all288_features,
        metrics$detected_features
    )

    all288_mito <- c(
        all288_mito,
        metrics$pct_mito
    )

    all288_mre11 <- c(
        all288_mre11,
        metrics$mre11_detected
    )

    all288_genes <- c(
        all288_genes,
        metrics$n_genes
    )

    if (nzchar(metrics$mre11_symbols)) {
        all288_symbols <- union(
            all288_symbols,
            strsplit(
                metrics$mre11_symbols,
                ";",
                fixed = TRUE
            )[[1]]
        )
    }

    all288_groups <- c(
        all288_groups,
        rep(meta$group, ncol(mat))
    )

    all288_subjects <- c(
        all288_subjects,
        rep(meta$subject_id, ncol(mat))
    )

    all288_samples <- c(
        all288_samples,
        rep(meta$sample_id, ncol(mat))
    )

    rm(
        matrix_object,
        mat,
        metrics
    )

    invisible(gc())
}

qc288 <- do.call(
    rbind,
    qc288_list
)

write.table(
    qc288,
    file.path(audit, "GSE288887_per_sample_qc.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

metrics288_all <- list(
    n_genes = max(all288_genes),
    total_counts = all288_counts,
    detected_features = all288_features,
    pct_mito = all288_mito,
    mre11_detected = all288_mre11,
    mre11_symbols = paste(
        all288_symbols,
        collapse = ";"
    )
)

summary288 <- summarize_metrics(
    dataset = "GSE288887",
    sample_id = "ALL_FILTERED",
    subject_id = paste0(
        length(unique(all288_subjects)),
        "_subjects"
    ),
    group = "CWP_vs_Control",
    metrics = metrics288_all,
    source_type = "25_downloadable_matrices"
)

cat("===== GSE174725 SAMPLE-MAPPING AUDIT =====\n")

matrix174 <- ReadMtx(
    mtx = file.path(
        raw174,
        "GSE174725_matrix.mtx.gz"
    ),
    features = file.path(
        raw174,
        "GSE174725_features.tsv.gz"
    ),
    cells = file.path(
        raw174,
        "GSE174725_barcodes.tsv.gz"
    ),
    feature.column = 2,
    cell.column = 1,
    unique.features = TRUE
)

barcodes174 <- colnames(matrix174)

csv_files174 <- list.files(
    ext174,
    pattern = "\\.csv\\.gz$",
    full.names = TRUE,
    recursive = TRUE
)

if (length(csv_files174) != 5) {
    stop(
        paste0(
            "Expected 5 GSE174725 CSV files, found ",
            length(csv_files174)
        )
    )
}

get_csv_cell_count <- function(path) {
    connection <- gzfile(
        path,
        open = "rt"
    )

    on.exit(
        close(connection),
        add = TRUE
    )

    header <- readLines(
        connection,
        n = 1,
        warn = FALSE
    )

    fields <- strsplit(
        header,
        ",",
        fixed = TRUE
    )[[1]]

    max(length(fields) - 1, 0)
}

sample174 <- data.frame(
    file_name = basename(csv_files174),
    n_cells_csv = vapply(
        csv_files174,
        get_csv_cell_count,
        numeric(1)
    ),
    stringsAsFactors = FALSE
)

sample174$sample_id <- ifelse(
    grepl("GSM5324866", sample174$file_name),
    "EXP1",
    ifelse(
        grepl("GSM5324867", sample174$file_name),
        "EXP2",
        ifelse(
            grepl("GSM5324868", sample174$file_name),
            "SIL1",
            ifelse(
                grepl("GSM5324869", sample174$file_name),
                "SIL2",
                ifelse(
                    grepl("GSM5324870", sample174$file_name),
                    "SIL3",
                    NA_character_
                )
            )
        )
    )
)

sample174$subject_id <- sample174$sample_id

sample174$group <- ifelse(
    grepl("^EXP", sample174$sample_id),
    "Exposed",
    "Silicosis"
)

candidate_suffixes <- list(
    underscore = sub("^.*_", "", barcodes174),
    hyphen = sub("^.*-", "", barcodes174),
    final_numeric = sub(
        "^.*[-_]",
        "",
        barcodes174
    )
)

candidate_report <- list()
selected_name <- NA_character_
selected_suffix <- NULL

for (candidate_name in names(candidate_suffixes)) {
    suffix <- candidate_suffixes[[candidate_name]]
    counts <- table(suffix)

    candidate_report[[candidate_name]] <- data.frame(
        candidate = candidate_name,
        suffix = names(counts),
        n_cells = as.integer(counts),
        stringsAsFactors = FALSE
    )

    exact_count_match <- (
        length(counts) == 5 &&
        identical(
            sort(as.integer(counts)),
            sort(sample174$n_cells_csv)
        )
    )

    if (
        is.na(selected_name) &&
        exact_count_match
    ) {
        selected_name <- candidate_name
        selected_suffix <- suffix
    }
}

candidate_report_df <- do.call(
    rbind,
    candidate_report
)

write.table(
    candidate_report_df,
    file.path(
        audit,
        "GSE174725_barcode_suffix_candidates.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

if (is.na(selected_name)) {
    stop(
        paste(
            "Could not map aggregated GSE174725 barcodes",
            "to five sample CSV files by exact cell counts"
        )
    )
}

suffix_counts <- as.data.frame(
    table(selected_suffix),
    stringsAsFactors = FALSE
)

names(suffix_counts) <- c(
    "suffix",
    "n_cells_matrix"
)

sample174$suffix <- NA_character_

for (i in seq_len(nrow(sample174))) {
    matched_suffix <- suffix_counts$suffix[
        suffix_counts$n_cells_matrix ==
            sample174$n_cells_csv[i]
    ]

    if (length(matched_suffix) != 1) {
        stop(
            paste0(
                "Ambiguous suffix mapping for ",
                sample174$sample_id,
                " with ",
                sample174$n_cells_csv[i],
                " cells"
            )
        )
    }

    sample174$suffix[i] <- matched_suffix
}

sample174$mapping_method <- paste0(
    selected_name,
    "_suffix_exact_cell_count"
)

write.table(
    sample174,
    file.path(
        audit,
        "GSE174725_sample_mapping.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

metrics174 <- calculate_metrics(
    matrix174
)

qc174_list <- vector(
    "list",
    nrow(sample174)
)

all174_group <- rep(
    NA_character_,
    ncol(matrix174)
)

all174_subject <- rep(
    NA_character_,
    ncol(matrix174)
)

for (i in seq_len(nrow(sample174))) {
    selected_cells <- (
        selected_suffix ==
        sample174$suffix[i]
    )

    sample_metrics <- list(
        n_genes = metrics174$n_genes,
        total_counts =
            metrics174$total_counts[selected_cells],
        detected_features =
            metrics174$detected_features[selected_cells],
        pct_mito =
            metrics174$pct_mito[selected_cells],
        mre11_detected =
            metrics174$mre11_detected[selected_cells],
        mre11_symbols =
            metrics174$mre11_symbols
    )

    qc174_list[[i]] <- summarize_metrics(
        dataset = "GSE174725",
        sample_id = sample174$sample_id[i],
        subject_id = sample174$subject_id[i],
        group = sample174$group[i],
        metrics = sample_metrics,
        source_type = "sample_processed_matrix"
    )

    all174_group[selected_cells] <-
        sample174$group[i]

    all174_subject[selected_cells] <-
        sample174$subject_id[i]
}

qc174 <- do.call(
    rbind,
    qc174_list
)

write.table(
    qc174,
    file.path(
        audit,
        "GSE174725_per_sample_qc.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

summary174 <- summarize_metrics(
    dataset = "GSE174725",
    sample_id = "ALL_FILTERED",
    subject_id = "5_subjects",
    group = "Silicosis_vs_Exposed",
    metrics = metrics174,
    source_type = "aggregated_processed_matrix"
)

dataset_summary <- rbind(
    summary288,
    summary174
)

dataset_summary$n_samples <- c(
    nrow(qc288),
    nrow(qc174)
)

dataset_summary$n_subjects <- c(
    length(unique(metadata288$subject_id)),
    length(unique(sample174$subject_id))
)

dataset_summary$n_case_subjects <- c(
    length(
        unique(
            metadata288$subject_id[
                metadata288$group == "CWP"
            ]
        )
    ),
    length(
        unique(
            sample174$subject_id[
                sample174$group == "Silicosis"
            ]
        )
    )
)

dataset_summary$n_control_subjects <- c(
    length(
        unique(
            metadata288$subject_id[
                metadata288$group == "Control"
            ]
        )
    ),
    length(
        unique(
            sample174$subject_id[
                sample174$group == "Exposed"
            ]
        )
    )
)

write.table(
    dataset_summary,
    file.path(
        audit,
        "corrected_dataset_qc_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

make_group_summary <- function(
    dataset,
    groups,
    subjects,
    counts,
    features,
    pct_mito,
    mre11,
    n_genes,
    symbols
) {
    output <- list()

    for (group_name in unique(groups)) {
        selected <- groups == group_name

        group_metrics <- list(
            n_genes = n_genes,
            total_counts = counts[selected],
            detected_features = features[selected],
            pct_mito = pct_mito[selected],
            mre11_detected = mre11[selected],
            mre11_symbols = symbols
        )

        output[[group_name]] <- summarize_metrics(
            dataset = dataset,
            sample_id = "GROUP_TOTAL",
            subject_id = paste0(
                length(unique(subjects[selected])),
                "_subjects"
            ),
            group = group_name,
            metrics = group_metrics,
            source_type = "group_summary"
        )
    }

    do.call(rbind, output)
}

group288 <- make_group_summary(
    dataset = "GSE288887",
    groups = all288_groups,
    subjects = all288_subjects,
    counts = all288_counts,
    features = all288_features,
    pct_mito = all288_mito,
    mre11 = all288_mre11,
    n_genes = max(all288_genes),
    symbols = paste(all288_symbols, collapse = ";")
)

group174 <- make_group_summary(
    dataset = "GSE174725",
    groups = all174_group,
    subjects = all174_subject,
    counts = metrics174$total_counts,
    features = metrics174$detected_features,
    pct_mito = metrics174$pct_mito,
    mre11 = metrics174$mre11_detected,
    n_genes = metrics174$n_genes,
    symbols = metrics174$mre11_symbols
)

group_summary <- rbind(
    group288,
    group174
)

write.table(
    group_summary,
    file.path(
        audit,
        "group_level_qc_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

coverage_check <- data.frame(
    item = c(
        "paper_reported_samples",
        "currently_downloadable_filtered_matrices",
        "paper_reported_cells",
        "corrected_downloadable_cells",
        "paper_reported_CWP_subjects",
        "downloadable_CWP_subjects",
        "paper_reported_control_subjects",
        "downloadable_control_subjects"
    ),
    value = c(
        26,
        nrow(qc288),
        214448,
        summary288$n_cells,
        9,
        length(
            unique(
                metadata288$subject_id[
                    metadata288$group == "CWP"
                ]
            )
        ),
        12,
        length(
            unique(
                metadata288$subject_id[
                    metadata288$group == "Control"
                ]
            )
        )
    ),
    stringsAsFactors = FALSE
)

write.table(
    coverage_check,
    file.path(
        audit,
        "GSE288887_coverage_check.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n===== CORRECTED DATASET SUMMARY =====\n")
print(
    dataset_summary,
    row.names = FALSE
)

cat("\n===== PER-SAMPLE GSE288887 QC =====\n")
print(
    qc288[
        ,
        c(
            "sample_id",
            "subject_id",
            "group",
            "n_cells",
            "median_umi",
            "median_features",
            "median_pct_mito",
            "common_pass_pct",
            "mre11_detected_pct"
        )
    ],
    row.names = FALSE
)

cat("\n===== GSE174725 SAMPLE MAPPING =====\n")
print(
    sample174,
    row.names = FALSE
)

cat("\n===== PER-SAMPLE GSE174725 QC =====\n")
print(
    qc174[
        ,
        c(
            "sample_id",
            "subject_id",
            "group",
            "n_cells",
            "median_umi",
            "median_features",
            "median_pct_mito",
            "common_pass_pct",
            "mre11_detected_pct"
        )
    ],
    row.names = FALSE
)

cat("\n===== GROUP-LEVEL QC =====\n")
print(
    group_summary[
        ,
        c(
            "dataset",
            "group",
            "subject_id",
            "n_cells",
            "median_umi",
            "median_features",
            "median_pct_mito",
            "common_pass_pct",
            "mre11_detected_pct"
        )
    ],
    row.names = FALSE
)

cat("\n===== GSE288887 COVERAGE CHECK =====\n")
print(
    coverage_check,
    row.names = FALSE
)

cat("\n===== IMPORTANT NOTES =====\n")
cat(
    "1. GSE288887 raw_feature_bc_matrix files were excluded.\n"
)
cat(
    "2. GSE288887 author QC excludes Scrublet because doublet scores are unavailable in count matrices.\n"
)
cat(
    "3. Common QC uses nFeature>=200 and MT<=15% for both datasets.\n"
)
cat(
    "4. MRE11 detection is reported but is not treated as a technical quality metric.\n"
)
cat(
    "5. Differential-expression replication units must be subjects, not cells or tissue regions.\n"
)

cat("\n===== STATUS =====\n")
cat("PASS_WITH_COVERAGE_WARNING\n")
RSCRIPT

echo "===== RUN CORRECTED R AUDIT ====="

Rscript \
    "${ROOT}/scripts/02_corrected_filtered_qc_audit.R"

cat > "${ROOT}/00_project/PROJECT_STATUS.md" <<'STATUS'
# Project status

- Project directory created.
- Official metadata comparison completed.
- Initial matrix audit completed but GSE288887 raw and filtered matrices were mixed.
- Corrected filtered-only audit completed.
- GSE288887 reused processed matrices were added.
- GSE174725 per-sample identities were reconstructed from sample CSV files and aggregated barcode suffixes.
- Formal filtering, annotation and differential-expression analysis have not started.
- Differential-expression replication units will be human subjects rather than individual cells.
STATUS

echo "===== COMPLETION ====="
echo "status=PASS_WITH_COVERAGE_WARNING"
date
