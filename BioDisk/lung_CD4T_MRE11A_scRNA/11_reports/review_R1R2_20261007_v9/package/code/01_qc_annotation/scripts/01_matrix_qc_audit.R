suppressPackageStartupMessages({
    library(Matrix)
    library(Seurat)
})

root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
ext288 <- file.path(root, "01_raw_data", "GSE288887", "extracted")
raw174 <- file.path(root, "01_raw_data", "GSE174725")
audit_dir <- file.path(root, "03_data_audit")

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

qvalue <- function(x, probability) {
    if (all(is.na(x))) {
        return(NA_real_)
    }
    as.numeric(quantile(x, probability, na.rm = TRUE, names = FALSE))
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
    cell_counts <- Matrix::colSums(mat)
    cell_features <- Matrix::colSums(mat > 0)

    mitochondrial <- grepl("^MT-", rownames(mat))

    if (any(mitochondrial)) {
        mitochondrial_counts <- Matrix::colSums(
            mat[mitochondrial, , drop = FALSE]
        )
        pct_mito <- ifelse(
            cell_counts > 0,
            100 * mitochondrial_counts / cell_counts,
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
        counts = cell_counts,
        features = cell_features,
        pct_mito = pct_mito,
        mre11_detected = mre11_detected,
        mre11_symbols = paste(mre11_symbols, collapse = ";")
    )
}

summarize_metrics <- function(dataset, sample_id, metrics, rule) {
    if (rule == "GSE288887") {
        basic_retain <- metrics$counts > 100 &
            metrics$pct_mito < 10
    } else {
        basic_retain <- metrics$features >= 200 &
            metrics$pct_mito <= 15
    }

    data.frame(
        dataset = dataset,
        sample_id = sample_id,
        n_genes = metrics$n_genes,
        n_cells = length(metrics$counts),
        total_umi = sum(metrics$counts),
        median_umi = median(metrics$counts),
        median_features = median(metrics$features),
        q05_features = qvalue(metrics$features, 0.05),
        q95_features = qvalue(metrics$features, 0.95),
        median_pct_mito = median(metrics$pct_mito, na.rm = TRUE),
        q95_pct_mito = qvalue(metrics$pct_mito, 0.95),
        basic_retain_n = sum(basic_retain, na.rm = TRUE),
        basic_retain_pct = 100 * mean(basic_retain, na.rm = TRUE),
        mre11_symbol_found = metrics$mre11_symbols,
        mre11_detected_n = sum(metrics$mre11_detected, na.rm = TRUE),
        mre11_detected_pct = 100 * mean(
            metrics$mre11_detected,
            na.rm = TRUE
        ),
        stringsAsFactors = FALSE
    )
}

cat("===== GSE288887 H5 AUDIT =====\n")

h5_files <- list.files(
    ext288,
    pattern = "\\.h5$",
    full.names = TRUE,
    recursive = TRUE
)

if (length(h5_files) == 0) {
    stop("No H5 files found for GSE288887")
}

qc288 <- list()
all288_counts <- numeric()
all288_features <- numeric()
all288_mito <- numeric()
all288_mre11 <- logical()
all288_symbols <- character()
genes288 <- integer()

for (i in seq_along(h5_files)) {
    path <- h5_files[[i]]
    sample_id <- basename(path)
    sample_id <- sub("^GSM[0-9]+_", "", sample_id)
    sample_id <- sub("_filtered_feature_bc_matrix\\.h5$", "", sample_id)
    sample_id <- sub("\\.h5$", "", sample_id)

    cat(
        sprintf(
            "[%d/%d] %s\n",
            i,
            length(h5_files),
            basename(path)
        )
    )

    matrix_object <- Read10X_h5(
        filename = path,
        use.names = TRUE,
        unique.features = TRUE
    )

    mat <- extract_expression_matrix(matrix_object)
    metrics <- calculate_metrics(mat)

    qc288[[i]] <- summarize_metrics(
        dataset = "GSE288887",
        sample_id = sample_id,
        metrics = metrics,
        rule = "GSE288887"
    )

    all288_counts <- c(all288_counts, metrics$counts)
    all288_features <- c(all288_features, metrics$features)
    all288_mito <- c(all288_mito, metrics$pct_mito)
    all288_mre11 <- c(all288_mre11, metrics$mre11_detected)
    all288_symbols <- union(
        all288_symbols,
        strsplit(metrics$mre11_symbols, ";", fixed = TRUE)[[1]]
    )
    genes288 <- c(genes288, metrics$n_genes)

    rm(matrix_object, mat, metrics)
    invisible(gc())
}

qc288 <- do.call(rbind, qc288)

write.table(
    qc288,
    file.path(audit_dir, "GSE288887_per_matrix_qc.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

metrics288_all <- list(
    n_genes = max(genes288),
    counts = all288_counts,
    features = all288_features,
    pct_mito = all288_mito,
    mre11_detected = all288_mre11,
    mre11_symbols = paste(
        setdiff(all288_symbols, ""),
        collapse = ";"
    )
)

summary288 <- summarize_metrics(
    dataset = "GSE288887",
    sample_id = "ALL_MATRICES",
    metrics = metrics288_all,
    rule = "GSE288887"
)

cat("===== GSE174725 MATRIX AUDIT =====\n")

matrix174 <- ReadMtx(
    mtx = file.path(raw174, "GSE174725_matrix.mtx.gz"),
    features = file.path(raw174, "GSE174725_features.tsv.gz"),
    cells = file.path(raw174, "GSE174725_barcodes.tsv.gz"),
    feature.column = 2,
    cell.column = 1,
    unique.features = TRUE
)

metrics174 <- calculate_metrics(matrix174)

summary174 <- summarize_metrics(
    dataset = "GSE174725",
    sample_id = "AGGREGATED_MATRIX",
    metrics = metrics174,
    rule = "GSE174725"
)

dataset_summary <- rbind(summary288, summary174)

write.table(
    dataset_summary,
    file.path(audit_dir, "dataset_qc_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

barcodes174 <- readLines(
    gzfile(file.path(raw174, "GSE174725_barcodes.tsv.gz"))
)

suffix_hyphen <- sub("^.*-", "", barcodes174)
suffix_underscore <- sub("^.*_", "", barcodes174)

hyphen_table <- as.data.frame(table(suffix_hyphen))
names(hyphen_table) <- c("suffix_after_hyphen", "n_cells")

underscore_table <- as.data.frame(table(suffix_underscore))
names(underscore_table) <- c("suffix_after_underscore", "n_cells")

write.table(
    hyphen_table,
    file.path(audit_dir, "GSE174725_barcode_hyphen_suffix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    underscore_table,
    file.path(audit_dir, "GSE174725_barcode_underscore_suffix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

gene_check <- data.frame(
    dataset = c("GSE288887", "GSE174725"),
    requested_name = c("MRE11A", "MRE11A"),
    official_symbol = c("MRE11", "MRE11"),
    symbol_found = c(
        summary288$mre11_symbol_found,
        summary174$mre11_symbol_found
    ),
    detected_cells = c(
        summary288$mre11_detected_n,
        summary174$mre11_detected_n
    ),
    detected_pct = c(
        summary288$mre11_detected_pct,
        summary174$mre11_detected_pct
    )
)

write.table(
    gene_check,
    file.path(audit_dir, "MRE11_gene_name_check.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n===== DATASET QC SUMMARY =====\n")
print(dataset_summary, row.names = FALSE)

cat("\n===== MRE11 GENE CHECK =====\n")
print(gene_check, row.names = FALSE)

cat("\n===== H5 FILE COUNT =====\n")
cat("GSE288887_h5_files=", length(h5_files), "\n", sep = "")

cat("\n===== STATUS =====\n")
cat("PASS\n")
