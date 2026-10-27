suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(SingleCellExperiment)
    library(scDblFinder)
    library(data.table)
    library(ggplot2)
})

set.seed(20260812)

root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

data_dirs <- c(
    file.path(
        root,
        "01_raw_data",
        "GSE288887",
        "extracted"
    ),
    file.path(
        root,
        "01_raw_data",
        "GSE288887",
        "reused_filtered"
    )
)

metadata_path <- file.path(
    root,
    "03_data_audit",
    "corrected_v2",
    "GSE288887_sample_metadata.tsv"
)

audit_dir <- file.path(
    root,
    "03_data_audit",
    "formal_preprocessing_v1"
)

object_dir <- file.path(
    root,
    "04_preprocessing"
)

figure_dir <- file.path(
    root,
    "09_figures",
    "formal_preprocessing_v1"
)

table_dir <- file.path(
    root,
    "10_tables",
    "formal_preprocessing_v1"
)

dir.create(
    audit_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    object_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    table_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

required_packages <- c(
    "Seurat",
    "Matrix",
    "SingleCellExperiment",
    "scDblFinder",
    "data.table",
    "ggplot2"
)

package_report <- data.frame(
    package = required_packages,
    available = vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    ),
    version = vapply(
        required_packages,
        function(package_name) {
            if (
                requireNamespace(
                    package_name,
                    quietly = TRUE
                )
            ) {
                as.character(
                    packageVersion(package_name)
                )
            } else {
                NA_character_
            }
        },
        character(1)
    ),
    stringsAsFactors = FALSE
)

fwrite(
    package_report,
    file.path(
        audit_dir,
        "package_versions.tsv"
    ),
    sep = "\t"
)

if (any(!package_report$available)) {
    stop(
        paste(
            "Missing packages:",
            paste(
                package_report$package[
                    !package_report$available
                ],
                collapse = ", "
            )
        )
    )
}

if (!file.exists(metadata_path)) {
    stop(
        paste(
            "Metadata file not found:",
            metadata_path
        )
    )
}

sample_metadata <- fread(
    metadata_path,
    data.table = FALSE
)

required_metadata_columns <- c(
    "file_name",
    "sample_id",
    "subject_id",
    "group",
    "region",
    "source_type"
)

missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(sample_metadata)
)

if (length(missing_metadata_columns) > 0) {
    stop(
        paste(
            "Missing metadata columns:",
            paste(
                missing_metadata_columns,
                collapse = ", "
            )
        )
    )
}

if (nrow(sample_metadata) != 25) {
    stop(
        paste0(
            "Expected 25 metadata rows, found ",
            nrow(sample_metadata)
        )
    )
}

all_h5_files <- unlist(
    lapply(
        data_dirs,
        function(directory) {
            list.files(
                directory,
                pattern = "filtered_feature_bc_matrix\\.h5$",
                full.names = TRUE,
                recursive = TRUE
            )
        }
    )
)

all_h5_files <- unique(
    normalizePath(
        all_h5_files,
        mustWork = TRUE
    )
)

sample_metadata$file_path <- vapply(
    sample_metadata$file_name,
    function(file_name) {
        matched <- all_h5_files[
            basename(all_h5_files) == file_name
        ]

        if (length(matched) != 1) {
            stop(
                paste0(
                    "Expected exactly one file for ",
                    file_name,
                    "; found ",
                    length(matched)
                )
            )
        }

        matched[[1]]
    },
    character(1)
)

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

    stop("Unsupported Read10X_h5 result")
}

safe_median <- function(x) {
    if (length(x) == 0 || all(is.na(x))) {
        return(NA_real_)
    }

    median(
        x,
        na.rm = TRUE
    )
}

safe_quantile <- function(x, probability) {
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

sample_objects <- vector(
    mode = "list",
    length = nrow(sample_metadata)
)

names(sample_objects) <- sample_metadata$sample_id

qc_records <- vector(
    mode = "list",
    length = nrow(sample_metadata)
)

doublet_records <- vector(
    mode = "list",
    length = nrow(sample_metadata)
)

mre11_records <- vector(
    mode = "list",
    length = nrow(sample_metadata)
)

cat("===== FORMAL QC START =====\n")
cat("samples=", nrow(sample_metadata), "\n", sep = "")

for (i in seq_len(nrow(sample_metadata))) {
    meta <- sample_metadata[i, , drop = FALSE]

    cat(
        sprintf(
            "\n[%d/%d] %s | %s | %s\n",
            i,
            nrow(sample_metadata),
            meta$sample_id,
            meta$subject_id,
            meta$group
        )
    )

    matrix_object <- Read10X_h5(
        filename = meta$file_path,
        use.names = TRUE,
        unique.features = TRUE
    )

    counts <- extract_expression_matrix(
        matrix_object
    )

    if (anyDuplicated(rownames(counts)) > 0) {
        stop(
            paste(
                "Duplicated feature names:",
                meta$sample_id
            )
        )
    }

    object <- CreateSeuratObject(
        counts = counts,
        project = "GSE288887",
        assay = "RNA",
        min.cells = 0,
        min.features = 0
    )

    object <- RenameCells(
        object,
        add.cell.id = meta$sample_id
    )

    object$dataset <- "GSE288887"
    object$sample_id <- meta$sample_id
    object$subject_id <- meta$subject_id
    object$group <- meta$group
    object$region <- meta$region
    object$source_type <- meta$source_type

    object[["percent.mt"]] <- PercentageFeatureSet(
        object,
        pattern = "^MT-"
    )

    object[["percent.ribo"]] <- PercentageFeatureSet(
        object,
        pattern = "^RP[SL]"
    )

    input_cells <- ncol(object)

    basic_pass <- (
        object$nFeature_RNA >= 200 &
        object$nCount_RNA > 100 &
        object$percent.mt <= 15
    )

    basic_pass_cells <- sum(
        basic_pass,
        na.rm = TRUE
    )

    qc_records[[i]] <- data.frame(
        sample_id = meta$sample_id,
        subject_id = meta$subject_id,
        group = meta$group,
        region = meta$region,
        input_cells = input_cells,
        basic_pass_cells = basic_pass_cells,
        basic_retention_pct = 100 *
            basic_pass_cells /
            input_cells,
        input_median_umi = safe_median(
            object$nCount_RNA
        ),
        input_median_features = safe_median(
            object$nFeature_RNA
        ),
        input_median_pct_mito = safe_median(
            object$percent.mt
        ),
        input_q95_pct_mito = safe_quantile(
            object$percent.mt,
            0.95
        ),
        stringsAsFactors = FALSE
    )

    object <- subset(
        object,
        cells = colnames(object)[basic_pass]
    )

    if (ncol(object) < 100) {
        stop(
            paste(
                "Fewer than 100 cells after basic QC:",
                meta$sample_id
            )
        )
    }

    cat(
        "basic_qc=",
        input_cells,
        "->",
        ncol(object),
        "\n",
        sep = ""
    )

    sce <- as.SingleCellExperiment(
        object,
        assay = "RNA"
    )

    set.seed(
        20260812 + i
    )

    sce <- scDblFinder(
        sce
    )

    object$doublet_class <- as.character(
        colData(sce)$scDblFinder.class
    )

    object$doublet_score <- as.numeric(
        colData(sce)$scDblFinder.score
    )

    class_table <- table(
        object$doublet_class
    )

    singlet_cells <- sum(
        object$doublet_class == "singlet"
    )

    doublet_cells <- sum(
        object$doublet_class == "doublet"
    )

    doublet_records[[i]] <- data.frame(
        sample_id = meta$sample_id,
        cells_entering_doublet_detection =
            ncol(object),
        singlet_cells = singlet_cells,
        predicted_doublets = doublet_cells,
        predicted_doublet_pct = 100 *
            doublet_cells /
            ncol(object),
        stringsAsFactors = FALSE
    )

    object <- subset(
        object,
        subset = doublet_class == "singlet"
    )

    counts_after <- GetAssayData(
        object,
        assay = "RNA",
        layer = "counts"
    )

    if ("MRE11" %in% rownames(counts_after)) {
        object$MRE11_detected <- as.logical(
            counts_after["MRE11", ] > 0
        )

        mre11_positive <- sum(
            object$MRE11_detected
        )
    } else {
        object$MRE11_detected <- FALSE
        mre11_positive <- 0
    }

    mre11_records[[i]] <- data.frame(
        sample_id = meta$sample_id,
        subject_id = meta$subject_id,
        group = meta$group,
        singlet_cells = ncol(object),
        MRE11_positive_cells = mre11_positive,
        MRE11_detection_pct = 100 *
            mre11_positive /
            ncol(object),
        stringsAsFactors = FALSE
    )

    cat(
        "doublet_filter=",
        singlet_cells + doublet_cells,
        "->",
        singlet_cells,
        "; predicted_doublets=",
        doublet_cells,
        "; MRE11_detection_pct=",
        round(
            mre11_records[[i]]$MRE11_detection_pct,
            3
        ),
        "\n",
        sep = ""
    )

    sample_objects[[meta$sample_id]] <- object

    rm(
        matrix_object,
        counts,
        object,
        sce,
        counts_after
    )

    invisible(gc())
}

qc_summary <- do.call(
    rbind,
    qc_records
)

doublet_summary <- do.call(
    rbind,
    doublet_records
)

mre11_summary <- do.call(
    rbind,
    mre11_records
)

sample_summary <- merge(
    qc_summary,
    doublet_summary,
    by = "sample_id",
    all = TRUE,
    sort = FALSE
)

sample_summary <- merge(
    sample_summary,
    mre11_summary[
        ,
        c(
            "sample_id",
            "singlet_cells",
            "MRE11_positive_cells",
            "MRE11_detection_pct"
        )
    ],
    by = c(
        "sample_id",
        "singlet_cells"
    ),
    all = TRUE,
    sort = FALSE
)

sample_summary$final_retention_pct <- 100 *
    sample_summary$singlet_cells /
    sample_summary$input_cells

sample_summary$sample_qc_flag <- "PASS"

sample_summary$sample_qc_flag[
    sample_summary$basic_retention_pct < 50
] <- "LOW_RETENTION"

sample_summary$sample_qc_flag[
    sample_summary$singlet_cells < 1000
] <- "LOW_CELL_COUNT"

sample_summary$sample_qc_flag[
    sample_summary$basic_retention_pct < 50 &
    sample_summary$singlet_cells < 1000
] <- "LOW_RETENTION_AND_COUNT"

sample_summary <- sample_summary[
    match(
        sample_metadata$sample_id,
        sample_summary$sample_id
    ),
    ,
    drop = FALSE
]

fwrite(
    sample_summary,
    file.path(
        audit_dir,
        "GSE288887_formal_sample_qc.tsv"
    ),
    sep = "\t"
)

fwrite(
    doublet_summary,
    file.path(
        audit_dir,
        "GSE288887_doublet_summary.tsv"
    ),
    sep = "\t"
)

fwrite(
    mre11_summary,
    file.path(
        table_dir,
        "GSE288887_MRE11_postQC_by_sample.tsv"
    ),
    sep = "\t"
)

cat("\n===== MERGE CLEANED SAMPLE OBJECTS =====\n")

merged <- merge(
    x = sample_objects[[1]],
    y = sample_objects[-1],
    project = "GSE288887",
    merge.data = FALSE
)

DefaultAssay(merged) <- "RNA"

if (
    "JoinLayers" %in%
    getNamespaceExports("SeuratObject")
) {
    merged <- JoinLayers(
        merged,
        assay = "RNA"
    )
}

merged$sample_id <- factor(
    merged$sample_id,
    levels = sample_metadata$sample_id
)

merged$subject_id <- factor(
    merged$subject_id,
    levels = unique(sample_metadata$subject_id)
)

merged$group <- factor(
    merged$group,
    levels = c(
        "Control",
        "CWP"
    )
)

merged$region <- factor(
    merged$region
)

cell_metadata <- merged[[]]

cell_metadata$cell_id <- rownames(
    cell_metadata
)

fwrite(
    cell_metadata,
    file.path(
        audit_dir,
        "GSE288887_clean_cell_metadata.tsv.gz"
    ),
    sep = "\t"
)

group_summary <- data.frame(
    group = levels(merged$group),
    subjects = vapply(
        levels(merged$group),
        function(group_name) {
            length(
                unique(
                    as.character(
                        merged$subject_id[
                            merged$group == group_name
                        ]
                    )
                )
            )
        },
        numeric(1)
    ),
    samples = vapply(
        levels(merged$group),
        function(group_name) {
            length(
                unique(
                    as.character(
                        merged$sample_id[
                            merged$group == group_name
                        ]
                    )
                )
            )
        },
        numeric(1)
    ),
    singlet_cells = vapply(
        levels(merged$group),
        function(group_name) {
            sum(
                merged$group == group_name
            )
        },
        numeric(1)
    ),
    MRE11_positive_cells = vapply(
        levels(merged$group),
        function(group_name) {
            sum(
                merged$MRE11_detected[
                    merged$group == group_name
                ]
            )
        },
        numeric(1)
    ),
    stringsAsFactors = FALSE
)

group_summary$MRE11_detection_pct <- 100 *
    group_summary$MRE11_positive_cells /
    group_summary$singlet_cells

fwrite(
    group_summary,
    file.path(
        table_dir,
        "GSE288887_postQC_group_summary.tsv"
    ),
    sep = "\t"
)

object_path <- file.path(
    object_dir,
    "GSE288887_formal_QC_singlets.rds"
)

cat(
    "saving_object=",
    object_path,
    "\n",
    sep = ""
)

saveRDS(
    merged,
    object_path,
    compress = FALSE
)

object_md5 <- tools::md5sum(
    object_path
)

write.table(
    data.frame(
        file = names(object_md5),
        md5 = unname(object_md5),
        size_bytes = file.info(
            object_path
        )$size,
        cells = ncol(merged),
        genes = nrow(merged),
        stringsAsFactors = FALSE
    ),
    file.path(
        audit_dir,
        "GSE288887_object_manifest.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

plot_data <- sample_summary

plot_data$sample_id <- factor(
    plot_data$sample_id,
    levels = plot_data$sample_id
)

retention_plot <- ggplot(
    plot_data,
    aes(
        x = sample_id,
        y = final_retention_pct,
        fill = group
    )
) +
    geom_col(
        width = 0.82
    ) +
    geom_hline(
        yintercept = 50,
        linetype = 2,
        colour = "#C0392B",
        linewidth = 0.6
    ) +
    scale_fill_manual(
        values = c(
            "Control" = "#4C78A8",
            "CWP" = "#E45756"
        )
    ) +
    labs(
        x = "Sample",
        y = "Retention (%)",
        fill = "Group",
        title = "Post-QC cell retention"
    ) +
    theme_classic(
        base_size = 11
    ) +
    theme(
        axis.text.x = element_text(
            angle = 60,
            hjust = 1,
            vjust = 1
        ),
        legend.position = "top"
    )

ggsave(
    filename = file.path(
        figure_dir,
        "GSE288887_postQC_retention.png"
    ),
    plot = retention_plot,
    width = 11,
    height = 5.5,
    dpi = 300
)

ggsave(
    filename = file.path(
        figure_dir,
        "GSE288887_postQC_retention.pdf"
    ),
    plot = retention_plot,
    width = 11,
    height = 5.5
)

writeLines(
    capture.output(
        sessionInfo()
    ),
    file.path(
        audit_dir,
        "sessionInfo.txt"
    )
)

cat("\n===== FORMAL QC SUMMARY =====\n")
print(
    sample_summary[
        ,
        c(
            "sample_id",
            "subject_id",
            "group",
            "input_cells",
            "basic_pass_cells",
            "singlet_cells",
            "basic_retention_pct",
            "predicted_doublet_pct",
            "final_retention_pct",
            "MRE11_detection_pct",
            "sample_qc_flag"
        )
    ],
    row.names = FALSE
)

cat("\n===== GROUP SUMMARY =====\n")
print(
    group_summary,
    row.names = FALSE
)

cat("\n===== OBJECT =====\n")
cat(
    "path=",
    object_path,
    "\n",
    sep = ""
)
cat(
    "genes=",
    nrow(merged),
    "\n",
    sep = ""
)
cat(
    "cells=",
    ncol(merged),
    "\n",
    sep = ""
)
cat(
    "size_bytes=",
    file.info(object_path)$size,
    "\n",
    sep = ""
)
cat(
    "md5=",
    unname(object_md5),
    "\n",
    sep = ""
)

cat("\n===== STATUS =====\n")
cat("PASS\n")
