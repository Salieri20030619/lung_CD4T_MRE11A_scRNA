suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(data.table)
    library(ggplot2)
})

set.seed(20260813)

root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

input_path <- file.path(
    root,
    "04_preprocessing",
    "GSE288887_formal_QC_singlets.rds"
)

output_path <- file.path(
    root,
    "04_preprocessing",
    "GSE288887_broad_clustered.rds"
)

audit_dir <- file.path(
    root,
    "03_data_audit",
    "broad_clustering_v1"
)

table_dir <- file.path(
    root,
    "10_tables",
    "broad_clustering_v1"
)

figure_dir <- file.path(
    root,
    "09_figures",
    "broad_clustering_v1"
)

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("harmony", quietly = TRUE)) {
    stop(
        paste(
            "R package harmony is unavailable.",
            "Do not continue without batch correction."
        )
    )
}

if (!file.exists(input_path)) {
    stop(
        paste(
            "Input object not found:",
            input_path
        )
    )
}

cat("===== READ QC OBJECT =====\n")

object <- readRDS(input_path)

DefaultAssay(object) <- "RNA"

required_metadata <- c(
    "sample_id",
    "subject_id",
    "group",
    "region",
    "doublet_class"
)

missing_metadata <- setdiff(
    required_metadata,
    colnames(object[[]])
)

if (length(missing_metadata) > 0) {
    stop(
        paste(
            "Missing metadata:",
            paste(missing_metadata, collapse = ",")
        )
    )
}

if (ncol(object) != 192575) {
    stop(
        paste0(
            "Expected 192575 cells, found ",
            ncol(object)
        )
    )
}

cat("genes=", nrow(object), "\n", sep = "")
cat("cells=", ncol(object), "\n", sep = "")
cat("samples=", length(unique(object$sample_id)), "\n", sep = "")
cat("subjects=", length(unique(object$subject_id)), "\n", sep = "")

cat("===== NORMALIZATION =====\n")

object <- NormalizeData(
    object,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = TRUE
)

object <- FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = TRUE
)

variable_features <- VariableFeatures(object)

writeLines(
    variable_features,
    file.path(
        table_dir,
        "GSE288887_variable_features.txt"
    )
)

cat("variable_features=", length(variable_features), "\n", sep = "")

cat("===== SCALE AND PCA =====\n")

object <- ScaleData(
    object,
    features = variable_features,
    verbose = TRUE
)

object <- RunPCA(
    object,
    features = variable_features,
    npcs = 50,
    approx = TRUE,
    seed.use = 20260813,
    verbose = TRUE
)

pca_sd <- Stdev(
    object,
    reduction = "pca"
)

pca_variance <- 100 * pca_sd^2 / sum(pca_sd^2)
pca_cumulative <- cumsum(pca_variance)

pca_report <- data.frame(
    PC = seq_along(pca_sd),
    standard_deviation = pca_sd,
    variance_pct = pca_variance,
    cumulative_pct = pca_cumulative
)

x <- pca_report$PC
y <- pca_report$variance_pct

line_y <- y[1] + (
    y[length(y)] - y[1]
) * (
    x - x[1]
) / (
    x[length(x)] - x[1]
)

raw_elbow <- which.max(
    line_y - y
)

n_dims <- max(
    20,
    min(
        40,
        raw_elbow
    )
)

pca_report$selected <- (
    pca_report$PC <= n_dims
)

fwrite(
    pca_report,
    file.path(
        audit_dir,
        "PCA_variance.tsv"
    ),
    sep = "\t"
)

cat("raw_elbow=", raw_elbow, "\n", sep = "")
cat("selected_dims=", n_dims, "\n", sep = "")

cat("===== HARMONY =====\n")

object <- harmony::RunHarmony(
    object = object,
    group.by.vars = "sample_id",
    reduction.use = "pca",
    dims.use = 1:50,
    reduction.save = "harmony",
    plot_convergence = FALSE,
    verbose = TRUE
)

if (!"harmony" %in% Reductions(object)) {
    stop("Harmony reduction was not created")
}

cat("===== NEIGHBORS AND UMAP =====\n")

object <- FindNeighbors(
    object,
    reduction = "harmony",
    dims = seq_len(n_dims),
    k.param = 30,
    graph.name = c(
        "RNA_nn",
        "RNA_snn"
    ),
    verbose = TRUE
)

object <- RunUMAP(
    object,
    reduction = "harmony",
    dims = seq_len(n_dims),
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = 20260813,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = TRUE
)

cat("===== MULTI-RESOLUTION CLUSTERING =====\n")

object <- FindClusters(
    object,
    graph.name = "RNA_snn",
    resolution = c(
        0.2,
        0.4,
        0.6,
        0.8
    ),
    algorithm = 1,
    random.seed = 20260813,
    verbose = TRUE
)

cluster_column <- "RNA_snn_res.0.4"

if (!cluster_column %in% colnames(object[[]])) {
    stop(
        paste(
            "Cluster column missing:",
            cluster_column
        )
    )
}

Idents(object) <- cluster_column

cluster_summary <- as.data.table(
    object[[]],
    keep.rownames = "cell_id"
)[
    ,
    .(
        cells = .N,
        subjects = uniqueN(subject_id),
        samples = uniqueN(sample_id),
        CWP_cells = sum(group == "CWP"),
        Control_cells = sum(group == "Control")
    ),
    by = .(
        cluster = get(cluster_column)
    )
][
    order(as.integer(as.character(cluster)))
]

cluster_summary[
    ,
    CWP_pct := 100 * CWP_cells / cells
]

cluster_summary[
    ,
    Control_pct := 100 * Control_cells / cells
]

fwrite(
    cluster_summary,
    file.path(
        table_dir,
        "GSE288887_cluster_summary_res04.tsv"
    ),
    sep = "\t"
)

mixing_summary <- as.data.table(
    object[[]]
)[
    ,
    .(
        cells = .N,
        subjects = uniqueN(subject_id),
        samples = uniqueN(sample_id)
    ),
    by = c(
        cluster_column,
        "group"
    )
]

setnames(
    mixing_summary,
    cluster_column,
    "cluster"
)

fwrite(
    mixing_summary,
    file.path(
        table_dir,
        "GSE288887_cluster_group_mixing.tsv"
    ),
    sep = "\t"
)

cat("===== SAVE CLUSTERED OBJECT =====\n")

saveRDS(
    object,
    output_path,
    compress = FALSE
)

object_md5 <- tools::md5sum(output_path)

write.table(
    data.frame(
        path = output_path,
        genes = nrow(object),
        cells = ncol(object),
        clusters_res04 = length(levels(Idents(object))),
        selected_dims = n_dims,
        size_bytes = file.info(output_path)$size,
        md5 = unname(object_md5),
        stringsAsFactors = FALSE
    ),
    file.path(
        audit_dir,
        "GSE288887_broad_object_manifest.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cluster_plot <- DimPlot(
    object,
    reduction = "umap",
    group.by = cluster_column,
    label = TRUE,
    repel = TRUE,
    raster = TRUE
) +
    labs(
        title = "Broad clusters",
        colour = "Cluster"
    ) +
    theme_classic(base_size = 11)

group_plot <- DimPlot(
    object,
    reduction = "umap",
    group.by = "group",
    raster = TRUE,
    cols = c(
        "Control" = "#4C78A8",
        "CWP" = "#E45756"
    )
) +
    labs(
        title = "Disease groups",
        colour = "Group"
    ) +
    theme_classic(base_size = 11)

sample_plot <- DimPlot(
    object,
    reduction = "umap",
    group.by = "sample_id",
    raster = TRUE
) +
    labs(
        title = "Sample mixing",
        colour = "Sample"
    ) +
    theme_classic(base_size = 10) +
    theme(
        legend.position = "none"
    )

ggsave(
    file.path(
        figure_dir,
        "GSE288887_UMAP_clusters_res04.png"
    ),
    cluster_plot,
    width = 8,
    height = 6.5,
    dpi = 300
)

ggsave(
    file.path(
        figure_dir,
        "GSE288887_UMAP_groups.png"
    ),
    group_plot,
    width = 8,
    height = 6.5,
    dpi = 300
)

ggsave(
    file.path(
        figure_dir,
        "GSE288887_UMAP_samples.png"
    ),
    sample_plot,
    width = 8,
    height = 6.5,
    dpi = 300
)

writeLines(
    capture.output(sessionInfo()),
    file.path(
        audit_dir,
        "sessionInfo_04a.txt"
    )
)

cat("\n===== CLUSTER SUMMARY =====\n")
print(cluster_summary)

cat("\n===== OBJECT =====\n")
cat("path=", output_path, "\n", sep = "")
cat("genes=", nrow(object), "\n", sep = "")
cat("cells=", ncol(object), "\n", sep = "")
cat("clusters_res04=", length(levels(Idents(object))), "\n", sep = "")
cat("selected_dims=", n_dims, "\n", sep = "")
cat("size_bytes=", file.info(output_path)$size, "\n", sep = "")
cat("md5=", unname(object_md5), "\n", sep = "")

cat("\n===== STATUS =====\n")
cat("PASS\n")
