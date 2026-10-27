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
    "GSE288887_broad_clustered.rds"
)

output_path <- file.path(
    root,
    "04_preprocessing",
    "GSE288887_broad_clustered_scored.rds"
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

audit_dir <- file.path(
    root,
    "03_data_audit",
    "broad_clustering_v1"
)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

cat("===== READ CLUSTERED OBJECT =====\n")

object <- readRDS(input_path)

DefaultAssay(object) <- "RNA"

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

marker_sets <- list(
    T = c(
        "CD3D",
        "CD3E",
        "TRAC",
        "LCK",
        "IL7R"
    ),
    NK = c(
        "NKG7",
        "GNLY",
        "PRF1",
        "KLRD1"
    ),
    B = c(
        "MS4A1",
        "CD79A",
        "CD37",
        "CD74"
    ),
    My = c(
        "LST1",
        "TYROBP",
        "FCER1G",
        "CTSS"
    ),
    Epi = c(
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19"
    ),
    Endo = c(
        "PECAM1",
        "VWF",
        "EMCN",
        "CLDN5"
    ),
    Str = c(
        "COL1A1",
        "COL1A2",
        "DCN",
        "LUM"
    ),
    Mast = c(
        "TPSAB1",
        "TPSB2",
        "KIT",
        "MS4A2"
    )
)

marker_sets <- lapply(
    marker_sets,
    function(genes) {
        intersect(
            genes,
            rownames(object)
        )
    }
)

if (any(lengths(marker_sets) < 2)) {
    stop(
        paste(
            "Insufficient marker genes for:",
            paste(
                names(marker_sets)[
                    lengths(marker_sets) < 2
                ],
                collapse = ","
            )
        )
    )
}

cat("===== LINEAGE SCORES =====\n")

expression_data <- GetAssayData(
    object,
    assay = "RNA",
    layer = "data"
)

score_columns <- character()

for (lineage in names(marker_sets)) {
    column_name <- paste0(
        "score_",
        lineage
    )

    object[[column_name]] <- Matrix::colMeans(
        expression_data[
            marker_sets[[lineage]],
            ,
            drop = FALSE
        ]
    )

    score_columns <- c(
        score_columns,
        column_name
    )
}

metadata <- as.data.table(
    object[[]],
    keep.rownames = "cell_id"
)

cluster_scores <- metadata[
    ,
    lapply(
        .SD,
        mean,
        na.rm = TRUE
    ),
    by = .(
        cluster = get(cluster_column)
    ),
    .SDcols = score_columns
]

score_matrix <- as.matrix(
    cluster_scores[
        ,
        ..score_columns
    ]
)

best_index <- max.col(
    score_matrix,
    ties.method = "first"
)

cluster_scores[
    ,
    lineage_prelim := names(marker_sets)[best_index]
]

score_margin <- apply(
    score_matrix,
    1,
    function(values) {
        ordered <- sort(
            values,
            decreasing = TRUE
        )

        ordered[1] - ordered[2]
    }
)

cluster_scores[
    ,
    score_margin := score_margin
]

fwrite(
    cluster_scores,
    file.path(
        table_dir,
        "GSE288887_preliminary_lineage_scores.tsv"
    ),
    sep = "\t"
)

lineage_map <- setNames(
    cluster_scores$lineage_prelim,
    as.character(cluster_scores$cluster)
)

object$lineage_prelim <- unname(
    lineage_map[
        as.character(
            object[[cluster_column]][, 1]
        )
    ]
)

object$lineage_prelim <- factor(
    object$lineage_prelim,
    levels = names(marker_sets)
)

lineage_summary <- as.data.table(
    object[[]]
)[
    ,
    .(
        cells = .N,
        subjects = uniqueN(subject_id),
        samples = uniqueN(sample_id),
        CWP_cells = sum(group == "CWP"),
        Control_cells = sum(group == "Control")
    ),
    by = lineage_prelim
][
    order(-cells)
]

fwrite(
    lineage_summary,
    file.path(
        table_dir,
        "GSE288887_preliminary_lineage_summary.tsv"
    ),
    sep = "\t"
)

t_candidate_clusters <- cluster_scores[
    lineage_prelim == "T",
    .(
        cluster,
        lineage_prelim,
        score_T,
        score_NK,
        score_margin
    )
]

fwrite(
    t_candidate_clusters,
    file.path(
        table_dir,
        "GSE288887_T_candidate_clusters.tsv"
    ),
    sep = "\t"
)

cat("===== CLUSTER MARKERS =====\n")

cluster_markers <- FindAllMarkers(
    object,
    assay = "RNA",
    only.pos = TRUE,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.10,
    max.cells.per.ident = 5000,
    random.seed = 20260813,
    verbose = TRUE
)

fwrite(
    cluster_markers,
    file.path(
        table_dir,
        "GSE288887_cluster_markers_res04.tsv.gz"
    ),
    sep = "\t"
)

top_markers <- as.data.table(
    cluster_markers
)[
    order(cluster, -avg_log2FC)
][
    ,
    head(.SD, 20),
    by = cluster
]

fwrite(
    top_markers,
    file.path(
        table_dir,
        "GSE288887_top20_markers_res04.tsv"
    ),
    sep = "\t"
)

cat("===== FIGURES =====\n")

marker_order <- unique(
    unlist(
        marker_sets,
        use.names = FALSE
    )
)

dot_plot <- DotPlot(
    object,
    features = marker_order,
    group.by = cluster_column,
    assay = "RNA"
) +
    RotatedAxis() +
    labs(
        x = "Marker",
        y = "Cluster",
        title = "Lineage markers"
    ) +
    theme_classic(base_size = 10)

lineage_plot <- DimPlot(
    object,
    reduction = "umap",
    group.by = "lineage_prelim",
    raster = TRUE,
    cols = c(
        "T" = "#4C78A8",
        "NK" = "#72B7B2",
        "B" = "#F2CF5B",
        "My" = "#E45756",
        "Epi" = "#B279A2",
        "Endo" = "#54A24B",
        "Str" = "#FF9DA6",
        "Mast" = "#9D755D"
    )
) +
    labs(
        title = "Prelim lineages",
        colour = "Lineage"
    ) +
    theme_classic(base_size = 11)

ggsave(
    file.path(
        figure_dir,
        "GSE288887_lineage_marker_dotplot.png"
    ),
    dot_plot,
    width = 12,
    height = 7,
    dpi = 300
)

ggsave(
    file.path(
        figure_dir,
        "GSE288887_preliminary_lineages.png"
    ),
    lineage_plot,
    width = 8,
    height = 6.5,
    dpi = 300
)

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
        clusters = length(levels(Idents(object))),
        marker_rows = nrow(cluster_markers),
        size_bytes = file.info(output_path)$size,
        md5 = unname(object_md5),
        stringsAsFactors = FALSE
    ),
    file.path(
        audit_dir,
        "GSE288887_scored_object_manifest.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n===== PRELIMINARY LINEAGES =====\n")
print(lineage_summary)

cat("\n===== T CANDIDATE CLUSTERS =====\n")
print(t_candidate_clusters)

cat("\n===== MARKERS =====\n")
cat("marker_rows=", nrow(cluster_markers), "\n", sep = "")
cat("clusters=", length(levels(Idents(object))), "\n", sep = "")

cat("\n===== OBJECT =====\n")
cat("path=", output_path, "\n", sep = "")
cat("cells=", ncol(object), "\n", sep = "")
cat("size_bytes=", file.info(output_path)$size, "\n", sep = "")
cat("md5=", unname(object_md5), "\n", sep = "")

cat("\n===== STATUS =====\n")
cat("PASS\n")
