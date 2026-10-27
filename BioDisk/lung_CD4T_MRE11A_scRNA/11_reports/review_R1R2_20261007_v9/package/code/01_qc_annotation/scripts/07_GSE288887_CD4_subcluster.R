suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

set.seed(20260812)
options(future.globals.maxSize = 16 * 1024^3)

if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

project_root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

input_rds <- file.path(
  project_root,
  "06_CD4T_analysis/TNK_subcluster_v1",
  "GSE288887_TNK_subclustered.rds"
)

out_dir <- file.path(
  project_root,
  "06_CD4T_analysis/CD4_subcluster_v1"
)

fig_dir <- file.path(
  project_root,
  "09_figures/CD4_subcluster_v1"
)

tab_dir <- file.path(
  project_root,
  "10_tables/CD4_subcluster_v1"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

message("===== READ T/NK OBJECT =====")
tnk <- readRDS(input_rds)
DefaultAssay(tnk) <- "RNA"

if (!"tnk_cluster" %in% colnames(tnk[[]])) {
  stop("tnk_cluster is absent")
}

tnk$tnk_cluster <- as.character(tnk$tnk_cluster)

candidate_clusters <- c("0", "6", "8")

candidate_cells <- rownames(tnk[[]])[
  tnk$tnk_cluster %in% candidate_clusters
]

message("Candidate clusters: ", paste(candidate_clusters, collapse = ","))
message("Candidate cells: ", length(candidate_cells))

print(
  table(
    cluster = tnk$tnk_cluster[candidate_cells],
    group = tnk$group[candidate_cells]
  )
)

message("===== EXTRACT CD4 CANDIDATES =====")
cd4 <- subset(tnk, cells = candidate_cells)
rm(tnk)
gc()

source_map <- c(
  "0" = "ConvCD4",
  "6" = "CD4like",
  "8" = "Treg"
)

cd4$cd4_source <- unname(
  source_map[as.character(cd4$tnk_cluster)]
)

DefaultAssay(cd4) <- "RNA"

if (inherits(cd4[["RNA"]], "Assay5")) {
  message("Joining RNA layers...")
  cd4 <- JoinLayers(cd4, assay = "RNA")
}

message("===== REPROCESS CD4 COMPARTMENT =====")
cd4 <- NormalizeData(
  cd4,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

cd4 <- FindVariableFeatures(
  cd4,
  selection.method = "vst",
  nfeatures = 2500,
  verbose = TRUE
)

pca_features <- VariableFeatures(cd4)

exclude_pattern <- paste0(
  "^(",
  "MT-|RPL|RPS|HBA|HBB|",
  "MALAT1$|XIST$|TSIX$|",
  "RPS4Y1$|DDX3Y$|KDM5D$|",
  "UTY$|EIF1AY$|ZFY$|",
  "MRE11$|MRE11A$",
  ")"
)

excluded_features <- grep(
  exclude_pattern,
  pca_features,
  value = TRUE
)

pca_features <- setdiff(
  pca_features,
  excluded_features
)

VariableFeatures(cd4) <- pca_features

message("PCA features: ", length(pca_features))
message("Excluded PCA features: ", length(excluded_features))

writeLines(
  excluded_features,
  file.path(tab_dir, "CD4_excluded_PCA_features.txt")
)

cd4 <- ScaleData(
  cd4,
  features = pca_features,
  verbose = TRUE
)

cd4 <- RunPCA(
  cd4,
  features = pca_features,
  npcs = 30,
  seed.use = 20260812,
  verbose = TRUE
)

if (!"sample_id" %in% colnames(cd4[[]])) {
  stop("sample_id is absent")
}

message("===== HARMONY =====")
cd4 <- harmony::RunHarmony(
  object = cd4,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = 1:25,
  reduction.save = "harmony",
  plot_convergence = FALSE,
  verbose = TRUE
)

cd4 <- FindNeighbors(
  cd4,
  reduction = "harmony",
  dims = 1:25,
  k.param = 25,
  verbose = TRUE
)

cd4 <- FindClusters(
  cd4,
  resolution = c(0.2, 0.3, 0.4, 0.5, 0.6),
  random.seed = 20260812,
  verbose = TRUE
)

cd4 <- RunUMAP(
  cd4,
  reduction = "harmony",
  dims = 1:25,
  n.neighbors = 25,
  min.dist = 0.30,
  metric = "cosine",
  seed.use = 20260812,
  verbose = TRUE
)

res_col <- "RNA_snn_res.0.4"

if (!res_col %in% colnames(cd4[[]])) {
  stop("Missing clustering column: ", res_col)
}

cd4$cd4_cluster <- as.character(
  cd4[[res_col, drop = TRUE]]
)

Idents(cd4) <- "cd4_cluster"

message("===== CD4 CLUSTER COUNTS =====")
print(table(Idents(cd4)))

output_rds <- file.path(
  out_dir,
  "GSE288887_CD4_candidate_subclustered.rds"
)

saveRDS(
  cd4,
  output_rds,
  compress = FALSE
)

message("Saved: ", output_rds)

message("===== COMPOSITION TABLES =====")
meta <- cd4[[]]

write_matrix_table <- function(tab, filename, row_name) {
  result <- as.data.frame.matrix(tab)
  result[[row_name]] <- rownames(result)
  rownames(result) <- NULL

  result <- result[
    ,
    c(row_name, setdiff(colnames(result), row_name)),
    drop = FALSE
  ]

  write.table(
    result,
    filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_matrix_table(
  table(meta$sample_id, meta$cd4_cluster),
  file.path(tab_dir, "CD4_cluster_counts_by_sample.tsv"),
  "sample_id"
)

group_candidates <- c(
  "group",
  "condition",
  "disease_group"
)

group_col <- group_candidates[
  group_candidates %in% colnames(meta)
][1]

if (!is.na(group_col)) {
  write_matrix_table(
    table(meta[[group_col]], meta$cd4_cluster),
    file.path(tab_dir, "CD4_cluster_counts_by_group.tsv"),
    group_col
  )
}

subject_candidates <- c(
  "subject_id",
  "patient_id",
  "donor_id"
)

subject_col <- subject_candidates[
  subject_candidates %in% colnames(meta)
][1]

if (!is.na(subject_col)) {
  write_matrix_table(
    table(meta[[subject_col]], meta$cd4_cluster),
    file.path(tab_dir, "CD4_cluster_counts_by_subject.tsv"),
    subject_col
  )
}

source_summary <- function(source_col, prefix) {
  tab <- table(
    meta[[source_col]],
    meta$cd4_cluster
  )

  data.frame(
    cluster = colnames(tab),
    cells = colSums(tab),
    source_n = colSums(tab > 0),
    source_ge10_n = colSums(tab >= 10),
    source_ge30_n = colSums(tab >= 30),
    max_fraction = round(
      apply(tab, 2, max) /
        pmax(colSums(tab), 1),
      4
    ),
    source_type = prefix,
    stringsAsFactors = FALSE
  )
}

prevalence_list <- list(
  source_summary("sample_id", "sample")
)

if (!is.na(subject_col)) {
  prevalence_list[[2]] <- source_summary(
    subject_col,
    "subject"
  )
}

prevalence <- do.call(rbind, prevalence_list)

write.table(
  prevalence,
  file.path(tab_dir, "CD4_cluster_prevalence.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("===== UMAP FIGURES =====")
p_cluster <- DimPlot(
  cd4,
  reduction = "umap",
  group.by = "cd4_cluster",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("CD4 clusters") +
  theme_classic(base_size = 13)

ggsave(
  file.path(fig_dir, "GSE288887_CD4_UMAP_clusters.png"),
  p_cluster,
  width = 10,
  height = 8,
  dpi = 240
)

p_source <- DimPlot(
  cd4,
  reduction = "umap",
  group.by = "cd4_source",
  raster = TRUE
) +
  ggtitle("Source") +
  theme_classic(base_size = 13)

ggsave(
  file.path(fig_dir, "GSE288887_CD4_UMAP_source.png"),
  p_source,
  width = 10,
  height = 8,
  dpi = 240
)

if (!is.na(group_col)) {
  p_group <- DimPlot(
    cd4,
    reduction = "umap",
    group.by = group_col,
    split.by = group_col,
    raster = TRUE
  ) +
    ggtitle("Groups") +
    theme_classic(base_size = 13)

  ggsave(
    file.path(fig_dir, "GSE288887_CD4_UMAP_groups.png"),
    p_group,
    width = 14,
    height = 7,
    dpi = 240
  )
}

marker_sets <- list(
  PanT = c(
    "CD3D", "CD3E", "CD3G",
    "TRAC", "TRBC1", "CD247", "LCK"
  ),
  Naive = c(
    "CCR7", "SELL", "TCF7",
    "LEF1", "MAL", "NOSIP"
  ),
  Mem = c(
    "CD4", "IL7R", "LTB",
    "AQP3", "CD40LG", "GPR183", "KLRB1"
  ),
  Treg = c(
    "FOXP3", "IL2RA", "CTLA4",
    "TIGIT", "CCR8", "TNFRSF4", "IKZF2"
  ),
  Th17 = c(
    "CCR6", "RORA",
    "IL23R", "IL17A", "CCL20"
  ),
  Th2 = c(
    "GATA3", "CCR4", "IL4R", "PTGDR2"
  ),
  Tfh = c(
    "CXCR5", "PDCD1", "ICOS",
    "MAF", "TOX", "CXCL13"
  ),
  Cyto = c(
    "CD8A", "CD8B", "CCL5",
    "GZMK", "GZMB", "PRF1", "NKG7", "GNLY"
  ),
  Act = c(
    "HLA-DRA", "HLA-DRB1", "CD69",
    "FOS", "JUN", "JUNB"
  ),
  Cycle = c(
    "MKI67", "TOP2A", "STMN1"
  ),
  Excl = c(
    "KLRD1", "FCGR3A", "TRDC",
    "TYROBP", "LST1", "EPCAM", "PECAM1"
  ),
  Target = c("MRE11")
)

marker_sets <- lapply(
  marker_sets,
  function(x) intersect(x, rownames(cd4))
)

seen_features <- character(0)
for (marker_group in names(marker_sets)) {
  marker_sets[[marker_group]] <- setdiff(
    unique(marker_sets[[marker_group]]),
    seen_features
  )
  seen_features <- c(
    seen_features,
    marker_sets[[marker_group]]
  )
}
stopifnot(
  anyDuplicated(
    unlist(marker_sets, use.names = FALSE)
  ) == 0
)

marker_sets <- marker_sets[
  lengths(marker_sets) > 0
]

p_dot <- DotPlot(
  cd4,
  features = marker_sets,
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle("CD4 markers") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  file.path(fig_dir, "GSE288887_CD4_marker_dotplot.png"),
  p_dot,
  width = 19,
  height = 8,
  dpi = 240
)

feature_genes <- intersect(
  c(
    "CD3D", "CD4", "IL7R", "CCR7",
    "KLRB1", "FOXP3", "IL2RA", "CD8A",
    "NKG7", "HLA-DRA", "MKI67", "MRE11"
  ),
  rownames(cd4)
)

p_feature <- FeaturePlot(
  cd4,
  features = feature_genes,
  reduction = "umap",
  order = TRUE,
  raster = TRUE,
  ncol = 4
)

ggsave(
  file.path(fig_dir, "GSE288887_CD4_key_features.png"),
  p_feature,
  width = 15,
  height = 12,
  dpi = 240
)

message("===== FIND CD4 MARKERS =====")
Idents(cd4) <- "cd4_cluster"

markers <- FindAllMarkers(
  cd4,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.20,
  min.pct = 0.10,
  max.cells.per.ident = 3000,
  random.seed = 20260812,
  verbose = TRUE
)

marker_file <- file.path(
  tab_dir,
  "GSE288887_CD4_subcluster_markers.tsv"
)

write.table(
  markers,
  marker_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

fc_col <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(markers)
)[1]

markers_sig <- markers[
  !is.na(markers$p_val_adj) &
    markers$p_val_adj < 0.05 &
    is.finite(markers[[fc_col]]),
  ,
  drop = FALSE
]

cluster_levels <- as.character(
  sort(unique(as.integer(markers_sig$cluster)))
)

top_list <- lapply(cluster_levels, function(cl) {
  z <- markers_sig[
    as.character(markers_sig$cluster) == cl,
    ,
    drop = FALSE
  ]

  z <- z[
    order(-z[[fc_col]], -z$pct.1),
    ,
    drop = FALSE
  ]

  head(z, 12)
})

top_markers <- do.call(rbind, top_list)

write.table(
  top_markers,
  file.path(tab_dir, "CD4_top12_markers_by_cluster.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("===== TOP 12 MARKERS =====")

for (cl in cluster_levels) {
  z <- top_markers[
    as.character(top_markers$cluster) == cl,
    ,
    drop = FALSE
  ]

  message(
    "Cluster ", cl, " | ",
    paste(z$gene, collapse = ", ")
  )
}

message("===== FINAL SUMMARY =====")
message("Candidate source clusters: 0,6,8")
message("CD4 candidate cells: ", ncol(cd4))
message("CD4 subclusters: ", length(unique(cd4$cd4_cluster)))
message("Object: ", output_rds)
message("Markers: ", marker_file)
message("STATUS=PASS")
