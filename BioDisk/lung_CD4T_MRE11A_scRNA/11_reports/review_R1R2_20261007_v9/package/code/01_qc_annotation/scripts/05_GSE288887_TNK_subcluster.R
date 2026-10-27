suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

set.seed(20260812)
options(future.globals.maxSize = 20 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

project_root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

input_rds <- file.path(
  project_root,
  "04_preprocessing",
  "GSE288887_broad_clustered.rds"
)

out_dir <- file.path(
  project_root,
  "06_CD4T_analysis",
  "TNK_subcluster_v1"
)

fig_dir <- file.path(
  project_root,
  "09_figures",
  "TNK_subcluster_v1"
)

tab_dir <- file.path(
  project_root,
  "10_tables",
  "TNK_subcluster_v1"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

message("===== READ BROAD OBJECT =====")
stopifnot(file.exists(input_rds))
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

cluster_candidates <- c(
  "broad_cluster",
  "RNA_snn_res.0.4",
  "seurat_clusters"
)

cluster_col <- cluster_candidates[
  cluster_candidates %in% colnames(obj[[]])
][1]

if (length(cluster_col) == 0 || is.na(cluster_col)) {
  stop(
    "Cannot find broad cluster column. Available metadata columns:\n",
    paste(colnames(obj[[]]), collapse = ", ")
  )
}

obj$broad_cluster <- as.character(obj[[cluster_col, drop = TRUE]])

message("Broad cluster column: ", cluster_col)
print(table(obj$broad_cluster))

required_clusters <- c("4", "9")
missing_clusters <- setdiff(required_clusters, unique(obj$broad_cluster))

if (length(missing_clusters) > 0) {
  stop(
    "Required broad clusters are absent: ",
    paste(missing_clusters, collapse = ", ")
  )
}

message("===== EXTRACT BROAD CLUSTERS 4 AND 9 =====")
keep_cells <- rownames(obj[[]])[
  obj$broad_cluster %in% required_clusters
]

tnk <- subset(obj, cells = keep_cells)
rm(obj)
gc()

message("Candidate T/NK cells: ", ncol(tnk))
print(table(tnk$broad_cluster))

DefaultAssay(tnk) <- "RNA"

# Seurat v5对象在差异分析前需要合并RNA layers
if (inherits(tnk[["RNA"]], "Assay5")) {
  message("Joining RNA layers...")
  tnk <- JoinLayers(tnk, assay = "RNA")
}

message("===== REPROCESS T/NK COMPARTMENT =====")
tnk <- NormalizeData(
  tnk,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

tnk <- FindVariableFeatures(
  tnk,
  selection.method = "vst",
  nfeatures = 2500,
  verbose = TRUE
)

tnk <- ScaleData(
  tnk,
  features = VariableFeatures(tnk),
  verbose = TRUE
)

tnk <- RunPCA(
  tnk,
  features = VariableFeatures(tnk),
  npcs = 40,
  verbose = TRUE,
  seed.use = 20260812
)

if (!"sample_id" %in% colnames(tnk[[]])) {
  stop("sample_id is absent from metadata")
}

message("===== HARMONY BY LIBRARY/SAMPLE =====")
tnk <- harmony::RunHarmony(
  object = tnk,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  plot_convergence = FALSE,
  verbose = TRUE
)

tnk <- FindNeighbors(
  tnk,
  reduction = "harmony",
  dims = 1:30,
  k.param = 30,
  verbose = TRUE
)

tnk <- FindClusters(
  tnk,
  resolution = c(0.2, 0.4, 0.6, 0.8),
  random.seed = 20260812,
  verbose = TRUE
)

tnk <- RunUMAP(
  tnk,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.30,
  metric = "cosine",
  seed.use = 20260812,
  verbose = TRUE
)

res_col <- "RNA_snn_res.0.4"

if (!res_col %in% colnames(tnk[[]])) {
  stop("Expected clustering column is absent: ", res_col)
}

tnk$tnk_cluster <- as.character(tnk[[res_col, drop = TRUE]])
Idents(tnk) <- "tnk_cluster"

message("===== T/NK CLUSTER COUNTS =====")
print(table(Idents(tnk)))

output_rds <- file.path(
  out_dir,
  "GSE288887_TNK_subclustered.rds"
)

saveRDS(tnk, output_rds, compress = FALSE)
message("Saved: ", output_rds)

message("===== COMPOSITION TABLES =====")
meta <- tnk[[]]

sample_counts <- as.data.frame.matrix(
  table(meta$sample_id, meta$tnk_cluster)
)
sample_counts$sample_id <- rownames(sample_counts)
rownames(sample_counts) <- NULL
sample_counts <- sample_counts[
  ,
  c("sample_id", setdiff(colnames(sample_counts), "sample_id")),
  drop = FALSE
]

write.table(
  sample_counts,
  file.path(tab_dir, "TNK_cluster_counts_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sample_matrix <- as.matrix(
  sample_counts[, setdiff(colnames(sample_counts), "sample_id"), drop = FALSE]
)

sample_props <- sweep(
  sample_matrix,
  1,
  pmax(rowSums(sample_matrix), 1),
  "/"
)

sample_props <- as.data.frame(sample_props)
sample_props$sample_id <- sample_counts$sample_id
sample_props <- sample_props[
  ,
  c("sample_id", setdiff(colnames(sample_props), "sample_id")),
  drop = FALSE
]

write.table(
  sample_props,
  file.path(tab_dir, "TNK_cluster_proportions_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

group_candidates <- c("group", "condition", "disease_group")
group_col <- group_candidates[
  group_candidates %in% colnames(meta)
][1]

if (length(group_col) > 0 && !is.na(group_col)) {
  group_counts <- as.data.frame.matrix(
    table(meta[[group_col]], meta$tnk_cluster)
  )
  group_counts[[group_col]] <- rownames(group_counts)
  rownames(group_counts) <- NULL
  group_counts <- group_counts[
    ,
    c(group_col, setdiff(colnames(group_counts), group_col)),
    drop = FALSE
  ]

  write.table(
    group_counts,
    file.path(tab_dir, "TNK_cluster_counts_by_group.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

message("===== UMAP FIGURES =====")
p_cluster <- DimPlot(
  tnk,
  reduction = "umap",
  group.by = "tnk_cluster",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("T/NK clusters") +
  theme_classic(base_size = 13)

ggsave(
  file.path(fig_dir, "GSE288887_TNK_UMAP_clusters.png"),
  p_cluster,
  width = 10,
  height = 8,
  dpi = 240
)

p_source <- DimPlot(
  tnk,
  reduction = "umap",
  group.by = "broad_cluster",
  raster = TRUE
) +
  ggtitle("Broad source") +
  theme_classic(base_size = 13)

ggsave(
  file.path(fig_dir, "GSE288887_TNK_UMAP_broad_source.png"),
  p_source,
  width = 10,
  height = 8,
  dpi = 240
)

if (length(group_col) > 0 && !is.na(group_col)) {
  p_group <- DimPlot(
    tnk,
    reduction = "umap",
    group.by = group_col,
    split.by = group_col,
    raster = TRUE
  ) +
    ggtitle("Groups") +
    theme_classic(base_size = 13)

  ggsave(
    file.path(fig_dir, "GSE288887_TNK_UMAP_groups.png"),
    p_group,
    width = 14,
    height = 7,
    dpi = 240
  )
}

marker_sets <- list(
  PanT = c(
    "CD3D", "CD3E", "CD3G", "TRAC",
    "TRBC1", "TRBC2", "CD247", "LCK"
  ),
  CD4N = c(
    "CD4", "IL7R", "CCR7", "SELL",
    "TCF7", "LEF1", "LTB", "MAL"
  ),
  CD8C = c(
    "CD8A", "CD8B", "CCL5", "GZMK",
    "GZMB", "PRF1"
  ),
  NK = c(
    "KLRD1", "FCGR3A", "TYROBP",
    "NKG7", "GNLY", "TRDC"
  ),
  Treg = c(
    "FOXP3", "IL2RA", "CTLA4",
    "TIGIT", "IKZF2"
  ),
  Th17 = c(
    "CCR6", "KLRB1", "RORA",
    "IL23R", "IL17A"
  ),
  Act = c(
    "PDCD1", "TOX", "LAG3",
    "HAVCR2", "CXCL13", "HLA-DRA"
  ),
  Cycle = c(
    "MKI67", "TOP2A", "STMN1"
  ),
  Contam = c(
    "LST1", "FCER1G", "MS4A1",
    "CD79A", "EPCAM", "COL1A1"
  ),
  Target = c("MRE11")
)

marker_sets <- lapply(
  marker_sets,
  function(x) intersect(x, rownames(tnk))
)
marker_sets <- marker_sets[lengths(marker_sets) > 0]

p_dot <- DotPlot(
  tnk,
  features = marker_sets,
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle("T/NK markers") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  file.path(fig_dir, "GSE288887_TNK_marker_dotplot.png"),
  p_dot,
  width = 18,
  height = 8,
  dpi = 240
)

feature_genes <- intersect(
  c(
    "CD3D", "CD4", "IL7R", "CCR7",
    "CD8A", "NKG7", "GNLY", "KLRD1",
    "FOXP3", "IL2RA", "MKI67", "MRE11"
  ),
  rownames(tnk)
)

p_feature <- FeaturePlot(
  tnk,
  features = feature_genes,
  reduction = "umap",
  order = TRUE,
  raster = TRUE,
  ncol = 4
)

ggsave(
  file.path(fig_dir, "GSE288887_TNK_key_featureplots.png"),
  p_feature,
  width = 15,
  height = 12,
  dpi = 240
)

message("===== FIND SUBCLUSTER MARKERS =====")
Idents(tnk) <- "tnk_cluster"

markers <- FindAllMarkers(
  tnk,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10,
  max.cells.per.ident = 5000,
  random.seed = 20260812,
  verbose = TRUE
)

write.table(
  markers,
  file.path(tab_dir, "GSE288887_TNK_subcluster_markers.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("===== FINAL SUMMARY =====")
message("Input broad clusters: 4,9")
message("T/NK cells: ", ncol(tnk))
message("T/NK subclusters: ", length(unique(tnk$tnk_cluster)))
message("Object: ", output_rds)
message("STATUS=PASS")
