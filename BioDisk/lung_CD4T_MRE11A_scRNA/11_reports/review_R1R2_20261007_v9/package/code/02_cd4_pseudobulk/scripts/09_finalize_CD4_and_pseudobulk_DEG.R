suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

set.seed(20260812)

project_root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

input_rds <- file.path(
  project_root,
  "06_CD4T_analysis/CD4_subcluster_v1",
  "GSE288887_CD4_candidate_subclustered.rds"
)

object_dir <- file.path(
  project_root,
  "06_CD4T_analysis/CD4_final"
)

table_dir <- file.path(
  project_root,
  "10_tables/CD4_final"
)

deg_dir <- file.path(
  project_root,
  "07_DEG/CD4_pseudobulk_v1"
)

figure_dir <- file.path(
  project_root,
  "09_figures/CD4_final"
)

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

message("===== READ CD4 CANDIDATE OBJECT =====")
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

obj$cd4_cluster <- as.character(obj$cd4_cluster)

label_map <- c(
  "0" = "CD4-GZMA",
  "1" = "CD4-IL7R",
  "2" = "CD4-CCR7",
  "3" = "CD4state",
  "4" = "Treg",
  "5" = "MAITgd",
  "6" = "NKlike"
)

decision_map <- c(
  "0" = "Retain",
  "1" = "Retain",
  "2" = "Retain",
  "3" = "Retain",
  "4" = "Retain",
  "5" = "Exclude",
  "6" = "Exclude"
)

obj$cd4_label_all <- unname(
  label_map[obj$cd4_cluster]
)

obj$cd4_decision <- unname(
  decision_map[obj$cd4_cluster]
)

if (anyNA(obj$cd4_label_all) || anyNA(obj$cd4_decision)) {
  stop("Unmapped CD4 clusters detected")
}

annotation_table <- data.frame(
  cluster = names(label_map),
  annotation = unname(label_map),
  decision = unname(decision_map),
  stringsAsFactors = FALSE
)

annotation_table$cells <- as.integer(
  table(
    factor(
      obj$cd4_cluster,
      levels = annotation_table$cluster
    )
  )
)

write.table(
  annotation_table,
  file.path(table_dir, "CD4_cluster_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("===== FINAL CD4 EXTRACTION =====")
final_cells <- rownames(obj[[]])[
  obj$cd4_decision == "Retain"
]

excluded_cells <- rownames(obj[[]])[
  obj$cd4_decision == "Exclude"
]

excluded_meta <- obj[[]][
  excluded_cells,
  ,
  drop = FALSE
]

write.table(
  excluded_meta,
  file.path(table_dir, "CD4_excluded_cells_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

cd4 <- subset(obj, cells = final_cells)
rm(obj)
gc()

cd4$cd4_label <- factor(
  cd4$cd4_label_all,
  levels = c(
    "CD4-GZMA",
    "CD4-IL7R",
    "CD4-CCR7",
    "CD4state",
    "Treg"
  )
)

Idents(cd4) <- "cd4_label"

message("Final CD4 cells: ", ncol(cd4))
print(table(cd4$cd4_label))
print(table(cd4$group))

final_rds <- file.path(
  object_dir,
  "GSE288887_CD4_final.rds"
)

saveRDS(
  cd4,
  final_rds,
  compress = FALSE
)

message("Saved: ", final_rds)

message("===== ANNOTATED UMAP =====")
p_umap <- DimPlot(
  cd4,
  reduction = "umap",
  group.by = "cd4_label",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("CD4 anno") +
  labs(
    color = "Subtype",
    caption = paste(
      "CD4-GZMA: GZMA-positive CD4;",
      "CD4-IL7R: IL7R-positive CD4;",
      "CD4-CCR7: CCR7-positive naive/central-memory CD4;",
      "CD4state: CWP-enriched CD4 state;",
      "Treg: regulatory T cell."
    )
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(figure_dir, "GSE288887_CD4_final_UMAP.png"),
  p_umap,
  width = 11,
  height = 8,
  dpi = 260
)

message("===== METADATA AUDIT =====")
meta <- cd4[[]]

subject_candidates <- c(
  "subject_id",
  "patient_id",
  "donor_id"
)

subject_col <- subject_candidates[
  subject_candidates %in% colnames(meta)
][1]

group_candidates <- c(
  "group",
  "condition",
  "disease_group"
)

group_col <- group_candidates[
  group_candidates %in% colnames(meta)
][1]

if (is.na(subject_col)) {
  stop("No subject identifier found")
}

if (is.na(group_col)) {
  stop("No group variable found")
}

subject_group <- unique(
  meta[, c(subject_col, group_col), drop = FALSE]
)

colnames(subject_group) <- c(
  "subject_id",
  "group"
)

if (anyDuplicated(subject_group$subject_id)) {
  stop("A subject maps to more than one group")
}

subject_group$group <- factor(
  subject_group$group,
  levels = c("Control", "CWP")
)

if (anyNA(subject_group$group)) {
  stop("Unexpected group labels detected")
}

message("Subject column: ", subject_col)
message("Group column: ", group_col)
print(table(subject_group$group))

message("===== SUBJECT-LEVEL COMPOSITION =====")
composition_counts <- table(
  subject_id = meta[[subject_col]],
  subtype = meta$cd4_label
)

composition_props <- sweep(
  composition_counts,
  1,
  pmax(rowSums(composition_counts), 1),
  "/"
)

composition_long <- as.data.frame(
  composition_props,
  stringsAsFactors = FALSE
)

colnames(composition_long) <- c(
  "subject_id",
  "subtype",
  "proportion"
)

composition_long <- merge(
  composition_long,
  subject_group,
  by = "subject_id",
  all.x = TRUE,
  sort = FALSE
)

subject_totals <- data.frame(
  subject_id = rownames(composition_counts),
  final_cd4_cells = rowSums(composition_counts),
  stringsAsFactors = FALSE
)

subject_totals <- merge(
  subject_totals,
  subject_group,
  by = "subject_id",
  all.x = TRUE,
  sort = FALSE
)

write.table(
  subject_totals,
  file.path(table_dir, "CD4_cells_by_subject.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  composition_long,
  file.path(table_dir, "CD4_subtype_proportions_by_subject.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

subtypes <- levels(cd4$cd4_label)

composition_tests <- do.call(
  rbind,
  lapply(subtypes, function(subtype_i) {
    z <- composition_long[
      composition_long$subtype == subtype_i,
      ,
      drop = FALSE
    ]

    wt <- wilcox.test(
      proportion ~ group,
      data = z,
      exact = FALSE
    )

    data.frame(
      subtype = subtype_i,
      control_median = median(
        z$proportion[z$group == "Control"]
      ),
      cwp_median = median(
        z$proportion[z$group == "CWP"]
      ),
      delta_median = median(
        z$proportion[z$group == "CWP"]
      ) - median(
        z$proportion[z$group == "Control"]
      ),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
)

composition_tests$fdr_bh <- p.adjust(
  composition_tests$p_value,
  method = "BH"
)

write.table(
  composition_tests,
  file.path(table_dir, "CD4_subtype_composition_tests.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_comp <- ggplot(
  composition_long,
  aes(
    x = subtype,
    y = proportion,
    color = group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.7)
  ) +
  geom_jitter(
    aes(shape = group),
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "CD4 comp",
    x = "Subtype",
    y = "Prop",
    color = "Group",
    shape = "Group",
    caption = paste(
      "Comp: composition; Prop: subject-level proportion.",
      "Each point represents one subject."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

ggsave(
  file.path(figure_dir, "GSE288887_CD4_subject_composition.png"),
  p_comp,
  width = 10,
  height = 7,
  dpi = 260
)

message("===== RAW COUNT MATRIX =====")
DefaultAssay(cd4) <- "RNA"

if (inherits(cd4[["RNA"]], "Assay5")) {
  available_layers <- Layers(cd4[["RNA"]])

  if (!"counts" %in% available_layers) {
    cd4 <- JoinLayers(cd4, assay = "RNA")
  }
}

counts <- GetAssayData(
  cd4,
  assay = "RNA",
  layer = "counts"
)

meta <- cd4[[]][
  colnames(counts),
  ,
  drop = FALSE
]

aggregate_subject_counts <- function(
    count_matrix,
    metadata,
    subject_column) {

  subject_vector <- as.character(
    metadata[[subject_column]]
  )

  subjects <- sort(unique(subject_vector))

  aggregation_matrix <- sparseMatrix(
    i = seq_along(subject_vector),
    j = match(subject_vector, subjects),
    x = 1,
    dims = c(
      length(subject_vector),
      length(subjects)
    ),
    dimnames = list(
      colnames(count_matrix),
      subjects
    )
  )

  pseudobulk <- count_matrix %*% aggregation_matrix
  colnames(pseudobulk) <- subjects

  list(
    counts = pseudobulk,
    cell_counts = table(
      factor(
        subject_vector,
        levels = subjects
      )
    )
  )
}

run_pseudobulk <- function(
    count_matrix,
    metadata,
    analysis_name,
    output_prefix,
    min_cells = 20) {

  message("----- ", analysis_name, " -----")

  aggregation <- aggregate_subject_counts(
    count_matrix,
    metadata,
    subject_col
  )

  pb_counts <- aggregation$counts

  pb_info <- data.frame(
    subject_id = colnames(pb_counts),
    cell_count = as.integer(
      aggregation$cell_counts
    ),
    stringsAsFactors = FALSE
  )

  pb_info <- merge(
    pb_info,
    subject_group,
    by = "subject_id",
    all.x = TRUE,
    sort = FALSE
  )

  pb_info <- pb_info[
    match(
      colnames(pb_counts),
      pb_info$subject_id
    ),
    ,
    drop = FALSE
  ]

  eligible <- (
    pb_info$cell_count >= min_cells &
      !is.na(pb_info$group)
  )

  eligibility_table <- pb_info
  eligibility_table$eligible <- eligible

  write.table(
    eligibility_table,
    file.path(
      deg_dir,
      paste0(
        output_prefix,
        "_subject_eligibility.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  pb_info <- pb_info[eligible, , drop = FALSE]
  pb_counts <- pb_counts[
    ,
    pb_info$subject_id,
    drop = FALSE
  ]

  pb_info$group <- factor(
    pb_info$group,
    levels = c("Control", "CWP")
  )

  if (any(table(pb_info$group) < 3)) {
    stop(
      analysis_name,
      ": fewer than three eligible subjects in a group"
    )
  }

  message("Eligible subjects:")
  print(table(pb_info$group))

  design <- model.matrix(
    ~ group,
    data = pb_info
  )

  rownames(design) <- pb_info$subject_id

  y <- DGEList(
    counts = as.matrix(pb_counts)
  )

  keep_gene <- filterByExpr(
    y,
    design = design
  )

  gene_filter <- data.frame(
    gene = rownames(y),
    retained = keep_gene,
    stringsAsFactors = FALSE
  )

  write.table(
    gene_filter,
    file.path(
      deg_dir,
      paste0(
        output_prefix,
        "_gene_filter.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  y <- y[
    keep_gene,
    ,
    keep.lib.sizes = FALSE
  ]

  y <- calcNormFactors(
    y,
    method = "TMM"
  )

  y <- estimateDisp(
    y,
    design,
    robust = TRUE
  )

  fit <- glmQLFit(
    y,
    design,
    robust = TRUE
  )

  coefficient <- grep(
    "^groupCWP$",
    colnames(design)
  )

  if (length(coefficient) != 1) {
    stop("Cannot identify groupCWP coefficient")
  }

  test <- glmQLFTest(
    fit,
    coef = coefficient
  )

  result <- topTags(
    test,
    n = Inf,
    sort.by = "PValue"
  )$table

  result$gene <- rownames(result)
  rownames(result) <- NULL

  result$direction <- "NS"

  result$direction[
    result$FDR < 0.05 &
      result$logFC >= 0.5
  ] <- "UpCWP"

  result$direction[
    result$FDR < 0.05 &
      result$logFC <= -0.5
  ] <- "DownCWP"

  result <- result[
    ,
    c(
      "gene",
      "logFC",
      "logCPM",
      "F",
      "PValue",
      "FDR",
      "direction"
    )
  ]

  result_file <- file.path(
    deg_dir,
    paste0(
      output_prefix,
      "_edgeR_QLF.tsv"
    )
  )

  write.table(
    result,
    result_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  log_cpm <- cpm(
    y,
    log = TRUE,
    prior.count = 2
  )

  gene_variance <- apply(
    log_cpm,
    1,
    var
  )

  top_variable <- names(
    sort(
      gene_variance,
      decreasing = TRUE
    )
  )

  top_variable <- head(
    top_variable,
    min(1000, length(top_variable))
  )

  pca <- prcomp(
    t(
      log_cpm[
        top_variable,
        ,
        drop = FALSE
      ]
    ),
    center = TRUE,
    scale. = FALSE
  )

  explained <- round(
    100 * summary(pca)$importance[2, 1:2],
    1
  )

  pca_df <- data.frame(
    subject = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = pb_info$group,
    stringsAsFactors = FALSE
  )

  p_pca <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      color = group,
      label = subject
    )
  ) +
    geom_point(size = 3) +
    geom_text(
      vjust = -0.7,
      size = 3,
      check_overlap = TRUE
    ) +
    labs(
      title = paste("PB PCA", analysis_name),
      x = paste0("PC1 ", explained[1], "%"),
      y = paste0("PC2 ", explained[2], "%"),
      color = "Group",
      caption = "PB: subject-level pseudobulk."
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(
      figure_dir,
      paste0(
        output_prefix,
        "_PB_PCA.png"
      )
    ),
    p_pca,
    width = 9,
    height = 7,
    dpi = 260
  )

  result$neglog_fdr <- -log10(
    pmax(result$FDR, 1e-50)
  )

  p_volcano <- ggplot(
    result,
    aes(
      x = logFC,
      y = neglog_fdr,
      color = direction
    )
  ) +
    geom_point(
      alpha = 0.65,
      size = 1.4
    ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      linetype = 2,
      color = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2,
      color = "grey50"
    ) +
    scale_color_manual(
      values = c(
        "DownCWP" = "#2166AC",
        "NS" = "grey75",
        "UpCWP" = "#B2182B"
      )
    ) +
    labs(
      title = paste("PB DEG", analysis_name),
      x = "log2FC",
      y = "-log10FDR",
      color = "Class",
      caption = paste(
        "PB: subject-level pseudobulk;",
        "positive log2FC indicates higher expression in CWP."
      )
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(
      figure_dir,
      paste0(
        output_prefix,
        "_PB_volcano.png"
      )
    ),
    p_volcano,
    width = 8,
    height = 7,
    dpi = 260
  )

  if ("MRE11" %in% rownames(pb_counts)) {
    mre11_table <- data.frame(
      subject_id = colnames(pb_counts),
      group = pb_info$group[
        match(
          colnames(pb_counts),
          pb_info$subject_id
        )
      ],
      cell_count = as.integer(
        aggregation$cell_counts[
          colnames(pb_counts)
        ]
      ),
      MRE11_count = as.numeric(
        pb_counts[
          "MRE11",
          ,
          drop = TRUE
        ]
      ),
      library_size = colSums(pb_counts),
      stringsAsFactors = FALSE
    )

    mre11_table$MRE11_CPM <- (
      1e6 * mre11_table$MRE11_count /
        pmax(mre11_table$library_size, 1)
    )

    write.table(
      mre11_table,
      file.path(
        deg_dir,
        paste0(
          output_prefix,
          "_MRE11_by_subject.tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  summary_table <- data.frame(
    analysis = analysis_name,
    control_subjects = sum(
      pb_info$group == "Control"
    ),
    cwp_subjects = sum(
      pb_info$group == "CWP"
    ),
    retained_genes = nrow(result),
    fdr_lt_005 = sum(
      result$FDR < 0.05,
      na.rm = TRUE
    ),
    up_cwp = sum(
      result$direction == "UpCWP"
    ),
    down_cwp = sum(
      result$direction == "DownCWP"
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    summary_table,
    file.path(
      deg_dir,
      paste0(
        output_prefix,
        "_summary.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Result: ", result_file)
  print(summary_table)

  invisible(summary_table)
}

message("===== PRIMARY PSEUDOBULK =====")
primary_summary <- run_pseudobulk(
  count_matrix = counts,
  metadata = meta,
  analysis_name = "Primary",
  output_prefix = "CD4_primary",
  min_cells = 20
)

message("===== CORE SENSITIVITY PSEUDOBULK =====")
core_cells <- meta$cd4_label != "CD4state"

core_summary <- run_pseudobulk(
  count_matrix = counts[, core_cells, drop = FALSE],
  metadata = meta[core_cells, , drop = FALSE],
  analysis_name = "Core",
  output_prefix = "CD4_core",
  min_cells = 20
)

combined_summary <- rbind(
  primary_summary,
  core_summary
)

write.table(
  combined_summary,
  file.path(
    deg_dir,
    "CD4_pseudobulk_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("===== FINAL OUTPUT =====")
message("Final CD4 cells: ", ncol(cd4))
message("Final object: ", final_rds)
message("DEG directory: ", deg_dir)
message("Figure directory: ", figure_dir)
message("STATUS=PASS")
