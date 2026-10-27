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
  "06_CD4T_analysis/CD4_final",
  "GSE288887_CD4_final.rds"
)

out_dir <- file.path(
  project_root,
  "07_DEG/CD4_pseudobulk_v1/source_adjusted"
)

figure_dir <- file.path(
  project_root,
  "09_figures/CD4_final/source_adjusted"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

message("===== READ FINAL CD4 OBJECT =====")
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

meta <- obj[[]]

subject_candidates <- c(
  "subject_id",
  "patient_id",
  "donor_id"
)

group_candidates <- c(
  "group",
  "condition",
  "disease_group"
)

subject_col <- subject_candidates[
  subject_candidates %in% colnames(meta)
][1]

group_col <- group_candidates[
  group_candidates %in% colnames(meta)
][1]

if (is.na(subject_col) || is.na(group_col)) {
  stop("Subject or group column is absent")
}

meta$subject_id2 <- as.character(
  meta[[subject_col]]
)

meta$group2 <- as.character(
  meta[[group_col]]
)

source_map <- c(
  "CWP85" = "G227",
  "CWP86" = "G227",

  "CTL092" = "G227",
  "CTL094" = "G227",
  "CTL098" = "G227",
  "CTL103" = "G227",
  "CTL105" = "G227",
  "CTL106" = "G227",
  "CTL108" = "G227",
  "CTL111" = "G227",
  "CTL117" = "G227",

  "CTL068" = "G135",
  "CTL070" = "G135",

  "CWP89" = "G288",
  "CWP94" = "G288",
  "CWP109" = "G288",
  "CWP127" = "G288",
  "CWP144" = "G288",
  "CWP149" = "G288",
  "CWP152" = "G288"
)

meta$series <- unname(
  source_map[meta$subject_id2]
)

if (anyNA(meta$series)) {
  stop(
    "Unmapped subjects: ",
    paste(
      unique(
        meta$subject_id2[
          is.na(meta$series)
        ]
      ),
      collapse = ","
    )
  )
}

message("===== SUBJECT DESIGN =====")
subject_info <- unique(
  meta[, c(
    "subject_id2",
    "group2",
    "series"
  )]
)

colnames(subject_info) <- c(
  "subject_id",
  "group",
  "series"
)

if (anyDuplicated(subject_info$subject_id)) {
  stop("A subject maps to multiple group/series combinations")
}

subject_info$group <- factor(
  subject_info$group,
  levels = c("Control", "CWP")
)

subject_info$series <- factor(
  subject_info$series,
  levels = c(
    "G227",
    "G288",
    "G135"
  )
)

print(
  addmargins(
    table(
      subject_info$series,
      subject_info$group
    )
  )
)

write.table(
  subject_info,
  file.path(
    out_dir,
    "subject_series_mapping.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("===== RAW COUNTS =====")
if (inherits(obj[["RNA"]], "Assay5")) {
  if (!"counts" %in% Layers(obj[["RNA"]])) {
    obj <- JoinLayers(obj, assay = "RNA")
  }
}

counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

meta <- meta[
  colnames(counts),
  ,
  drop = FALSE
]

build_pseudobulk <- function(
    cell_selector,
    analysis_name,
    min_cells = 20) {

  selected_counts <- counts[
    ,
    cell_selector,
    drop = FALSE
  ]

  selected_meta <- meta[
    cell_selector,
    ,
    drop = FALSE
  ]

  subjects <- sort(
    unique(selected_meta$subject_id2)
  )

  aggregation_matrix <- sparseMatrix(
    i = seq_len(ncol(selected_counts)),
    j = match(
      selected_meta$subject_id2,
      subjects
    ),
    x = 1,
    dims = c(
      ncol(selected_counts),
      length(subjects)
    ),
    dimnames = list(
      colnames(selected_counts),
      subjects
    )
  )

  pb_counts <- selected_counts %*% aggregation_matrix
  colnames(pb_counts) <- subjects

  cell_counts <- table(
    factor(
      selected_meta$subject_id2,
      levels = subjects
    )
  )

  info <- subject_info[
    match(
      subjects,
      subject_info$subject_id
    ),
    ,
    drop = FALSE
  ]

  info$cell_count <- as.integer(
    cell_counts
  )

  info$eligible <- (
    info$cell_count >= min_cells
  )

  write.table(
    info,
    file.path(
      out_dir,
      paste0(
        analysis_name,
        "_eligibility.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  info <- info[
    info$eligible,
    ,
    drop = FALSE
  ]

  pb_counts <- pb_counts[
    ,
    info$subject_id,
    drop = FALSE
  ]

  info$group <- droplevels(info$group)
  info$series <- droplevels(info$series)

  list(
    counts = pb_counts,
    info = info
  )
}

primary_pb <- build_pseudobulk(
  rep(TRUE, ncol(counts)),
  "Primary"
)

core_selector <- (
  as.character(meta$cd4_label) !=
    "CD4state"
)

core_pb <- build_pseudobulk(
  core_selector,
  "Core"
)

fit_edger <- function(
    pb_counts,
    info,
    design,
    analysis_name,
    model_name) {

  if (qr(design)$rank != ncol(design)) {
    stop(
      analysis_name,
      " ",
      model_name,
      ": design matrix is not full rank"
    )
  }

  rownames(design) <- info$subject_id

  y <- DGEList(
    counts = as.matrix(pb_counts)
  )

  keep_gene <- filterByExpr(
    y,
    design = design
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
    stop(
      "Cannot identify groupCWP coefficient in ",
      model_name
    )
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

  write.table(
    result,
    file.path(
      out_dir,
      paste0(
        analysis_name,
        "_",
        model_name,
        "_edgeR.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  result
}

run_analysis <- function(
    pb_object,
    analysis_name) {

  info <- pb_object$info
  pb_counts <- pb_object$counts

  message("===== ", analysis_name, " DESIGN =====")
  print(
    addmargins(
      table(
        info$series,
        info$group
      )
    )
  )

  design_unadjusted <- model.matrix(
    ~ group,
    data = info
  )

  design_adjusted <- model.matrix(
    ~ series + group,
    data = info
  )

  design_audit <- data.frame(
    model = c(
      "Unadj",
      "SrcAdj"
    ),
    samples = nrow(info),
    columns = c(
      ncol(design_unadjusted),
      ncol(design_adjusted)
    ),
    rank = c(
      qr(design_unadjusted)$rank,
      qr(design_adjusted)$rank
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    design_audit,
    file.path(
      out_dir,
      paste0(
        analysis_name,
        "_design_audit.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  unadjusted <- fit_edger(
    pb_counts,
    info,
    design_unadjusted,
    analysis_name,
    "Unadj"
  )

  adjusted <- fit_edger(
    pb_counts,
    info,
    design_adjusted,
    analysis_name,
    "SrcAdj"
  )

  bridge_selector <- (
    info$series == "G227"
  )

  bridge_info <- droplevels(
    info[
      bridge_selector,
      ,
      drop = FALSE
    ]
  )

  bridge_counts <- pb_counts[
    ,
    bridge_info$subject_id,
    drop = FALSE
  ]

  if (
    sum(bridge_info$group == "CWP") < 2 ||
      sum(bridge_info$group == "Control") < 3
  ) {
    stop(
      analysis_name,
      ": insufficient G227 bridge subjects"
    )
  }

  design_bridge <- model.matrix(
    ~ group,
    data = bridge_info
  )

  bridge <- fit_edger(
    bridge_counts,
    bridge_info,
    design_bridge,
    analysis_name,
    "G227"
  )

  adjusted_small <- adjusted[
    ,
    c(
      "gene",
      "logFC",
      "FDR",
      "direction"
    )
  ]

  colnames(adjusted_small)[2:4] <- c(
    "logFC_SrcAdj",
    "FDR_SrcAdj",
    "direction_SrcAdj"
  )

  unadjusted_small <- unadjusted[
    ,
    c(
      "gene",
      "logFC",
      "FDR",
      "direction"
    )
  ]

  colnames(unadjusted_small)[2:4] <- c(
    "logFC_Unadj",
    "FDR_Unadj",
    "direction_Unadj"
  )

  bridge_small <- bridge[
    ,
    c(
      "gene",
      "logFC",
      "FDR",
      "direction"
    )
  ]

  colnames(bridge_small)[2:4] <- c(
    "logFC_G227",
    "FDR_G227",
    "direction_G227"
  )

  comparison <- merge(
    adjusted_small,
    unadjusted_small,
    by = "gene",
    all = TRUE,
    sort = FALSE
  )

  comparison <- merge(
    comparison,
    bridge_small,
    by = "gene",
    all = TRUE,
    sort = FALSE
  )

  comparison$adjusted_significant <- (
    !is.na(comparison$FDR_SrcAdj) &
      comparison$FDR_SrcAdj < 0.05 &
      abs(comparison$logFC_SrcAdj) >= 0.5
  )

  comparison$direction_concordant <- (
    !is.na(comparison$logFC_SrcAdj) &
      !is.na(comparison$logFC_Unadj) &
      !is.na(comparison$logFC_G227) &
      sign(comparison$logFC_SrcAdj) ==
        sign(comparison$logFC_Unadj) &
      sign(comparison$logFC_SrcAdj) ==
        sign(comparison$logFC_G227)
  )

  comparison$final_candidate <- (
    comparison$adjusted_significant &
      comparison$direction_concordant
  )

  comparison <- comparison[
    order(
      -comparison$final_candidate,
      comparison$FDR_SrcAdj
    ),
    ,
    drop = FALSE
  ]

  write.table(
    comparison,
    file.path(
      out_dir,
      paste0(
        analysis_name,
        "_source_adjustment_comparison.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  summary_table <- data.frame(
    analysis = analysis_name,
    eligible_control = sum(
      info$group == "Control"
    ),
    eligible_cwp = sum(
      info$group == "CWP"
    ),
    unadjusted_deg = sum(
      unadjusted$FDR < 0.05 &
        abs(unadjusted$logFC) >= 0.5
    ),
    source_adjusted_deg = sum(
      adjusted$FDR < 0.05 &
        abs(adjusted$logFC) >= 0.5
    ),
    g227_deg = sum(
      bridge$FDR < 0.05 &
        abs(bridge$logFC) >= 0.5
    ),
    final_candidates = sum(
      comparison$final_candidate,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )

  if ("MRE11" %in% comparison$gene) {
    mre11_row <- comparison[
      comparison$gene == "MRE11",
      ,
      drop = FALSE
    ]

    write.table(
      mre11_row,
      file.path(
        out_dir,
        paste0(
          analysis_name,
          "_MRE11_DEG_record.tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  list(
    summary = summary_table,
    comparison = comparison,
    adjusted = adjusted
  )
}

primary_result <- run_analysis(
  primary_pb,
  "Primary"
)

core_result <- run_analysis(
  core_pb,
  "Core"
)

combined_summary <- rbind(
  primary_result$summary,
  core_result$summary
)

write.table(
  combined_summary,
  file.path(
    out_dir,
    "source_adjusted_DEG_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("===== PCA WITH SOURCE SERIES =====")
make_source_pca <- function(
    pb_object,
    analysis_name) {

  info <- pb_object$info
  pb_counts <- pb_object$counts

  design <- model.matrix(
    ~ series + group,
    data = info
  )

  y <- DGEList(
    counts = as.matrix(pb_counts)
  )

  keep_gene <- filterByExpr(
    y,
    design = design
  )

  y <- y[
    keep_gene,
    ,
    keep.lib.sizes = FALSE
  ]

  y <- calcNormFactors(y)

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
    group = info$group,
    series = info$series,
    stringsAsFactors = FALSE
  )

  p <- ggplot(
    pca_df,
    aes(
      PC1,
      PC2,
      color = group,
      shape = series,
      label = subject
    )
  ) +
    geom_point(size = 3.5) +
    geom_text(
      vjust = -0.7,
      size = 3,
      check_overlap = TRUE
    ) +
    labs(
      title = paste("PB src", analysis_name),
      x = paste0(
        "PC1 ",
        explained[1],
        "%"
      ),
      y = paste0(
        "PC2 ",
        explained[2],
        "%"
      ),
      color = "Group",
      shape = "Series",
      caption = paste(
        "PB: subject-level pseudobulk;",
        "G288: new GSE288887;",
        "G227: GSE227136;",
        "G135: GSE135893."
      )
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(
      figure_dir,
      paste0(
        analysis_name,
        "_PB_source_PCA.png"
      )
    ),
    p,
    width = 10,
    height = 7.5,
    dpi = 260
  )
}

make_source_pca(
  primary_pb,
  "Primary"
)

make_source_pca(
  core_pb,
  "Core"
)

message("===== SUMMARY =====")
print(combined_summary, row.names = FALSE)

message("===== CORE FINAL CANDIDATES =====")
core_candidates <- core_result$comparison[
  core_result$comparison$final_candidate,
  c(
    "gene",
    "logFC_SrcAdj",
    "FDR_SrcAdj",
    "logFC_Unadj",
    "logFC_G227"
  ),
  drop = FALSE
]

print(
  core_candidates,
  row.names = FALSE
)

message("===== FINAL OUTPUT =====")
message("Output: ", out_dir)
message("Figures: ", figure_dir)
message("STATUS=PASS")
