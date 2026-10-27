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
  "07_DEG/CD4_pseudobulk_v1/robustness"
)

table_dir <- file.path(
  project_root,
  "10_tables/CD4_final"
)

figure_dir <- file.path(
  project_root,
  "09_figures/CD4_final/robustness"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
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

meta$subject_audit <- as.character(
  meta[[subject_col]]
)

meta$group_audit <- as.character(
  meta[[group_col]]
)

if (
  anyNA(meta$subject_audit) ||
  anyNA(meta$group_audit)
) {
  stop("Actual cells contain missing subject/group values")
}

observed_subjects <- sort(
  unique(meta$subject_audit)
)

subject_group <- unique(
  meta[, c(
    "subject_audit",
    "group_audit"
  )]
)

colnames(subject_group) <- c(
  "subject_id",
  "group"
)

if (anyDuplicated(subject_group$subject_id)) {
  stop("A subject maps to multiple groups")
}

subject_group <- subject_group[
  match(
    observed_subjects,
    subject_group$subject_id
  ),
  ,
  drop = FALSE
]

subject_group$group <- factor(
  subject_group$group,
  levels = c("Control", "CWP")
)

stopifnot(!anyNA(subject_group$group))

message("Observed subjects:")
print(table(subject_group$group))

message("===== CORRECT COMPOSITION STATISTICS =====")
subtype_levels <- c(
  "CD4-GZMA",
  "CD4-IL7R",
  "CD4-CCR7",
  "CD4state",
  "Treg"
)

meta$cd4_label <- factor(
  as.character(meta$cd4_label),
  levels = subtype_levels
)

composition_counts <- table(
  factor(
    meta$subject_audit,
    levels = observed_subjects
  ),
  meta$cd4_label
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

composition_long$subject_id <- as.character(
  composition_long$subject_id
)

composition_long$subtype <- as.character(
  composition_long$subtype
)

group_map <- setNames(
  as.character(subject_group$group),
  subject_group$subject_id
)

composition_long$group <- unname(
  group_map[composition_long$subject_id]
)

if (anyNA(composition_long$group)) {
  stop("Composition table contains unmapped subjects")
}

composition_long$group <- factor(
  composition_long$group,
  levels = c("Control", "CWP")
)

composition_tests <- do.call(
  rbind,
  lapply(subtype_levels, function(subtype_i) {
    z <- composition_long[
      composition_long$subtype == subtype_i,
      ,
      drop = FALSE
    ]

    control_values <- z$proportion[
      z$group == "Control"
    ]

    cwp_values <- z$proportion[
      z$group == "CWP"
    ]

    wt <- wilcox.test(
      control_values,
      cwp_values,
      exact = FALSE
    )

    data.frame(
      subtype = subtype_i,
      control_median = median(
        control_values,
        na.rm = TRUE
      ),
      cwp_median = median(
        cwp_values,
        na.rm = TRUE
      ),
      delta_median = median(
        cwp_values,
        na.rm = TRUE
      ) - median(
        control_values,
        na.rm = TRUE
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
  file.path(
    table_dir,
    "CD4_subtype_composition_tests_corrected.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  composition_long,
  file.path(
    table_dir,
    "CD4_subtype_proportions_corrected.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_comp <- ggplot(
  composition_long,
  aes(
    x = subtype,
    y = 100 * proportion,
    color = group,
    shape = group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.65)
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.2,
    alpha = 0.85
  ) +
  labs(
    title = "CD4 comp",
    x = "Subtype",
    y = "Prop (%)",
    color = "Group",
    shape = "Group",
    caption = paste(
      "Comp: composition; Prop: subject-level proportion.",
      "Each point represents one observed subject."
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
  file.path(
    figure_dir,
    "CD4_subject_composition_corrected.png"
  ),
  p_comp,
  width = 10,
  height = 7,
  dpi = 260
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

aggregate_counts <- function(
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

  subject_vector <- selected_meta$subject_audit
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
      colnames(selected_counts),
      subjects
    )
  )

  pb_counts <- selected_counts %*% aggregation_matrix
  colnames(pb_counts) <- subjects

  cell_counts <- table(
    factor(
      subject_vector,
      levels = subjects
    )
  )

  info <- data.frame(
    subject_id = subjects,
    group = unname(group_map[subjects]),
    cell_count = as.integer(cell_counts),
    library_size = as.numeric(
      colSums(pb_counts)
    ),
    detected_genes = as.numeric(
      colSums(pb_counts > 0)
    ),
    stringsAsFactors = FALSE
  )

  info$group <- factor(
    info$group,
    levels = c("Control", "CWP")
  )

  info$eligible <- (
    info$cell_count >= min_cells
  )

  if ("MRE11" %in% rownames(pb_counts)) {
    info$MRE11_count <- as.numeric(
      pb_counts["MRE11", ]
    )

    info$MRE11_CPM <- (
      1e6 * info$MRE11_count /
        pmax(info$library_size, 1)
    )
  }

  subtype_tab <- table(
    factor(
      selected_meta$subject_audit,
      levels = subjects
    ),
    factor(
      selected_meta$cd4_label,
      levels = subtype_levels
    )
  )

  subtype_prop <- sweep(
    subtype_tab,
    1,
    pmax(rowSums(subtype_tab), 1),
    "/"
  )

  subtype_prop <- as.data.frame.matrix(
    subtype_prop
  )

  colnames(subtype_prop) <- paste0(
    "prop_",
    colnames(subtype_prop)
  )

  info <- cbind(
    info,
    subtype_prop[
      info$subject_id,
      ,
      drop = FALSE
    ]
  )

  write.table(
    info,
    file.path(
      out_dir,
      paste0(
        analysis_name,
        "_subject_QC.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  eligible_subjects <- info$subject_id[
    info$eligible
  ]

  list(
    counts = pb_counts[
      ,
      eligible_subjects,
      drop = FALSE
    ],
    info = info[
      info$eligible,
      ,
      drop = FALSE
    ],
    all_info = info
  )
}

primary_pb <- aggregate_counts(
  rep(TRUE, ncol(counts)),
  "Primary"
)

core_pb <- aggregate_counts(
  as.character(meta$cd4_label) != "CD4state",
  "Core"
)

message("===== ELIGIBILITY =====")
cat("Primary:\n")
print(
  primary_pb$all_info[
    ,
    c(
      "subject_id",
      "group",
      "cell_count",
      "eligible"
    )
  ],
  row.names = FALSE
)

cat("\nCore:\n")
print(
  core_pb$all_info[
    ,
    c(
      "subject_id",
      "group",
      "cell_count",
      "eligible"
    )
  ],
  row.names = FALSE
)

run_edger <- function(
    pb_object,
    analysis_name,
    model_name,
    excluded_subjects = character(0),
    make_pca = FALSE) {

  keep_subject <- !(
    pb_object$info$subject_id %in%
      excluded_subjects
  )

  info <- pb_object$info[
    keep_subject,
    ,
    drop = FALSE
  ]

  pb_counts <- pb_object$counts[
    ,
    info$subject_id,
    drop = FALSE
  ]

  info$group <- factor(
    info$group,
    levels = c("Control", "CWP")
  )

  if (any(table(info$group) < 3)) {
    stop(
      analysis_name,
      " ",
      model_name,
      ": insufficient subjects"
    )
  }

  design <- model.matrix(
    ~ group,
    data = info
  )

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

  test <- glmQLFTest(
    fit,
    coef = "groupCWP"
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

  prefix <- paste0(
    analysis_name,
    "_",
    model_name
  )

  write.table(
    result,
    file.path(
      out_dir,
      paste0(prefix, "_edgeR.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (make_pca) {
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
      stringsAsFactors = FALSE
    )

    write.table(
      pca_df,
      file.path(
        out_dir,
        paste0(prefix, "_PCA_scores.tsv")
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    p <- ggplot(
      pca_df,
      aes(
        PC1,
        PC2,
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
        title = paste("PB PCA", model_name),
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
        caption = paste(
          "PB: subject-level pseudobulk;",
          paste(
            excluded_subjects,
            collapse = ","
          ),
          "excluded."
        )
      ) +
      theme_classic(base_size = 12)

    ggsave(
      file.path(
        figure_dir,
        paste0(prefix, "_PCA.png")
      ),
      p,
      width = 9,
      height = 7,
      dpi = 260
    )
  }

  summary_row <- data.frame(
    analysis = analysis_name,
    model = model_name,
    control_n = sum(
      info$group == "Control"
    ),
    cwp_n = sum(
      info$group == "CWP"
    ),
    retained_genes = nrow(result),
    fdr_lt_005 = sum(
      result$FDR < 0.05
    ),
    up_cwp = sum(
      result$direction == "UpCWP"
    ),
    down_cwp = sum(
      result$direction == "DownCWP"
    ),
    stringsAsFactors = FALSE
  )

  list(
    result = result,
    summary = summary_row
  )
}

model_definitions <- list(
  All = character(0),
  No068 = "CTL068",
  No106 = "CTL106",
  NoBoth = c(
    "CTL068",
    "CTL106"
  )
)

all_results <- list()
all_summaries <- list()

for (analysis_name in c("Primary", "Core")) {
  pb_object <- if (
    analysis_name == "Primary"
  ) {
    primary_pb
  } else {
    core_pb
  }

  for (model_name in names(model_definitions)) {
    key <- paste(
      analysis_name,
      model_name,
      sep = "__"
    )

    model_result <- run_edger(
      pb_object = pb_object,
      analysis_name = analysis_name,
      model_name = model_name,
      excluded_subjects =
        model_definitions[[model_name]],
      make_pca = model_name == "NoBoth"
    )

    all_results[[key]] <- model_result$result
    all_summaries[[key]] <- model_result$summary
  }
}

summary_table <- do.call(
  rbind,
  all_summaries
)

write.table(
  summary_table,
  file.path(
    out_dir,
    "DEG_sensitivity_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

build_stability <- function(analysis_name) {
  model_names <- names(model_definitions)

  base_result <- all_results[[
    paste(
      analysis_name,
      "All",
      sep = "__"
    )
  ]]

  base_significant <- base_result[
    base_result$FDR < 0.05 &
      abs(base_result$logFC) >= 0.5,
    c(
      "gene",
      "logFC",
      "FDR",
      "direction"
    ),
    drop = FALSE
  ]

  colnames(base_significant)[2:4] <- c(
    "logFC_All",
    "FDR_All",
    "direction_All"
  )

  stability <- base_significant

  for (model_name in setdiff(
    model_names,
    "All"
  )) {
    z <- all_results[[
      paste(
        analysis_name,
        model_name,
        sep = "__"
      )
    ]]

    z <- z[
      ,
      c(
        "gene",
        "logFC",
        "FDR",
        "direction"
      ),
      drop = FALSE
    ]

    colnames(z)[2:4] <- c(
      paste0("logFC_", model_name),
      paste0("FDR_", model_name),
      paste0("direction_", model_name)
    )

    stability <- merge(
      stability,
      z,
      by = "gene",
      all.x = TRUE,
      sort = FALSE
    )
  }

  logfc_columns <- grep(
    "^logFC_",
    colnames(stability),
    value = TRUE
  )

  fdr_columns <- grep(
    "^FDR_",
    colnames(stability),
    value = TRUE
  )

  logfc_matrix <- as.matrix(
    stability[, logfc_columns, drop = FALSE]
  )

  fdr_matrix <- as.matrix(
    stability[, fdr_columns, drop = FALSE]
  )

  base_sign <- sign(
    stability$logFC_All
  )

  stability$models_present <- rowSums(
    !is.na(logfc_matrix)
  )

  stability$sign_consistent_n <- rowSums(
    sweep(
      sign(logfc_matrix),
      1,
      base_sign,
      "=="
    ),
    na.rm = TRUE
  )

  stability$sig_model_n <- rowSums(
    !is.na(fdr_matrix) &
      fdr_matrix < 0.05 &
      abs(logfc_matrix) >= 0.5,
    na.rm = TRUE
  )

  stability$stable_all4 <- (
    stability$models_present == 4 &
      stability$sign_consistent_n == 4 &
      stability$sig_model_n == 4
  )

  stability$max_fdr <- apply(
    fdr_matrix,
    1,
    max,
    na.rm = TRUE
  )

  stability <- stability[
    order(
      -stability$stable_all4,
      stability$max_fdr
    ),
    ,
    drop = FALSE
  ]

  write.table(
    stability,
    file.path(
      out_dir,
      paste0(
        analysis_name,
        "_DEG_stability.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  stability
}

primary_stability <- build_stability(
  "Primary"
)

core_stability <- build_stability(
  "Core"
)

message("===== DEG SENSITIVITY SUMMARY =====")
print(summary_table, row.names = FALSE)

message("===== PRIMARY STABLE GENES =====")
print(
  primary_stability[
    primary_stability$stable_all4,
    c(
      "gene",
      "logFC_All",
      "FDR_All",
      "max_fdr"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

message("===== CORE STABLE GENES =====")
print(
  core_stability[
    core_stability$stable_all4,
    c(
      "gene",
      "logFC_All",
      "FDR_All",
      "max_fdr"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

message("===== CORRECTED COMPOSITION =====")
print(composition_tests, row.names = FALSE)

message("===== FINAL OUTPUT =====")
message("Output: ", out_dir)
message("Figures: ", figure_dir)
message("STATUS=PASS")
