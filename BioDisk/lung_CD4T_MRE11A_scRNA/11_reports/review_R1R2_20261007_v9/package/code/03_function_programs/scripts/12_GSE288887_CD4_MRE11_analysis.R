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
  "08_MRE11A_analysis/GSE288887_CD4"
)

figure_dir <- file.path(
  project_root,
  "09_figures/MRE11_GSE288887"
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
  stop("Subject or group metadata is absent")
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
  stop("Some subjects could not be mapped to source series")
}

meta$group2 <- factor(
  meta$group2,
  levels = c("Control", "CWP")
)

meta$series <- factor(
  meta$series,
  levels = c("G227", "G288", "G135")
)

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

if (!"MRE11" %in% rownames(counts)) {
  stop("Official symbol MRE11 is absent")
}

meta <- meta[
  colnames(counts),
  ,
  drop = FALSE
]

message("===== MRE11 DETECTION =====")
mre11_detected <- as.numeric(
  counts["MRE11", ] > 0
)

mre11_umi <- as.numeric(
  counts["MRE11", ]
)

run_mre11_analysis <- function(
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

  selected_detected <- mre11_detected[
    cell_selector
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

  cell_n <- table(
    factor(
      selected_meta$subject_id2,
      levels = subjects
    )
  )

  detected_n <- tapply(
    selected_detected,
    factor(
      selected_meta$subject_id2,
      levels = subjects
    ),
    sum
  )

  subject_group <- tapply(
    as.character(selected_meta$group2),
    factor(
      selected_meta$subject_id2,
      levels = subjects
    ),
    function(x) unique(x)[1]
  )

  subject_series <- tapply(
    as.character(selected_meta$series),
    factor(
      selected_meta$subject_id2,
      levels = subjects
    ),
    function(x) unique(x)[1]
  )

  subject_table <- data.frame(
    subject_id = subjects,
    group = unname(subject_group),
    series = unname(subject_series),
    cell_count = as.integer(cell_n),
    detected_cells = as.integer(detected_n),
    MRE11_count = as.numeric(
      pb_counts["MRE11", ]
    ),
    library_size = as.numeric(
      colSums(pb_counts)
    ),
    stringsAsFactors = FALSE
  )

  subject_table$detection_pct <- (
    100 * subject_table$detected_cells /
      pmax(subject_table$cell_count, 1)
  )

  subject_table$MRE11_CPM <- (
    1e6 * subject_table$MRE11_count /
      pmax(subject_table$library_size, 1)
  )

  subject_table$eligible <- (
    subject_table$cell_count >= min_cells
  )

  write.table(
    subject_table,
    file.path(
      out_dir,
      paste0(
        analysis_name,
        "_MRE11_by_subject.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  eligible_table <- subject_table[
    subject_table$eligible,
    ,
    drop = FALSE
  ]

  eligible_table$group <- factor(
    eligible_table$group,
    levels = c("Control", "CWP")
  )

  eligible_table$series <- factor(
    eligible_table$series,
    levels = c("G227", "G288", "G135")
  )

  eligible_table$logit_detection <- qlogis(
    (
      eligible_table$detected_cells + 0.5
    ) /
      (
        eligible_table$cell_count + 1
      )
  )

  detection_lm <- lm(
    logit_detection ~ series + group,
    data = eligible_table
  )

  detection_coef <- summary(
    detection_lm
  )$coefficients

  detection_result <- data.frame(
    analysis = analysis_name,
    method = "DetectionLM",
    estimate = detection_coef[
      "groupCWP",
      "Estimate"
    ],
    std_error = detection_coef[
      "groupCWP",
      "Std. Error"
    ],
    statistic = detection_coef[
      "groupCWP",
      "t value"
    ],
    p_value = detection_coef[
      "groupCWP",
      "Pr(>|t|)"
    ],
    stringsAsFactors = FALSE
  )

  wilcox_result <- wilcox.test(
    detection_pct ~ group,
    data = eligible_table,
    exact = FALSE
  )

  group_summary <- data.frame(
    analysis = analysis_name,
    control_n = sum(
      eligible_table$group == "Control"
    ),
    cwp_n = sum(
      eligible_table$group == "CWP"
    ),
    control_detection_median = median(
      eligible_table$detection_pct[
        eligible_table$group == "Control"
      ]
    ),
    cwp_detection_median = median(
      eligible_table$detection_pct[
        eligible_table$group == "CWP"
      ]
    ),
    control_cpm_median = median(
      eligible_table$MRE11_CPM[
        eligible_table$group == "Control"
      ]
    ),
    cwp_cpm_median = median(
      eligible_table$MRE11_CPM[
        eligible_table$group == "CWP"
      ]
    ),
    wilcox_detection_p = wilcox_result$p.value,
    stringsAsFactors = FALSE
  )

  pb_eligible <- pb_counts[
    ,
    eligible_table$subject_id,
    drop = FALSE
  ]

  design <- model.matrix(
    ~ series + group,
    data = eligible_table
  )

  if (qr(design)$rank != ncol(design)) {
    stop(
      analysis_name,
      ": source-adjusted design is not full rank"
    )
  }

  rownames(design) <- eligible_table$subject_id

  y <- DGEList(
    counts = as.matrix(pb_eligible)
  )

  keep_gene <- filterByExpr(
    y,
    design = design
  )

  keep_gene[
    rownames(y) == "MRE11"
  ] <- TRUE

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

  edge_result <- topTags(
    test,
    n = Inf,
    sort.by = "none"
  )$table

  edge_result$gene <- rownames(edge_result)

  mre11_edge <- edge_result[
    edge_result$gene == "MRE11",
    ,
    drop = FALSE
  ]

  mre11_edge$analysis <- analysis_name

  write.table(
    mre11_edge,
    file.path(
      out_dir,
      paste0(
        analysis_name,
        "_MRE11_forced_edgeR.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p_det <- ggplot(
    eligible_table,
    aes(
      x = group,
      y = detection_pct,
      color = group,
      shape = series,
      label = subject_id
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.55
    ) +
    geom_jitter(
      width = 0.10,
      size = 3
    ) +
    geom_text(
      vjust = -0.7,
      size = 3,
      check_overlap = TRUE
    ) +
    labs(
      title = paste("MRE11 det", analysis_name),
      x = "Group",
      y = "Det (%)",
      color = "Group",
      shape = "Series",
      caption = paste(
        "Det: percentage of CD4 cells with non-zero MRE11 counts;",
        "each point represents one subject."
      )
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(
      figure_dir,
      paste0(
        analysis_name,
        "_MRE11_detection.png"
      )
    ),
    p_det,
    width = 9,
    height = 7,
    dpi = 260
  )

  p_cpm <- ggplot(
    eligible_table,
    aes(
      x = group,
      y = log10(MRE11_CPM + 0.01),
      color = group,
      shape = series,
      label = subject_id
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.55
    ) +
    geom_jitter(
      width = 0.10,
      size = 3
    ) +
    geom_text(
      vjust = -0.7,
      size = 3,
      check_overlap = TRUE
    ) +
    labs(
      title = paste("MRE11 PB", analysis_name),
      x = "Group",
      y = "log10 CPM",
      color = "Group",
      shape = "Series",
      caption = paste(
        "PB: subject-level pseudobulk;",
        "CPM: counts per million."
      )
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(
      figure_dir,
      paste0(
        analysis_name,
        "_MRE11_PB_CPM.png"
      )
    ),
    p_cpm,
    width = 9,
    height = 7,
    dpi = 260
  )

  list(
    subject_table = subject_table,
    group_summary = group_summary,
    detection_result = detection_result,
    edge_result = mre11_edge
  )
}

primary_result <- run_mre11_analysis(
  rep(TRUE, ncol(counts)),
  "Primary"
)

core_result <- run_mre11_analysis(
  as.character(meta$cd4_label) !=
    "CD4state",
  "Core"
)

message("===== SUBTYPE DETECTION =====")
subtype_levels <- c(
  "CD4-GZMA",
  "CD4-IL7R",
  "CD4-CCR7",
  "CD4state",
  "Treg"
)

subtype_subject <- do.call(
  rbind,
  lapply(subtype_levels, function(subtype_i) {
    selector <- (
      as.character(meta$cd4_label) ==
        subtype_i
    )

    z_meta <- meta[
      selector,
      ,
      drop = FALSE
    ]

    z_detected <- mre11_detected[
      selector
    ]

    subjects <- sort(
      unique(z_meta$subject_id2)
    )

    if (length(subjects) == 0) {
      return(NULL)
    }

    cell_n <- table(
      factor(
        z_meta$subject_id2,
        levels = subjects
      )
    )

    detected_n <- tapply(
      z_detected,
      factor(
        z_meta$subject_id2,
        levels = subjects
      ),
      sum
    )

    data.frame(
      subtype = subtype_i,
      subject_id = subjects,
      group = unname(
        tapply(
          as.character(z_meta$group2),
          factor(
            z_meta$subject_id2,
            levels = subjects
          ),
          function(x) unique(x)[1]
        )
      ),
      series = unname(
        tapply(
          as.character(z_meta$series),
          factor(
            z_meta$subject_id2,
            levels = subjects
          ),
          function(x) unique(x)[1]
        )
      ),
      cell_count = as.integer(cell_n),
      detected_cells = as.integer(
        detected_n
      ),
      stringsAsFactors = FALSE
    )
  })
)

subtype_subject$detection_pct <- (
  100 * subtype_subject$detected_cells /
    pmax(subtype_subject$cell_count, 1)
)

write.table(
  subtype_subject,
  file.path(
    out_dir,
    "MRE11_detection_by_subtype_subject.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

subtype_tests <- do.call(
  rbind,
  lapply(subtype_levels, function(subtype_i) {
    z <- subtype_subject[
      subtype_subject$subtype == subtype_i &
        subtype_subject$cell_count >= 10,
      ,
      drop = FALSE
    ]

    control_values <- z$detection_pct[
      z$group == "Control"
    ]

    cwp_values <- z$detection_pct[
      z$group == "CWP"
    ]

    if (
      length(control_values) < 3 ||
        length(cwp_values) < 3
    ) {
      return(
        data.frame(
          subtype = subtype_i,
          control_n = length(control_values),
          cwp_n = length(cwp_values),
          control_median = if (
            length(control_values) > 0
          ) median(control_values) else NA,
          cwp_median = if (
            length(cwp_values) > 0
          ) median(cwp_values) else NA,
          p_value = NA,
          stringsAsFactors = FALSE
        )
      )
    }

    wt <- wilcox.test(
      control_values,
      cwp_values,
      exact = FALSE
    )

    data.frame(
      subtype = subtype_i,
      control_n = length(control_values),
      cwp_n = length(cwp_values),
      control_median = median(
        control_values
      ),
      cwp_median = median(
        cwp_values
      ),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
)

subtype_tests$fdr_bh <- p.adjust(
  subtype_tests$p_value,
  method = "BH"
)

write.table(
  subtype_tests,
  file.path(
    out_dir,
    "MRE11_subtype_detection_tests.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

combined_summary <- rbind(
  primary_result$group_summary,
  core_result$group_summary
)

combined_detection <- rbind(
  primary_result$detection_result,
  core_result$detection_result
)

write.table(
  combined_summary,
  file.path(
    out_dir,
    "MRE11_group_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  combined_detection,
  file.path(
    out_dir,
    "MRE11_source_adjusted_detection_model.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("===== GROUP SUMMARY =====")
print(combined_summary, row.names = FALSE)

message("===== SOURCE-ADJUSTED DETECTION =====")
print(combined_detection, row.names = FALSE)

message("===== FORCED EDGER =====")
print(
  rbind(
    primary_result$edge_result,
    core_result$edge_result
  ),
  row.names = FALSE
)

message("===== SUBTYPE RESULTS =====")
print(subtype_tests, row.names = FALSE)

message("===== FINAL OUTPUT =====")
message("Output: ", out_dir)
message("Figures: ", figure_dir)
message("STATUS=PASS")
