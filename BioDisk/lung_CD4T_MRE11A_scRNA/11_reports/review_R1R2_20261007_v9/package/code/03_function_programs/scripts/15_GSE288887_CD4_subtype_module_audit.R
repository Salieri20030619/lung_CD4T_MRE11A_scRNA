suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260813)

root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

object_file <- file.path(
  root,
  "06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds"
)

score_file <- file.path(
  root,
  paste0(
    "08_MRE11A_analysis/GSE288887_CD4/module_audit_v1/",
    "GSE288887_CD4_module_cell_scores.tsv.gz"
  )
)

module_file <- file.path(
  root,
  paste0(
    "08_MRE11A_analysis/GSE288887_CD4/module_audit_v1/",
    "GSE288887_CD4_module_audit_results.rds"
  )
)

tab_dir <- file.path(
  root,
  "10_tables/CD4_final/subtype_module_audit_v1"
)

fig_dir <- file.path(
  root,
  "09_figures/CD4_final/subtype_module_audit_v1"
)

out_dir <- file.path(
  root,
  "08_MRE11A_analysis/GSE288887_CD4/subtype_module_audit_v1"
)

dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

stopifnot(
  file.exists(object_file),
  file.exists(score_file),
  file.exists(module_file)
)

cat("===== READ INPUTS =====\n")

obj <- readRDS(object_file)
obj <- tryCatch(
  JoinLayers(obj, assay="RNA"),
  error=function(e) obj
)

DefaultAssay(obj) <- "RNA"

counts <- GetAssayData(
  obj,
  assay="RNA",
  layer="counts"
)

cell_scores <- read.delim(
  score_file,
  check.names=FALSE,
  stringsAsFactors=FALSE
)

module_results <- readRDS(module_file)
signatures <- module_results$signatures

module_order <- c(
  "MitoResp", "MitoISR", "MitoROS", "Mitophagy",
  "PyroCore", "PyroSens", "PyroExec", "PyroNC", "GSDMAlt"
)

required_columns <- c(
  "cell_id", "Subject", "Group", "Series", "Subtype",
  module_order
)

missing_columns <- setdiff(
  required_columns,
  colnames(cell_scores)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse=", ")
  )
}

cell_scores <- cell_scores[
  cell_scores$cell_id %in% colnames(obj),
]

cell_scores <- cell_scores[
  match(colnames(obj), cell_scores$cell_id),
]

stopifnot(
  identical(cell_scores$cell_id, colnames(obj))
)

cell_scores$Group <- factor(
  cell_scores$Group,
  levels=c("Control", "CWP")
)

cell_scores$Series <- factor(
  cell_scores$Series,
  levels=c("G227", "G288", "G135")
)

subtype_order <- c(
  "CD4-CCR7",
  "CD4-GZMA",
  "CD4-IL7R",
  "CD4state",
  "Treg"
)

cell_scores$Subtype <- factor(
  cell_scores$Subtype,
  levels=subtype_order
)

cat("===== SUBJECT-SUBTYPE COUNTS =====\n")

subject_counts <- cell_scores %>%
  count(
    Subtype,
    Subject,
    Group,
    Series,
    name="n_cells"
  )

write.table(
  subject_counts,
  file.path(tab_dir, "subtype_subject_cell_counts.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

eligibility <- subject_counts %>%
  filter(n_cells >= 10) %>%
  group_by(Subtype, Group, Series) %>%
  summarise(
    subject_n=n_distinct(Subject),
    .groups="drop"
  ) %>%
  complete(
    Subtype,
    Group,
    Series,
    fill=list(subject_n=0)
  )

write.table(
  eligibility,
  file.path(tab_dir, "subtype_subject_eligibility.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

print(eligibility)

cat("===== SUBJECT-SUBTYPE MODULE SCORES =====\n")

subject_module <- cell_scores %>%
  group_by(
    Subtype,
    Subject,
    Group,
    Series
  ) %>%
  summarise(
    n_cells=n(),
    across(
      all_of(module_order),
      ~mean(.x, na.rm=TRUE)
    ),
    .groups="drop"
  ) %>%
  filter(n_cells >= 10)

write.table(
  subject_module,
  file.path(tab_dir, "subtype_subject_module_scores.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

test_one_module <- function(dat, module_name) {
  control_n <- sum(dat$Group == "Control")
  cwp_n <- sum(dat$Group == "CWP")

  g227_control_n <- sum(
    dat$Group == "Control" &
    dat$Series == "G227"
  )

  g227_cwp_n <- sum(
    dat$Group == "CWP" &
    dat$Series == "G227"
  )

  estimable <- (
    control_n >= 3 &&
    cwp_n >= 3 &&
    g227_control_n >= 3 &&
    g227_cwp_n >= 2
  )

  if (!estimable) {
    return(data.frame(
      estimate=NA_real_,
      std_error=NA_real_,
      statistic=NA_real_,
      p_value=NA_real_,
      control_n=control_n,
      cwp_n=cwp_n,
      g227_control_n=g227_control_n,
      g227_cwp_n=g227_cwp_n,
      status="INSUFFICIENT"
    ))
  }

  dat$Group <- factor(
    dat$Group,
    levels=c("Control", "CWP")
  )

  dat$Series <- droplevels(factor(dat$Series))
  model_rhs <- if (nlevels(dat$Series) > 1) "Series + Group" else "Group"

  fit <- lm(
    as.formula(
      paste(module_name, "~", model_rhs)
    ),
    data=dat
  )

  cc <- coef(summary(fit))

  if (!"GroupCWP" %in% rownames(cc)) {
    return(data.frame(
      estimate=NA_real_,
      std_error=NA_real_,
      statistic=NA_real_,
      p_value=NA_real_,
      control_n=control_n,
      cwp_n=cwp_n,
      g227_control_n=g227_control_n,
      g227_cwp_n=g227_cwp_n,
      status="NOT_ESTIMABLE"
    ))
  }

  data.frame(
    estimate=cc["GroupCWP", "Estimate"],
    std_error=cc["GroupCWP", "Std. Error"],
    statistic=cc["GroupCWP", "t value"],
    p_value=cc["GroupCWP", "Pr(>|t|)"],
    control_n=control_n,
    cwp_n=cwp_n,
    g227_control_n=g227_control_n,
    g227_cwp_n=g227_cwp_n,
    status="PASS"
  )
}

module_tests <- bind_rows(
  lapply(subtype_order, function(st) {
    dat <- subject_module[
      subject_module$Subtype == st,
    ]

    bind_rows(
      lapply(module_order, function(mod) {
        result <- test_one_module(
          dat,
          mod
        )

        cbind(
          data.frame(
            subtype=st,
            module=mod,
            stringsAsFactors=FALSE
          ),
          result
        )
      })
    )
  })
)

module_tests <- module_tests %>%
  group_by(subtype) %>%
  mutate(
    fdr_bh=p.adjust(
      p_value,
      method="BH"
    )
  ) %>%
  ungroup()

write.table(
  module_tests,
  file.path(tab_dir, "subtype_module_tests.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("===== MODULE-GENE PSEUDOBULK =====\n")

module_genes <- unique(
  intersect(
    unlist(signatures),
    rownames(counts)
  )
)

run_subtype_edger <- function(subtype_name) {
  meta <- cell_scores[
    cell_scores$Subtype == subtype_name,
    c(
      "cell_id",
      "Subject",
      "Group",
      "Series"
    )
  ]

  cell_number <- table(meta$Subject)
  eligible_subjects <- names(
    cell_number[cell_number >= 10]
  )

  meta <- meta[
    meta$Subject %in% eligible_subjects,
  ]

  subject_meta <- meta %>%
    distinct(
      Subject,
      Group,
      Series
    ) %>%
    arrange(Subject)

  control_n <- sum(
    subject_meta$Group == "Control"
  )

  cwp_n <- sum(
    subject_meta$Group == "CWP"
  )

  g227_control_n <- sum(
    subject_meta$Group == "Control" &
    subject_meta$Series == "G227"
  )

  g227_cwp_n <- sum(
    subject_meta$Group == "CWP" &
    subject_meta$Series == "G227"
  )

  estimable <- (
    control_n >= 3 &&
    cwp_n >= 3 &&
    g227_control_n >= 3 &&
    g227_cwp_n >= 2
  )

  if (!estimable) {
    return(data.frame(
      subtype=subtype_name,
      gene=module_genes,
      logFC=NA_real_,
      logCPM=NA_real_,
      PValue=NA_real_,
      FDR=NA_real_,
      control_n=control_n,
      cwp_n=cwp_n,
      status="INSUFFICIENT"
    ))
  }

  pb_list <- lapply(
    subject_meta$Subject,
    function(subject_name) {
      selected_cells <- meta$cell_id[
        meta$Subject == subject_name
      ]

      Matrix::rowSums(
        counts[
          ,
          selected_cells,
          drop=FALSE
        ]
      )
    }
  )

  pb <- do.call(cbind, pb_list)
  rownames(pb) <- rownames(counts)
  colnames(pb) <- subject_meta$Subject

  subject_meta$Group <- factor(
    subject_meta$Group,
    levels=c("Control", "CWP")
  )

  subject_meta$Group <- droplevels(subject_meta$Group)
  subject_meta$Series <- droplevels(factor(subject_meta$Series))
  design_formula <- if (nlevels(subject_meta$Series) > 1) {
    ~ Series + Group
  } else {
    ~ Group
  }

  design <- model.matrix(
    design_formula,
    data=subject_meta
  )

  y <- DGEList(
    counts=pb,
    samples=subject_meta
  )

  keep <- filterByExpr(y, design=design)

  if (sum(keep) < 2) {
    return(data.frame(
      subtype=subtype_name,
      gene=module_genes,
      logFC=NA_real_,
      logCPM=NA_real_,
      PValue=NA_real_,
      FDR=NA_real_,
      control_n=control_n,
      cwp_n=cwp_n,
      status="LOW_EXPRESSION"
    ))
  }

  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust=TRUE)

  fit <- glmQLFit(
    y,
    design,
    robust=TRUE
  )

  qlf <- glmQLFTest(
    fit,
    coef="GroupCWP"
  )

  result <- topTags(
    qlf,
    n=Inf,
    sort.by="none"
  )$table

  result$gene <- rownames(result)
  result <- result[result$gene %in% module_genes, , drop=FALSE]

  result <- result[, c(
    "gene",
    "logFC",
    "logCPM",
    "PValue",
    "FDR"
  )]

  result$subtype <- subtype_name
  result$control_n <- control_n
  result$cwp_n <- cwp_n
  result$status <- "PASS"

  result[, c(
    "subtype",
    "gene",
    "logFC",
    "logCPM",
    "PValue",
    "FDR",
    "control_n",
    "cwp_n",
    "status"
  )]
}

gene_tests <- bind_rows(
  lapply(
    subtype_order,
    run_subtype_edger
  )
)

write.table(
  gene_tests,
  file.path(tab_dir, "subtype_module_gene_edger.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("===== FIGURES =====\n")

plot_dat <- subject_module %>%
  select(
    Subtype,
    Subject,
    Group,
    Series,
    MitoResp,
    PyroCore
  ) %>%
  pivot_longer(
    cols=c(MitoResp, PyroCore),
    names_to="Module",
    values_to="Score"
  )

plot_dat$Subtype <- factor(
  plot_dat$Subtype,
  levels=subtype_order
)

plot_dat$Module <- factor(
  plot_dat$Module,
  levels=c("MitoResp", "PyroCore")
)

p_box <- ggplot(
  plot_dat,
  aes(
    x=Group,
    y=Score,
    colour=Group
  )
) +
  geom_boxplot(
    width=0.60,
    outlier.shape=NA
  ) +
  geom_point(
    aes(shape=Series),
    position=position_jitter(
      width=0.10,
      height=0
    ),
    size=2.3
  ) +
  facet_grid(
    Module ~ Subtype,
    scales="free_y"
  ) +
  scale_colour_manual(
    values=c(
      Control="#F8766D",
      CWP="#00BFC4"
    )
  ) +
  labs(
    title="Subtype module audit",
    x=NULL,
    y="UCell",
    caption=paste0(
      "MitoResp: mitochondrial response; ",
      "PyroCore: core pyroptosis-related program. ",
      "Each point represents one subject with at least 10 cells."
    )
  ) +
  theme_bw(base_size=11) +
  theme(
    legend.position="right",
    strip.background=element_rect(
      fill="white"
    ),
    axis.text.x=element_text(
      angle=35,
      hjust=1
    )
  )

ggsave(
  file.path(
    fig_dir,
    "GSE288887_CD4_subtype_module_audit.png"
  ),
  p_box,
  width=18,
  height=7,
  dpi=300
)

heat_dat <- module_tests %>%
  filter(status == "PASS") %>%
  mutate(
    subtype=factor(
      subtype,
      levels=subtype_order
    ),
    module=factor(
      module,
      levels=module_order
    ),
    label=case_when(
      !is.na(fdr_bh) & fdr_bh < 0.05 ~ "*",
      !is.na(p_value) & p_value < 0.05 ~ "+",
      TRUE ~ ""
    )
  )

if (nrow(heat_dat) > 0) {
  p_heat <- ggplot(
    heat_dat,
    aes(
      x=module,
      y=subtype,
      fill=estimate
    )
  ) +
    geom_tile(
      colour="white"
    ) +
    geom_text(
      aes(label=label),
      size=5
    ) +
    scale_fill_gradient2(
      low="#2166AC",
      mid="white",
      high="#B2182B",
      midpoint=0
    ) +
    labs(
      title="Adjusted module effects",
      x=NULL,
      y=NULL,
      fill="Effect",
      caption=paste0(
        "Positive values indicate higher scores in CWP. ",
        "*: FDR<0.05; +: nominal P<0.05."
      )
    ) +
    theme_bw(base_size=11) +
    theme(
      panel.grid=element_blank(),
      axis.text.x=element_text(
        angle=40,
        hjust=1
      )
    )

  ggsave(
    file.path(
      fig_dir,
      "GSE288887_CD4_subtype_module_effects.png"
    ),
    p_heat,
    width=11,
    height=5,
    dpi=300
  )
}

saveRDS(
  list(
    subject_counts=subject_counts,
    eligibility=eligibility,
    subject_module=subject_module,
    module_tests=module_tests,
    gene_tests=gene_tests
  ),
  file.path(
    out_dir,
    "GSE288887_CD4_subtype_module_audit.rds"
  )
)

cat("===== SUBTYPE MODULE SUMMARY =====\n")

print(
  module_tests %>%
    filter(
      module %in% c(
        "MitoResp",
        "PyroCore"
      )
    )
)

cat("===== NOMINAL MODULE RESULTS =====\n")

print(
  module_tests %>%
    filter(
      status == "PASS",
      !is.na(p_value),
      p_value < 0.10
    ) %>%
    arrange(p_value)
)

cat("===== MODULE GENE RESULTS =====\n")

print(
  gene_tests %>%
    filter(
      status == "PASS",
      !is.na(PValue),
      PValue < 0.05
    ) %>%
    arrange(FDR, PValue) %>%
    select(
      subtype,
      gene,
      logFC,
      PValue,
      FDR
    )
)

cat("===== FINAL OUTPUT =====\n")
cat("Tables:", tab_dir, "\n")
cat("Figures:", fig_dir, "\n")
cat("Results:", out_dir, "\n")
cat("STATUS=PASS\n")
