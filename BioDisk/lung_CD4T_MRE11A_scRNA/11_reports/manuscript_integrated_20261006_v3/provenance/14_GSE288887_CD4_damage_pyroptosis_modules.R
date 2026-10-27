suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260812)

root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

input_rds <- file.path(
  root, "06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds"
)

out_dir <- file.path(
  root, "08_MRE11A_analysis/GSE288887_CD4/module_audit_v1"
)

fig_dir <- file.path(
  root, "09_figures/MRE11_GSE288887/module_audit_v1"
)

tab_dir <- file.path(
  root, "10_tables/CD4_final/module_audit_v1"
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)

cat("===== READ FINAL CD4 OBJECT =====\n")
stopifnot(file.exists(input_rds))

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

obj <- tryCatch(
  JoinLayers(obj, assay="RNA"),
  error=function(e) obj
)

counts <- GetAssayData(obj, assay="RNA", layer="counts")

if (!"MRE11" %in% rownames(counts)) {
  stop("MRE11 not found in RNA counts")
}

md0 <- obj[[]]

pick_col <- function(candidates, label) {
  found <- candidates[candidates %in% colnames(md0)]
  if (length(found) == 0) {
    stop(
      paste0(
        "Cannot identify ", label, " column. Available: ",
        paste(colnames(md0), collapse=", ")
      )
    )
  }
  found[1]
}

subject_col <- pick_col(
  c("subject_id", "subject", "Subject", "donor_id", "sample_id"),
  "subject"
)

group_col <- pick_col(
  c("group", "Group", "condition", "Condition"),
  "group"
)

subtype_col <- pick_col(
  c("cd4_label", "CD4_anno", "cd4_anno", "Subtype", "subtype"),
  "CD4 subtype"
)

ncount_col <- pick_col(
  c("nCount_RNA", "nCount"),
  "nCount"
)

nfeature_col <- pick_col(
  c("nFeature_RNA", "nFeature"),
  "nFeature"
)

pctmt_col <- pick_col(
  c("percent.mt", "percent_mito", "pct_mito"),
  "mitochondrial percentage"
)

normalize_subject <- function(x) {
  x <- toupper(as.character(x))
  m <- regexpr("(CTL|CWP)[0-9]+", x)
  ans <- rep(NA_character_, length(x))
  ok <- m > 0
  ans[ok] <- regmatches(x, m)[ok]
  ans
}

md0$Subject <- normalize_subject(md0[[subject_col]])

if (anyNA(md0$Subject)) {
  bad <- unique(as.character(md0[[subject_col]])[is.na(md0$Subject)])
  stop("Cannot normalize subjects: ", paste(bad, collapse=", "))
}

group_raw <- tolower(as.character(md0[[group_col]]))

md0$Group <- ifelse(
  grepl("cwp", group_raw), "CWP",
  ifelse(
    grepl("control|ctl", group_raw), "Control",
    as.character(md0[[group_col]])
  )
)

if (!all(md0$Group %in% c("Control", "CWP"))) {
  stop("Unexpected group values: ", paste(unique(md0$Group), collapse=", "))
}

g227 <- c(
  "CWP85", "CWP86",
  "CTL092", "CTL094", "CTL098", "CTL103", "CTL105",
  "CTL106", "CTL108", "CTL111", "CTL117"
)

g135 <- c("CTL068", "CTL070")

g288 <- c(
  "CWP89", "CWP94", "CWP109", "CWP127",
  "CWP144", "CWP149", "CWP152"
)

md0$Series <- ifelse(
  md0$Subject %in% g227, "G227",
  ifelse(
    md0$Subject %in% g135, "G135",
    ifelse(md0$Subject %in% g288, "G288", NA_character_)
  )
)

if (anyNA(md0$Series)) {
  stop(
    "Unmapped subjects: ",
    paste(unique(md0$Subject[is.na(md0$Series)]), collapse=", ")
  )
}

md0$Subtype <- as.character(md0[[subtype_col]])
rownames(md0) <- colnames(obj)

obj$Subject <- md0$Subject
obj$Group <- md0$Group
obj$Series <- md0$Series
obj$Subtype <- md0$Subtype

cat("===== DEFINE MODULES =====\n")

signatures <- list(
  MitoISR = c(
    "ATF4", "ATF5", "DDIT3", "DDIT4", "CHAC1", "ASNS",
    "TRIB3", "PPP1R15A", "SESN2", "DELE1", "EIF2AK1",
    "HSPA9", "HSPD1", "HSPE1", "LONP1", "CLPP", "OMA1"
  ),

  MitoROS = c(
    "SOD2", "PRDX3", "PRDX5", "GPX4", "TXN2", "TXNRD2",
    "HMOX1", "NQO1", "GCLC", "GCLM", "SRXN1", "FTH1", "FTL"
  ),

  Mitophagy = c(
    "PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1",
    "SQSTM1", "OPTN", "TBK1", "ULK1", "BECN1", "MAP1LC3B"
  ),

  PyroSens = c(
    "NLRP1", "NLRP3", "AIM2", "NLRC4",
    "NAIP", "MEFV", "PYCARD"
  ),

  PyroExec = c(
    "PYCARD", "CASP1", "GSDMD", "IL1B", "IL18"
  ),

  PyroNC = c(
    "CASP4", "CASP5", "GSDMD"
  ),

  GSDMAlt = c(
    "GSDME", "CASP3", "GSDMB", "GZMA", "GSDMC"
  )
)

signatures$MitoResp <- unique(c(
  signatures$MitoISR,
  signatures$MitoROS,
  signatures$Mitophagy
))

signatures$PyroCore <- unique(c(
  signatures$PyroSens,
  signatures$PyroExec,
  signatures$PyroNC
))

module_order <- c(
  "MitoResp", "MitoISR", "MitoROS", "Mitophagy",
  "PyroCore", "PyroSens", "PyroExec", "PyroNC", "GSDMAlt"
)

signatures <- signatures[module_order]

coverage <- bind_rows(lapply(names(signatures), function(nm) {
  expected <- unique(signatures[[nm]])
  present <- expected[expected %in% rownames(counts)]
  missing <- setdiff(expected, present)

  data.frame(
    module=nm,
    expected_n=length(expected),
    present_n=length(present),
    coverage_pct=100 * length(present) / length(expected),
    present_genes=paste(present, collapse=","),
    missing_genes=paste(missing, collapse=","),
    stringsAsFactors=FALSE
  )
}))

write.table(
  coverage,
  file.path(tab_dir, "module_gene_coverage.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

print(coverage[, c(
  "module", "expected_n", "present_n", "coverage_pct"
)])

signatures_use <- lapply(
  signatures,
  function(x) intersect(unique(x), rownames(counts))
)

too_small <- names(signatures_use)[lengths(signatures_use) < 3]

if (length(too_small) > 0) {
  stop(
    "Fewer than three detected genes in module(s): ",
    paste(too_small, collapse=", ")
  )
}

cat("===== UCELL SCORING =====\n")

obj <- UCell::AddModuleScore_UCell(
  obj,
  features=signatures_use,
  name=NULL
)

missing_scores <- setdiff(module_order, colnames(obj[[]]))

if (length(missing_scores) > 0) {
  stop(
    "UCell score columns missing: ",
    paste(missing_scores, collapse=", ")
  )
}

md <- obj[[]]

mre11_count <- as.numeric(counts["MRE11", colnames(obj)])

md$MRE11_count <- mre11_count
md$MRE11_detect <- as.numeric(mre11_count > 0)
md$cell_id <- rownames(md)

cell_output <- md[, c(
  "cell_id", "Subject", "Group", "Series", "Subtype",
  module_order, "MRE11_count", "MRE11_detect",
  ncount_col, nfeature_col, pctmt_col
)]

cell_con <- gzfile(
  file.path(out_dir, "GSE288887_CD4_module_cell_scores.tsv.gz"),
  open="wt"
)

write.table(
  cell_output, cell_con,
  sep="\t", quote=FALSE, row.names=FALSE
)

close(cell_con)

cat("===== SUBJECT AGGREGATION =====\n")

make_subject_table <- function(dat, analysis_name) {
  dat %>%
    group_by(Subject, Group, Series) %>%
    summarise(
      n_cells=n(),
      across(all_of(module_order), ~mean(.x, na.rm=TRUE)),
      MRE11_detection=100 * mean(MRE11_detect > 0),
      MRE11_mean_count=mean(MRE11_count),
      median_nCount=median(.data[[ncount_col]], na.rm=TRUE),
      median_nFeature=median(.data[[nfeature_col]], na.rm=TRUE),
      median_pct_mt=median(.data[[pctmt_col]], na.rm=TRUE),
      .groups="drop"
    ) %>%
    mutate(analysis=analysis_name)
}

subject_primary <- make_subject_table(md, "Primary")

subject_core <- make_subject_table(
  md[md$Subtype != "CD4state", , drop=FALSE],
  "Core"
)

subject_scores <- bind_rows(subject_primary, subject_core)

write.table(
  subject_scores,
  file.path(tab_dir, "module_subject_scores.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("===== SOURCE-ADJUSTED MODULE TESTS =====\n")

extract_group <- function(dat, response, model_name) {
  dat <- dat[dat$n_cells >= 20, , drop=FALSE]

  dat$Group <- factor(dat$Group, levels=c("Control", "CWP"))
  dat$Series <- factor(dat$Series, levels=c("G227", "G288", "G135"))
  dat$logN <- log10(dat$median_nCount + 1)

  formula_use <- if (model_name == "SourceAdj") {
    as.formula(paste(response, "~ Series + Group"))
  } else {
    as.formula(
      paste(
        response,
        "~ Series + Group + scale(logN) + scale(median_pct_mt)"
      )
    )
  }

  fit <- tryCatch(
    lm(formula_use, data=dat),
    error=function(e) NULL
  )

  if (is.null(fit)) {
    return(data.frame(
      estimate=NA_real_, std_error=NA_real_,
      statistic=NA_real_, p_value=NA_real_,
      control_n=sum(dat$Group == "Control"),
      cwp_n=sum(dat$Group == "CWP")
    ))
  }

  cc <- coef(summary(fit))

  if (!"GroupCWP" %in% rownames(cc)) {
    return(data.frame(
      estimate=NA_real_, std_error=NA_real_,
      statistic=NA_real_, p_value=NA_real_,
      control_n=sum(dat$Group == "Control"),
      cwp_n=sum(dat$Group == "CWP")
    ))
  }

  data.frame(
    estimate=cc["GroupCWP", "Estimate"],
    std_error=cc["GroupCWP", "Std. Error"],
    statistic=cc["GroupCWP", "t value"],
    p_value=cc["GroupCWP", "Pr(>|t|)"],
    control_n=sum(dat$Group == "Control"),
    cwp_n=sum(dat$Group == "CWP")
  )
}

group_tests <- bind_rows(lapply(c("Primary", "Core"), function(an) {
  dat <- subject_scores[subject_scores$analysis == an, ]

  bind_rows(lapply(c("SourceAdj", "QCAdj"), function(model_name) {
    bind_rows(lapply(module_order, function(mod) {
      z <- extract_group(dat, mod, model_name)

      cbind(
        data.frame(
          analysis=an,
          model=model_name,
          module=mod,
          stringsAsFactors=FALSE
        ),
        z
      )
    }))
  }))
}))

group_tests <- group_tests %>%
  group_by(analysis, model) %>%
  mutate(fdr_bh=p.adjust(p_value, method="BH")) %>%
  ungroup()

write.table(
  group_tests,
  file.path(tab_dir, "module_group_tests.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("===== WITHIN-G227 DIRECTION =====\n")

g227_results <- bind_rows(lapply(c("Primary", "Core"), function(an) {
  dat <- subject_scores %>%
    filter(
      analysis == an,
      Series == "G227",
      n_cells >= 20
    )

  bind_rows(lapply(module_order, function(mod) {
    ctl <- dat[[mod]][dat$Group == "Control"]
    cwp <- dat[[mod]][dat$Group == "CWP"]

    p <- if (length(ctl) > 0 && length(cwp) > 0) {
      wilcox.test(cwp, ctl, exact=FALSE)$p.value
    } else {
      NA_real_
    }

    data.frame(
      analysis=an,
      module=mod,
      control_n=length(ctl),
      cwp_n=length(cwp),
      control_mean=mean(ctl),
      cwp_mean=mean(cwp),
      delta_mean=mean(cwp) - mean(ctl),
      wilcox_p=p
    )
  }))
}))

g227_results <- g227_results %>%
  group_by(analysis) %>%
  mutate(fdr_bh=p.adjust(wilcox_p, method="BH")) %>%
  ungroup()

write.table(
  g227_results,
  file.path(tab_dir, "module_within_G227.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("===== MRE11 ASSOCIATIONS =====\n")

mre11_assoc <- bind_rows(lapply(c("Primary", "Core"), function(an) {
  dat <- subject_scores %>%
    filter(analysis == an, n_cells >= 20)

  bind_rows(lapply(c("MitoResp", "PyroCore"), function(mod) {
    ct <- suppressWarnings(
      cor.test(
        dat[[mod]],
        dat$MRE11_detection,
        method="spearman",
        exact=FALSE
      )
    )

    dat$Group <- factor(dat$Group, levels=c("Control", "CWP"))
    dat$Series <- factor(dat$Series, levels=c("G227", "G288", "G135"))
    dat$logN <- log10(dat$median_nCount + 1)

    fit <- lm(
      as.formula(
        paste(
          mod,
          "~ Series + Group + MRE11_detection + scale(logN)"
        )
      ),
      data=dat
    )

    cc <- coef(summary(fit))

    if ("MRE11_detection" %in% rownames(cc)) {
      adj_est <- cc["MRE11_detection", "Estimate"]
      adj_p <- cc["MRE11_detection", "Pr(>|t|)"]
    } else {
      adj_est <- NA_real_
      adj_p <- NA_real_
    }

    data.frame(
      analysis=an,
      module=mod,
      subject_n=nrow(dat),
      spearman_rho=unname(ct$estimate),
      spearman_p=ct$p.value,
      adjusted_estimate=adj_est,
      adjusted_p=adj_p
    )
  }))
}))

mre11_assoc$fdr_bh <- p.adjust(
  mre11_assoc$adjusted_p,
  method="BH"
)

write.table(
  mre11_assoc,
  file.path(tab_dir, "module_MRE11_association.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("===== MODULE GENE DETECTION =====\n")

all_module_genes <- unique(unlist(signatures_use))

gene_detection <- bind_rows(lapply(c("Control", "CWP"), function(g) {
  cells <- rownames(md)[md$Group == g]

  data.frame(
    group=g,
    gene=all_module_genes,
    detection_pct=100 * Matrix::rowMeans(
      counts[all_module_genes, cells, drop=FALSE] > 0
    ),
    mean_count=Matrix::rowMeans(
      counts[all_module_genes, cells, drop=FALSE]
    )
  )
}))

write.table(
  gene_detection,
  file.path(tab_dir, "module_gene_detection.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("===== FIGURES =====\n")

plot_dat <- subject_scores %>%
  filter(n_cells >= 20) %>%
  pivot_longer(
    cols=all_of(module_order),
    names_to="Module",
    values_to="Score"
  )

plot_dat$Module <- factor(
  plot_dat$Module,
  levels=module_order
)

plot_dat$Group <- factor(
  plot_dat$Group,
  levels=c("Control", "CWP")
)

plot_dat$Series <- factor(
  plot_dat$Series,
  levels=c("G227", "G288", "G135")
)

p_group <- ggplot(
  plot_dat,
  aes(x=Group, y=Score, colour=Group)
) +
  geom_boxplot(
    width=0.60,
    outlier.shape=NA
  ) +
  geom_point(
    aes(shape=Series),
    position=position_jitter(width=0.10, height=0),
    size=2.4,
    alpha=0.9
  ) +
  facet_grid(
    analysis ~ Module,
    scales="free_y"
  ) +
  scale_colour_manual(
    values=c(Control="#F8766D", CWP="#00BFC4")
  ) +
  labs(
    title="CD4 module scores",
    x=NULL,
    y="UCell"
  ) +
  theme_bw(base_size=12) +
  theme(
    legend.position="right",
    strip.background=element_rect(fill="white"),
    axis.text.x=element_text(angle=35, hjust=1)
  )

ggsave(
  file.path(fig_dir, "GSE288887_CD4_module_subject_scores.png"),
  p_group,
  width=20,
  height=8,
  dpi=300
)

if ("umap" %in% names(obj@reductions)) {
  p_umap <- FeaturePlot(
    obj,
    features=c("MitoResp", "PyroCore"),
    reduction="umap",
    order=TRUE,
    raster=TRUE,
    cols=c("grey92", "#5E3C99"),
    min.cutoff="q05",
    max.cutoff="q95",
    ncol=2
  ) &
    theme_bw(base_size=12)

  ggsave(
    file.path(fig_dir, "GSE288887_CD4_module_UMAP.png"),
    p_umap,
    width=13,
    height=6,
    dpi=300
  )
}

subtype_subject <- md %>%
  group_by(Subject, Group, Series, Subtype) %>%
  summarise(
    n_cells=n(),
    MitoResp=mean(MitoResp),
    PyroCore=mean(PyroCore),
    .groups="drop"
  ) %>%
  filter(n_cells >= 10) %>%
  pivot_longer(
    cols=c(MitoResp, PyroCore),
    names_to="Module",
    values_to="Score"
  )

subtype_subject$Group <- factor(
  subtype_subject$Group,
  levels=c("Control", "CWP")
)

p_subtype <- ggplot(
  subtype_subject,
  aes(x=Group, y=Score, colour=Group)
) +
  geom_boxplot(
    width=0.60,
    outlier.shape=NA
  ) +
  geom_point(
    aes(shape=Series),
    position=position_jitter(width=0.10),
    size=2
  ) +
  facet_grid(
    Module ~ Subtype,
    scales="free_y"
  ) +
  scale_colour_manual(
    values=c(Control="#F8766D", CWP="#00BFC4")
  ) +
  labs(
    title="Modules by CD4 subtype",
    x=NULL,
    y="UCell"
  ) +
  theme_bw(base_size=11) +
  theme(
    legend.position="right",
    strip.background=element_rect(fill="white"),
    axis.text.x=element_text(angle=35, hjust=1)
  )

ggsave(
  file.path(fig_dir, "GSE288887_CD4_module_by_subtype.png"),
  p_subtype,
  width=18,
  height=7,
  dpi=300
)

saveRDS(
  list(
    signatures=signatures,
    coverage=coverage,
    subject_scores=subject_scores,
    group_tests=group_tests,
    within_G227=g227_results,
    MRE11_association=mre11_assoc
  ),
  file.path(out_dir, "GSE288887_CD4_module_audit_results.rds")
)

cat("===== MAIN GROUP RESULTS =====\n")
print(
  group_tests %>%
    filter(
      model == "SourceAdj",
      module %in% c("MitoResp", "PyroCore")
    )
)

cat("===== MRE11 ASSOCIATIONS =====\n")
print(mre11_assoc)

cat("===== FINAL OUTPUT =====\n")
cat("Tables:", tab_dir, "\n")
cat("Figures:", fig_dir, "\n")
cat("Results:", out_dir, "\n")
cat("STATUS=PASS\n")

sessionInfo()