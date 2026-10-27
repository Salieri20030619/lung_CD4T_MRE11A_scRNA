#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)
set.seed(20260814)

project_root <- Sys.getenv(
  "PROJECT_ROOT",
  unset = "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
)

rank_dir <- file.path(
  project_root,
  "07_DEG", "CD4_pseudobulk_v1", "source_adjusted"
)
fig_dir <- file.path(
  project_root,
  "09_figures", "CD4_final", "pathway_v1", "multicurve_v2"
)
tab_dir <- file.path(
  project_root,
  "10_tables", "CD4_final", "pathway_v1", "multicurve_v2"
)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

rank_files <- c(
  Primary = file.path(rank_dir, "Primary_SrcAdj_edgeR.tsv"),
  Core = file.path(rank_dir, "Core_SrcAdj_edgeR.tsv")
)

stopifnot(all(file.exists(rank_files)))

pick_col <- function(nms, candidates, label) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0L) {
    stop(label, " column not found. Available: ", paste(nms, collapse = ", "))
  }
  hit[[1]]
}

read_rank <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  gene_col <- pick_col(names(x), c("gene", "Gene", "symbol", "SYMBOL"), "gene")
  lfc_col <- pick_col(names(x), c("logFC", "log2FC", "estimate"), "logFC")
  stat_col <- intersect(c("F", "LR", "statistic", "t", "Wald"), names(x))

  if (length(stat_col) > 0L) {
    s <- suppressWarnings(as.numeric(x[[stat_col[[1]]]]))
    lfc <- suppressWarnings(as.numeric(x[[lfc_col]]))
    metric <- sign(lfc) * sqrt(pmax(s, 0))
    rank_method <- paste0("signed_sqrt_", stat_col[[1]])
  } else {
    metric <- suppressWarnings(as.numeric(x[[lfc_col]]))
    rank_method <- lfc_col
  }

  out <- tibble(
    gene = toupper(as.character(x[[gene_col]])),
    metric = metric
  ) %>%
    filter(!is.na(gene), nzchar(gene), is.finite(metric)) %>%
    arrange(desc(abs(metric))) %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(desc(metric))

  ranks <- out$metric
  names(ranks) <- out$gene
  attr(ranks, "rank_method") <- rank_method
  ranks
}

message("===== LOAD LOCAL MSIGDB =====")
msig <- msigdbr(db_species = "HS", species = "Homo sapiens")

collection_col <- pick_col(
  names(msig), c("gs_collection", "gs_cat"), "MSigDB collection"
)
subcollection_col <- pick_col(
  names(msig), c("gs_subcollection", "gs_subcat"), "MSigDB subcollection"
)
term_col <- pick_col(names(msig), c("gs_name", "gs_id"), "MSigDB term")
symbol_col <- pick_col(
  names(msig), c("gene_symbol", "human_gene_symbol"), "MSigDB gene symbol"
)

make_term2gene <- function(collection, subcollection) {
  msig %>%
    filter(
      .data[[collection_col]] == collection,
      .data[[subcollection_col]] == subcollection
    ) %>%
    transmute(
      term = as.character(.data[[term_col]]),
      gene = toupper(as.character(.data[[symbol_col]]))
    ) %>%
    filter(nzchar(term), nzchar(gene)) %>%
    distinct()
}

go_t2g <- make_term2gene("C5", "GO:BP")
kegg_t2g <- make_term2gene("C2", "CP:KEGG_LEGACY")
message("GO pathways: ", n_distinct(go_t2g$term))
message("KEGG pathways: ", n_distinct(kegg_t2g$term))

run_gsea <- function(ranks, term2gene) {
  clusterProfiler::GSEA(
    geneList = sort(ranks, decreasing = TRUE),
    TERM2GENE = term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  )
}

running_es <- function(ranks, genes, exponent = 1) {
  n <- length(ranks)
  hit <- names(ranks) %in% genes
  nh <- sum(hit)
  if (nh == 0L || nh == n) return(NULL)

  weights <- abs(ranks)^exponent
  hit_denom <- sum(weights[hit])
  if (!is.finite(hit_denom) || hit_denom <= 0) return(NULL)

  increment <- rep(-1 / (n - nh), n)
  increment[hit] <- weights[hit] / hit_denom
  es <- cumsum(increment)
  tibble(
    rank = seq_len(n),
    es = es,
    hit = hit
  )
}

short_code <- function(database, n) {
  prefix <- ifelse(database == "GO_BP", "G", "K")
  paste0(prefix, seq_len(n))
}

select_top <- function(gsea_obj, database, scenario, top_n = 5L) {
  x <- as.data.frame(gsea_obj@result) %>%
    filter(is.finite(NES), is.finite(p.adjust))

  if (scenario == "NonRP") {
    x <- x %>%
      filter(!str_detect(
        toupper(paste(ID, Description)),
        "RIBOSOM|CYTOPLASMIC_TRANSLATION|TRANSLATION_INITIATION|TRANSLATION_ELONGATION"
      ))
  }

  x %>%
    arrange(desc(abs(NES)), p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(code = short_code(database, n()))
}

make_multicurve <- function(
  ranks,
  gsea_obj,
  term2gene,
  analysis,
  database,
  scenario,
  top_n = 5L
) {
  selected <- select_top(gsea_obj, database, scenario, top_n)
  if (nrow(selected) == 0L) {
    stop("No pathways selected: ", analysis, "/", database, "/", scenario)
  }

  palette <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00")
  names(palette) <- selected$code

  curve_list <- list()
  hit_list <- list()
  for (i in seq_len(nrow(selected))) {
    term <- selected$ID[[i]]
    code <- selected$code[[i]]
    genes <- term2gene$gene[term2gene$term == term]
    curve <- running_es(ranks, genes)
    if (is.null(curve)) next
    curve$code <- code
    curve_list[[code]] <- curve
    hit_list[[code]] <- curve %>% filter(hit) %>% transmute(rank, code = code)
  }

  curves <- bind_rows(curve_list)
  hits <- bind_rows(hit_list)
  selected <- selected %>% filter(code %in% unique(curves$code))
  selected$code <- factor(selected$code, levels = selected$code)
  curves$code <- factor(curves$code, levels = levels(selected$code))
  hits$code <- factor(hits$code, levels = rev(levels(selected$code)))

  legend_labels <- setNames(
    sprintf(
      "%s  NES=%+.2f; FDR=%.2g",
      selected$code, selected$NES, selected$p.adjust
    ),
    as.character(selected$code)
  )

  p_es <- ggplot(curves, aes(rank, es, color = code)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_line(linewidth = 0.85, alpha = 0.95) +
    scale_color_manual(values = palette, labels = legend_labels, drop = FALSE) +
    labs(
      title = paste(database, "GSEA"),
      subtitle = paste(analysis, scenario, "| positive ranks indicate CWP higher"),
      x = NULL,
      y = "Running ES",
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.key.width = grid::unit(12, "pt")
    )

  hit_levels <- levels(hits$code)
  hits$code_y <- as.numeric(hits$code)
  p_hit <- ggplot(hits) +
    geom_segment(
      aes(
        x = rank,
        xend = rank,
        y = code_y - 0.34,
        yend = code_y + 0.34,
        color = code
      ),
      linewidth = 0.32,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_y_continuous(
      breaks = seq_along(hit_levels),
      labels = hit_levels,
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    labs(x = NULL, y = "Hits") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.y = element_blank()
    )

  rank_df <- tibble(
    rank = seq_along(ranks),
    metric = as.numeric(ranks),
    direction = ifelse(ranks >= 0, "CWP", "Control")
  )
  p_rank <- ggplot(rank_df, aes(rank, metric, fill = direction)) +
    geom_area(linewidth = 0, alpha = 0.85) +
    scale_fill_manual(values = c(CWP = "#D55E00", Control = "#0072B2")) +
    labs(x = "Rank in ordered dataset", y = "Rank metric") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")

  combined <- p_es / p_hit / p_rank +
    plot_layout(heights = c(3.1, 0.9, 1.25), guides = "collect")

  stem <- paste(analysis, database, "GSEA", scenario, sep = "_")
  ggsave(
    file.path(fig_dir, paste0(stem, ".png")),
    combined,
    width = 11,
    height = 7.5,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    file.path(fig_dir, paste0(stem, ".pdf")),
    combined,
    width = 11,
    height = 7.5,
    device = grDevices::cairo_pdf
  )

  mapping <- selected %>%
    transmute(
      analysis = analysis,
      database = database,
      scenario = scenario,
      code = as.character(code),
      ID,
      Description,
      NES,
      pvalue,
      p.adjust,
      core_enrichment
    )
  write_tsv(mapping, file.path(tab_dir, paste0(stem, "_mapping.tsv")))
  mapping
}

all_maps <- list()
for (analysis in names(rank_files)) {
  message("===== ", analysis, " =====")
  ranks <- read_rank(rank_files[[analysis]])
  message("Ranked genes: ", length(ranks), "; method: ", attr(ranks, "rank_method"))

  go_obj <- run_gsea(ranks, go_t2g)
  kegg_obj <- run_gsea(ranks, kegg_t2g)

  saveRDS(go_obj, file.path(tab_dir, paste0(analysis, "_GO_GSEA.rds")))
  saveRDS(kegg_obj, file.path(tab_dir, paste0(analysis, "_KEGG_GSEA.rds")))
  write_tsv(as.data.frame(go_obj@result), file.path(tab_dir, paste0(analysis, "_GO_GSEA.tsv")))
  write_tsv(as.data.frame(kegg_obj@result), file.path(tab_dir, paste0(analysis, "_KEGG_GSEA.tsv")))

  for (scenario in c("Full", "NonRP")) {
    all_maps[[paste(analysis, "GO_BP", scenario, sep = "_")]] <-
      make_multicurve(ranks, go_obj, go_t2g, analysis, "GO_BP", scenario)
    all_maps[[paste(analysis, "KEGG", scenario, sep = "_")]] <-
      make_multicurve(ranks, kegg_obj, kegg_t2g, analysis, "KEGG", scenario)
  }
}

write_tsv(bind_rows(all_maps), file.path(tab_dir, "multicurve_pathway_mapping_all.tsv"))

message("===== OUTPUT =====")
message("Figures: ", fig_dir)
message("Tables: ", tab_dir)
message("STATUS=PASS")
