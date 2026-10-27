suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(GSVA)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

set.seed(20260813)

root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

input_rds <- file.path(
  root,
  "06_CD4T_analysis/CD4_final/GSE288887_CD4_final.rds"
)

out_dir <- file.path(
  root,
  "07_DEG/CD4_pseudobulk_v1/pathway_v1"
)

tab_dir <- file.path(
  root,
  "10_tables/CD4_final/pathway_v1"
)

fig_dir <- file.path(
  root,
  "09_figures/CD4_final/pathway_v1"
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(input_rds))

cat("===== READ FINAL CD4 OBJECT =====\n")

obj <- readRDS(input_rds)
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

md <- obj[[]]
md$cell_id <- rownames(md)

required_columns <- c(
  "subject_id",
  "group",
  "cd4_label"
)

missing_columns <- setdiff(
  required_columns,
  colnames(md)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse=", ")
  )
}

normalize_subject <- function(x) {
  x <- toupper(as.character(x))
  m <- regexpr("(CTL|CWP)[0-9]+", x)

  result <- rep(NA_character_, length(x))
  valid <- m > 0

  result[valid] <- regmatches(
    x,
    m
  )[valid]

  result
}

md$Subject <- normalize_subject(
  md$subject_id
)

group_raw <- tolower(
  as.character(md$group)
)

md$Group <- ifelse(
  grepl("cwp", group_raw),
  "CWP",
  ifelse(
    grepl("control|ctl", group_raw),
    "Control",
    NA_character_
  )
)

g227 <- c(
  "CWP85", "CWP86",
  "CTL092", "CTL094", "CTL098",
  "CTL103", "CTL105", "CTL106",
  "CTL108", "CTL111", "CTL117"
)

g135 <- c(
  "CTL068",
  "CTL070"
)

g288 <- c(
  "CWP89", "CWP94", "CWP109",
  "CWP127", "CWP144",
  "CWP149", "CWP152"
)

md$Series <- ifelse(
  md$Subject %in% g227,
  "G227",
  ifelse(
    md$Subject %in% g135,
    "G135",
    ifelse(
      md$Subject %in% g288,
      "G288",
      NA_character_
    )
  )
)

md$Subtype <- as.character(
  md$cd4_label
)

if (
  anyNA(md$Subject) ||
  anyNA(md$Group) ||
  anyNA(md$Series)
) {
  stop("Subject, group or series mapping contains NA")
}

cat("===== LOAD MSIGDB GENE SETS =====\n")

msig_all <- msigdbr(
  species="Homo sapiens"
)

collection_col <- if (
  "gs_collection" %in% colnames(msig_all)
) {
  "gs_collection"
} else {
  "gs_cat"
}

subcollection_col <- if (
  "gs_subcollection" %in% colnames(msig_all)
) {
  "gs_subcollection"
} else {
  "gs_subcat"
}

available_subcollections <- sort(
  unique(
    as.character(
      msig_all[[subcollection_col]]
    )
  )
)

go_subcollections <- available_subcollections[
  grepl(
    "GO[:_]?BP|GOBP",
    available_subcollections,
    ignore.case=TRUE
  )
]

kegg_subcollections <- available_subcollections[
  grepl(
    "KEGG",
    available_subcollections,
    ignore.case=TRUE
  )
]

if (length(go_subcollections) == 0) {
  stop("Cannot identify GO BP collection in msigdbr")
}

if (length(kegg_subcollections) == 0) {
  stop("Cannot identify KEGG collection in msigdbr")
}

go_subcollection <- go_subcollections[1]

legacy_kegg <- kegg_subcollections[
  grepl(
    "LEGACY",
    kegg_subcollections,
    ignore.case=TRUE
  )
]

kegg_subcollection <- if (
  length(legacy_kegg) > 0
) {
  legacy_kegg[1]
} else {
  kegg_subcollections[1]
}

cat(
  "GO_SUBCOLLECTION=",
  go_subcollection,
  "\n",
  sep=""
)

cat(
  "KEGG_SUBCOLLECTION=",
  kegg_subcollection,
  "\n",
  sep=""
)

description_col <- if (
  "gs_description" %in% colnames(msig_all)
) {
  "gs_description"
} else {
  "gs_name"
}

go_msig <- msig_all[
  as.character(msig_all[[subcollection_col]]) ==
    go_subcollection,
]

kegg_msig <- msig_all[
  as.character(msig_all[[subcollection_col]]) ==
    kegg_subcollection,
]

make_term_tables <- function(dat, prefix_pattern) {
  term2gene <- dat %>%
    transmute(
      term=as.character(gs_name),
      gene=as.character(gene_symbol)
    ) %>%
    filter(
      !is.na(term),
      !is.na(gene),
      term != "",
      gene != ""
    ) %>%
    distinct()

  term2name <- dat %>%
    transmute(
      term=as.character(gs_name),
      raw_name=as.character(.data[[description_col]])
    ) %>%
    filter(
      !is.na(term),
      !is.na(raw_name)
    ) %>%
    distinct(term, .keep_all=TRUE) %>%
    mutate(
      name=gsub(
        prefix_pattern,
        "",
        raw_name,
        ignore.case=TRUE
      ),
      name=gsub("_", " ", name),
      name=str_to_sentence(
        str_to_lower(name)
      )
    ) %>%
    select(term, name)

  list(
    term2gene=term2gene,
    term2name=term2name
  )
}

go_terms <- make_term_tables(
  go_msig,
  "^GOBP[ _:]?"
)

kegg_terms <- make_term_tables(
  kegg_msig,
  "^KEGG[ _:]?"
)

write.table(
  go_terms$term2name,
  file.path(tab_dir, "GO_term_names.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

write.table(
  kegg_terms$term2name,
  file.path(tab_dir, "KEGG_term_names.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat(
  "GO pathways:",
  n_distinct(go_terms$term2gene$term),
  "\n"
)

cat(
  "KEGG pathways:",
  n_distinct(kegg_terms$term2gene$term),
  "\n"
)

cat("===== PSEUDOBULK EDGE-R =====\n")

run_pseudobulk <- function(
  analysis_name,
  excluded_subtype=NULL
) {
  selected <- if (is.null(excluded_subtype)) {
    rep(TRUE, nrow(md))
  } else {
    md$Subtype != excluded_subtype
  }

  meta <- md[selected, ]

  subject_cells <- table(meta$Subject)

  eligible_subjects <- names(
    subject_cells[
      subject_cells >= 20
    ]
  )

  meta <- meta[
    meta$Subject %in% eligible_subjects,
  ]

  sample_meta <- meta %>%
    distinct(
      Subject,
      Group,
      Series
    ) %>%
    arrange(Subject)

  pb_list <- lapply(
    sample_meta$Subject,
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

  pb <- do.call(
    cbind,
    pb_list
  )

  rownames(pb) <- rownames(counts)
  colnames(pb) <- sample_meta$Subject

  sample_meta$Group <- factor(
    sample_meta$Group,
    levels=c("Control", "CWP")
  )

  sample_meta$Series <- factor(
    sample_meta$Series,
    levels=c("G227", "G288", "G135")
  )

  design <- model.matrix(
    ~ Series + Group,
    data=sample_meta
  )

  rownames(design) <- sample_meta$Subject

  y <- DGEList(
    counts=pb,
    samples=sample_meta
  )

  keep <- filterByExpr(
    y,
    design=design
  )

  y <- y[
    keep,
    ,
    keep.lib.sizes=FALSE
  ]

  y <- calcNormFactors(y)
  y <- estimateDisp(
    y,
    design,
    robust=TRUE
  )

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

  result$rank_score <- sign(
    result$logFC
  ) * sqrt(
    pmax(result$F, 0)
  )

  result$analysis <- analysis_name

  result <- result[, c(
    "analysis",
    "gene",
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR",
    "rank_score"
  )]

  log_cpm <- cpm(
    y,
    log=TRUE,
    prior.count=2
  )

  list(
    result=result,
    log_cpm=log_cpm,
    sample_meta=sample_meta,
    design=design
  )
}

primary <- run_pseudobulk(
  "Primary",
  excluded_subtype=NULL
)

core <- run_pseudobulk(
  "Core",
  excluded_subtype="CD4state"
)

deg_results <- bind_rows(
  primary$result,
  core$result
)

write.table(
  deg_results,
  file.path(tab_dir, "source_adjusted_full_DEG.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

deg_summary <- deg_results %>%
  group_by(analysis) %>%
  summarise(
    retained_genes=n(),
    fdr_lt_005=sum(FDR < 0.05),
    up_cwp=sum(
      FDR < 0.05 &
      logFC > 0
    ),
    down_cwp=sum(
      FDR < 0.05 &
      logFC < 0
    ),
    .groups="drop"
  )

write.table(
  deg_summary,
  file.path(tab_dir, "pathway_DEG_summary.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

print(deg_summary)

cat("===== GO AND KEGG ORA =====\n")

run_ora <- function(
  result,
  term_object,
  analysis_name,
  database_name
) {
  significant_genes <- result$gene[
    result$FDR < 0.05
  ]

  universe_genes <- result$gene

  if (length(significant_genes) < 3) {
    return(data.frame())
  }

  ora <- enricher(
    gene=significant_genes,
    universe=universe_genes,
    TERM2GENE=term_object$term2gene,
    TERM2NAME=term_object$term2name,
    pvalueCutoff=1,
    pAdjustMethod="BH",
    qvalueCutoff=1,
    minGSSize=5,
    maxGSSize=500
  )

  if (
    is.null(ora) ||
    nrow(as.data.frame(ora)) == 0
  ) {
    return(data.frame())
  }

  output <- as.data.frame(ora)
  output$analysis <- analysis_name
  output$database <- database_name
  output
}

ora_results <- bind_rows(
  run_ora(
    primary$result,
    go_terms,
    "Primary",
    "GO_BP"
  ),
  run_ora(
    primary$result,
    kegg_terms,
    "Primary",
    "KEGG"
  ),
  run_ora(
    core$result,
    go_terms,
    "Core",
    "GO_BP"
  ),
  run_ora(
    core$result,
    kegg_terms,
    "Core",
    "KEGG"
  )
)

write.table(
  ora_results,
  file.path(tab_dir, "GO_KEGG_ORA.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat(
  "ORA rows:",
  nrow(ora_results),
  "\n"
)

cat("===== RANKED GO AND KEGG GSEA =====\n")

run_gsea <- function(
  result,
  term_object,
  analysis_name,
  database_name
) {
  gene_list <- result$rank_score
  names(gene_list) <- result$gene

  gene_list <- gene_list[
    is.finite(gene_list)
  ]

  gene_list <- sort(
    gene_list,
    decreasing=TRUE
  )

  gsea <- GSEA(
    geneList=gene_list,
    TERM2GENE=term_object$term2gene,
    TERM2NAME=term_object$term2name,
    minGSSize=10,
    maxGSSize=500,
    pvalueCutoff=1,
    pAdjustMethod="BH",
    verbose=FALSE,
    seed=TRUE,
    eps=0
  )

  result_table <- as.data.frame(gsea)

  if (nrow(result_table) > 0) {
    result_table$analysis <- analysis_name
    result_table$database <- database_name
  }

  list(
    object=gsea,
    table=result_table
  )
}

gsea_objects <- list(
  Primary_GO=run_gsea(
    primary$result,
    go_terms,
    "Primary",
    "GO_BP"
  ),
  Primary_KEGG=run_gsea(
    primary$result,
    kegg_terms,
    "Primary",
    "KEGG"
  ),
  Core_GO=run_gsea(
    core$result,
    go_terms,
    "Core",
    "GO_BP"
  ),
  Core_KEGG=run_gsea(
    core$result,
    kegg_terms,
    "Core",
    "KEGG"
  )
)

gsea_results <- bind_rows(
  lapply(
    gsea_objects,
    function(x) x$table
  )
)

write.table(
  gsea_results,
  file.path(tab_dir, "GO_KEGG_GSEA.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

safe_gseaplot2 <- function(
    plot_object,
    geneSetID,
    title = "",
    color = "#2166AC",
    ...
) {
    pathway_id <- as.character(geneSetID)[1]

    if (is.na(pathway_id) || !nzchar(pathway_id)) {
        stop("safe_gseaplot2: invalid geneSetID")
    }

    result_df <- as.data.frame(plot_object@result)

    if (!pathway_id %in% result_df$ID &&
        pathway_id %in% rownames(result_df)) {
        pathway_id <- as.character(
            result_df[pathway_id, "ID", drop = TRUE]
        )
    }

    if (!pathway_id %in% result_df$ID) {
        stop(
            "safe_gseaplot2: pathway ID not found: ",
            pathway_id
        )
    }

    gs <- enrichplot:::gsInfo(
        object = plot_object,
        geneSetID = pathway_id
    )

    required_columns <- c(
        "x",
        "runningScore",
        "position",
        "geneList"
    )

    missing_columns <- setdiff(required_columns, colnames(gs))

    if (length(missing_columns) > 0) {
        stop(
            "safe_gseaplot2: missing columns: ",
            paste(missing_columns, collapse = ", ")
        )
    }

    gs$x <- as.numeric(gs$x)
    gs$runningScore <- as.numeric(gs$runningScore)
    gs$position <- as.numeric(gs$position)
    gs$geneList <- as.numeric(gs$geneList)

    gs <- gs[
        is.finite(gs$x) &
        is.finite(gs$runningScore) &
        is.finite(gs$geneList),
        ,
        drop = FALSE
    ]

    if (nrow(gs) == 0) {
        stop(
            "safe_gseaplot2: no finite plotting data for ",
            pathway_id
        )
    }

    hits <- gs[
        !is.na(gs$position) & gs$position == 1,
        ,
        drop = FALSE
    ]

    if (nrow(hits) == 0) {
        stop(
            "safe_gseaplot2: no gene-set hits for ",
            pathway_id
        )
    }

    result_row <- result_df[
        match(pathway_id, result_df$ID),
        ,
        drop = FALSE
    ]

    nes_value <- if ("NES" %in% colnames(result_row)) {
        suppressWarnings(as.numeric(result_row$NES[1]))
    } else {
        NA_real_
    }

    fdr_value <- if ("p.adjust" %in% colnames(result_row)) {
        suppressWarnings(as.numeric(result_row$p.adjust[1]))
    } else {
        NA_real_
    }

    peak_index <- which.max(abs(gs$runningScore))
    peak_x <- gs$x[peak_index]
    peak_es <- gs$runningScore[peak_index]

    line_color <- as.character(color)[1]

    if (is.na(line_color) || !nzchar(line_color)) {
        line_color <- "#2166AC"
    }

    title_text <- as.character(title)[1]

    if (is.na(title_text) || !nzchar(title_text)) {
        if ("Description" %in% colnames(result_row)) {
            title_text <- as.character(result_row$Description[1])
        } else {
            title_text <- pathway_id
        }
    }

    stat_text <- paste0(
        if (is.finite(nes_value)) {
            sprintf("NES = %.2f", nes_value)
        } else {
            "NES = NA"
        },
        "; ",
        if (is.finite(fdr_value)) {
            sprintf("FDR = %.3g", fdr_value)
        } else {
            "FDR = NA"
        }
    )

    p_es <- ggplot2::ggplot(
        gs,
        ggplot2::aes(
            x = .data$x,
            y = .data$runningScore
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linewidth = 0.35,
            color = "grey60"
        ) +
        ggplot2::geom_line(
            linewidth = 0.9,
            color = line_color
        ) +
        ggplot2::geom_vline(
            xintercept = peak_x,
            linetype = "dashed",
            linewidth = 0.4,
            color = "grey45"
        ) +
        ggplot2::annotate(
            geom = "text",
            x = peak_x,
            y = peak_es,
            label = sprintf("ES = %.2f", peak_es),
            hjust = ifelse(
                peak_x < max(gs$x, na.rm = TRUE) / 2,
                -0.05,
                1.05
            ),
            vjust = ifelse(peak_es >= 0, -0.5, 1.5),
            size = 3.0
        ) +
        ggplot2::labs(
            title = title_text,
            subtitle = stat_text,
            x = NULL,
            y = "Running ES"
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(
                face = "bold",
                size = 11
            ),
            plot.subtitle = ggplot2::element_text(
                size = 9
            ),
            plot.margin = ggplot2::margin(
                t = 5,
                r = 6,
                b = 0,
                l = 6
            )
        )

    p_hits <- ggplot2::ggplot(
        hits,
        ggplot2::aes(x = .data$x)
    ) +
        ggplot2::geom_segment(
            ggplot2::aes(
                xend = .data$x,
                y = 0,
                yend = 1
            ),
            linewidth = 0.3,
            color = line_color
        ) +
        ggplot2::coord_cartesian(
            xlim = range(gs$x, na.rm = TRUE),
            ylim = c(0, 1),
            expand = FALSE
        ) +
        ggplot2::labs(
            x = NULL,
            y = "Hits"
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(
                t = 0,
                r = 6,
                b = 0,
                l = 6
            )
        )

    p_rank <- ggplot2::ggplot(
        gs,
        ggplot2::aes(
            x = .data$x,
            y = .data$geneList
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linewidth = 0.35,
            color = "grey60"
        ) +
        ggplot2::geom_area(
            ggplot2::aes(
                fill = .data$geneList >= 0
            ),
            alpha = 0.75,
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                "TRUE" = "#D73027",
                "FALSE" = "#4575B4"
            )
        ) +
        ggplot2::labs(
            x = "Rank in Ordered Dataset",
            y = "Rank Metric"
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(
                t = 0,
                r = 6,
                b = 5,
                l = 6
            )
        )

    patchwork::wrap_plots(
        list(p_es, p_hits, p_rank),
        ncol = 1,
        heights = c(3.0, 0.55, 1.4)
    )
}

make_gsea_plot <- function(
  gsea_result,
  prefix,
  title_text
) {
  gsea_object <- gsea_result$object
  result_table <- as.data.frame(gsea_object)

  if (nrow(result_table) == 0) {
    return(NULL)
  }

  selected <- result_table %>%
    filter(
      is.finite(NES),
      is.finite(pvalue)
    ) %>%
    arrange(
      p.adjust,
      pvalue,
      desc(abs(NES))
    ) %>%
    slice_head(n=5)

  if (nrow(selected) == 0) {
    return(NULL)
  }

  selected$code <- paste0(
    prefix,
    seq_len(nrow(selected))
  )

  pathway_map <- selected %>%
    transmute(
      code=code,
      pathway_id=ID,
      pathway_name=Description,
      NES=NES,
      p_value=pvalue,
      fdr_bh=p.adjust
    )

  write.table(
    pathway_map,
    file.path(
      tab_dir,
      paste0(
        title_text,
        "_pathway_legend.tsv"
      )
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

  plot_object <- gsea_object

  matched_rows <- match(
    selected$ID,
    plot_object@result$ID
  )

  plot_object@result$Description[
    matched_rows
  ] <- selected$code

  caption_text <- paste(
    paste0(
      pathway_map$code,
      "=",
      pathway_map$pathway_name
    ),
    collapse="; "
  )

  p <- safe_gseaplot2(
    plot_object,
    geneSetID=selected$ID,
    title=title_text,
    base_size=12,
    subplots=1:3,
    pvalue_table=FALSE,
    ES_geom="line"
  ) +
    plot_annotation(
      caption=str_wrap(
        caption_text,
        width=180
      )
    )

  ggsave(
    file.path(
      fig_dir,
      paste0(
        title_text,
        "_multicurve.png"
      )
    ),
    p,
    width=13,
    height=9,
    dpi=300
  )

  pathway_map
}

gsea_legends <- bind_rows(
  make_gsea_plot(
    gsea_objects$Core_GO,
    "G",
    "Core_GO_GSEA"
  ),
  make_gsea_plot(
    gsea_objects$Core_KEGG,
    "K",
    "Core_KEGG_GSEA"
  ),
  make_gsea_plot(
    gsea_objects$Primary_GO,
    "G",
    "Primary_GO_GSEA"
  ),
  make_gsea_plot(
    gsea_objects$Primary_KEGG,
    "K",
    "Primary_KEGG_GSEA"
  )
)

cat("===== SUBJECT-LEVEL GSVA =====\n")

make_gene_sets <- function(
  term2gene,
  available_genes
) {
  gene_sets <- split(
    term2gene$gene,
    term2gene$term
  )

  gene_sets <- lapply(
    gene_sets,
    function(x) {
      intersect(
        unique(x),
        available_genes
      )
    }
  )

  gene_sets[
    lengths(gene_sets) >= 10 &
    lengths(gene_sets) <= 500
  ]
}

run_gsva_matrix <- function(
  expression_matrix,
  gene_sets
) {
  expression_matrix <- as.matrix(
    expression_matrix
  )

  if (
    exists(
      "gsvaParam",
      where=asNamespace("GSVA"),
      inherits=FALSE
    )
  ) {
    parameter <- GSVA::gsvaParam(
      exprData=expression_matrix,
      geneSets=gene_sets,
      kcdf="Gaussian",
      minSize=10,
      maxSize=500
    )

    GSVA::gsva(
      parameter,
      verbose=FALSE,
      BPPARAM=SerialParam()
    )
  } else {
    GSVA::gsva(
      expression_matrix,
      gene_sets,
      method="gsva",
      kcdf="Gaussian",
      min.sz=10,
      max.sz=500,
      parallel.sz=1,
      verbose=FALSE
    )
  }
}

run_gsva_analysis <- function(
  pseudobulk_object,
  analysis_name
) {
  expression_matrix <- pseudobulk_object$log_cpm
  sample_meta <- pseudobulk_object$sample_meta
  design <- pseudobulk_object$design

  go_sets <- make_gene_sets(
    go_terms$term2gene,
    rownames(expression_matrix)
  )

  kegg_sets <- make_gene_sets(
    kegg_terms$term2gene,
    rownames(expression_matrix)
  )

  cat(
    analysis_name,
    "GO GSVA pathways:",
    length(go_sets),
    "\n"
  )

  cat(
    analysis_name,
    "KEGG GSVA pathways:",
    length(kegg_sets),
    "\n"
  )

  go_scores <- run_gsva_matrix(
    expression_matrix,
    go_sets
  )

  kegg_scores <- run_gsva_matrix(
    expression_matrix,
    kegg_sets
  )

  test_scores <- function(
    score_matrix,
    database_name
  ) {
    fit <- lmFit(
      score_matrix,
      design
    )

    fit <- eBayes(
      fit,
      trend=TRUE,
      robust=TRUE
    )

    result <- topTable(
      fit,
      coef="GroupCWP",
      number=Inf,
      sort.by="none"
    )

    result$pathway_id <- rownames(result)
    result$analysis <- analysis_name
    result$database <- database_name

    result[, c(
      "analysis",
      "database",
      "pathway_id",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val"
    )]
  }

  list(
    GO_scores=go_scores,
    KEGG_scores=kegg_scores,
    results=bind_rows(
      test_scores(
        go_scores,
        "GO_BP"
      ),
      test_scores(
        kegg_scores,
        "KEGG"
      )
    ),
    sample_meta=sample_meta
  )
}

primary_gsva <- run_gsva_analysis(
  primary,
  "Primary"
)

core_gsva <- run_gsva_analysis(
  core,
  "Core"
)

gsva_results <- bind_rows(
  primary_gsva$results,
  core_gsva$results
)

gsva_results <- gsva_results %>%
  left_join(
    bind_rows(
      go_terms$term2name %>%
        mutate(database="GO_BP"),
      kegg_terms$term2name %>%
        mutate(database="KEGG")
    ),
    by=c(
      "pathway_id"="term",
      "database"="database"
    )
  )

write.table(
  gsva_results,
  file.path(tab_dir, "GO_KEGG_GSVA_tests.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("===== GSVA MAIN FIGURES =====\n")

select_top_gsva <- function(
  result_table,
  database_name,
  number=4
) {
  result_table %>%
    filter(
      database == database_name,
      is.finite(P.Value)
    ) %>%
    arrange(
      adj.P.Val,
      P.Value,
      desc(abs(logFC))
    ) %>%
    slice_head(n=number)
}

top_core_gsva <- bind_rows(
  select_top_gsva(
    core_gsva$results,
    "GO_BP",
    4
  ),
  select_top_gsva(
    core_gsva$results,
    "KEGG",
    4
  )
)

top_core_gsva <- top_core_gsva %>%
  group_by(database) %>%
  mutate(
    code=ifelse(
      database == "GO_BP",
      paste0("G", row_number()),
      paste0("K", row_number())
    )
  ) %>%
  ungroup() %>%
  left_join(
    bind_rows(
      go_terms$term2name %>%
        mutate(database="GO_BP"),
      kegg_terms$term2name %>%
        mutate(database="KEGG")
    ),
    by=c(
      "pathway_id"="term",
      "database"="database"
    )
  )

write.table(
  top_core_gsva,
  file.path(
    tab_dir,
    "Core_GSVA_pathway_legend.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

get_score_rows <- function(
  selected_table,
  gsva_object
) {
  bind_rows(
    lapply(
      seq_len(nrow(selected_table)),
      function(i) {
        pathway_id <- selected_table$pathway_id[i]
        database <- selected_table$database[i]

        score_matrix <- if (
          database == "GO_BP"
        ) {
          gsva_object$GO_scores
        } else {
          gsva_object$KEGG_scores
        }

        data.frame(
          pathway_id=pathway_id,
          database=database,
          code=selected_table$code[i],
          Subject=colnames(score_matrix),
          Score=as.numeric(
            score_matrix[
              pathway_id,
            ]
          ),
          stringsAsFactors=FALSE
        )
      }
    )
  )
}

gsva_plot_data <- get_score_rows(
  top_core_gsva,
  core_gsva
) %>%
  left_join(
    core_gsva$sample_meta %>%
      transmute(
        Subject,
        Group,
        Series
      ),
    by="Subject"
  ) %>%
  group_by(pathway_id) %>%
  mutate(
    ZScore=as.numeric(
      scale(Score)
    )
  ) %>%
  ungroup()

caption_text <- paste(
  paste0(
    top_core_gsva$code,
    "=",
    top_core_gsva$name
  ),
  collapse="; "
)

gsva_plot_data$code <- factor(
  gsva_plot_data$code,
  levels=rev(
    top_core_gsva$code
  )
)

gsva_plot_data$Group <- factor(
  gsva_plot_data$Group,
  levels=c("Control", "CWP")
)

p_heat <- ggplot(
  gsva_plot_data,
  aes(
    x=Subject,
    y=code,
    fill=ZScore
  )
) +
  geom_tile(
    colour="white",
    linewidth=0.3
  ) +
  facet_grid(
    . ~ Group,
    scales="free_x",
    space="free_x"
  ) +
  scale_fill_gradient2(
    low="#2166AC",
    mid="white",
    high="#B2182B",
    midpoint=0
  ) +
  labs(
    title="Core GSVA",
    x=NULL,
    y=NULL,
    fill="Z",
    caption=str_wrap(
      caption_text,
      width=180
    )
  ) +
  theme_bw(base_size=11) +
  theme(
    panel.grid=element_blank(),
    strip.background=element_rect(
      fill="white"
    ),
    axis.text.x=element_text(
      angle=45,
      hjust=1
    )
  )

ggsave(
  file.path(
    fig_dir,
    "Core_GSVA_heatmap.png"
  ),
  p_heat,
  width=13,
  height=6,
  dpi=300
)

p_box <- ggplot(
  gsva_plot_data,
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
    size=2.1
  ) +
  facet_wrap(
    ~ code,
    scales="free_y",
    ncol=4
  ) +
  scale_colour_manual(
    values=c(
      Control="#F8766D",
      CWP="#00BFC4"
    )
  ) +
  labs(
    title="Core GSVA pathways",
    x=NULL,
    y="GSVA",
    caption=str_wrap(
      caption_text,
      width=180
    )
  ) +
  theme_bw(base_size=11) +
  theme(
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
    "Core_GSVA_boxplots.png"
  ),
  p_box,
  width=14,
  height=8,
  dpi=300
)

saveRDS(
  list(
    primary=primary,
    core=core,
    GO_terms=go_terms,
    KEGG_terms=kegg_terms,
    GSEA=gsea_objects,
    Primary_GSVA=primary_gsva,
    Core_GSVA=core_gsva
  ),
  file.path(
    out_dir,
    "GSE288887_CD4_pathway_results.rds"
  )
)

cat("===== GSEA SUMMARY =====\n")

print(
  gsea_results %>%
    filter(
      analysis == "Core"
    ) %>%
    arrange(
      p.adjust,
      pvalue
    ) %>%
    select(
      database,
      ID,
      Description,
      NES,
      pvalue,
      p.adjust
    ) %>%
    head(20)
)

cat("===== GSVA SUMMARY =====\n")

print(
  gsva_results %>%
    filter(
      analysis == "Core"
    ) %>%
    arrange(
      adj.P.Val,
      P.Value
    ) %>%
    select(
      database,
      pathway_id,
      name,
      logFC,
      P.Value,
      adj.P.Val
    ) %>%
    head(20)
)

cat("===== FINAL OUTPUT =====\n")
cat("Tables:", tab_dir, "\n")
cat("Figures:", fig_dir, "\n")
cat("Results:", out_dir, "\n")
cat("STATUS=PASS\n")