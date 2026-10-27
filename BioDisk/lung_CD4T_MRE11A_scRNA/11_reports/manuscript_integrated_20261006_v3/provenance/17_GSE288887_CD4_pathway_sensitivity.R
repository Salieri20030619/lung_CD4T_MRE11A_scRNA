suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(clusterProfiler)
    library(msigdbr)
})

options(stringsAsFactors = FALSE)
set.seed(20260813)

project_root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"

table_dir <- file.path(
    project_root,
    "10_tables/CD4_final/pathway_sensitivity_v1"
)

figure_dir <- file.path(
    project_root,
    "09_figures/CD4_final/pathway_sensitivity_v1"
)

result_dir <- file.path(
    project_root,
    "07_DEG/CD4_pseudobulk_v1/pathway_sensitivity_v1"
)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cat("===== PATHWAY SENSITIVITY ANALYSIS =====\n")
cat("Project:", project_root, "\n")

pick_col <- function(nms, candidates) {
    idx <- match(tolower(candidates), tolower(nms))
    idx <- idx[!is.na(idx)]

    if (length(idx) == 0) {
        return(NA_character_)
    }

    nms[idx[1]]
}

safe_numeric <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
}

input_dirs <- c(
    file.path(
        project_root,
        "07_DEG/CD4_pseudobulk_v1/source_adjusted"
    ),
    file.path(
        project_root,
        "07_DEG/CD4_pseudobulk_v1/pathway_v1"
    ),
    file.path(
        project_root,
        "10_tables/CD4_final/pathway_v1"
    )
)

input_files <- unique(unlist(lapply(
    input_dirs,
    function(x) {
        if (!dir.exists(x)) {
            return(character())
        }

        list.files(
            x,
            pattern = "\\.(tsv|csv)$",
            recursive = TRUE,
            full.names = TRUE,
            ignore.case = TRUE
        )
    }
)))

if (length(input_files) == 0) {
    stop("No DEG/ranking tables found")
}

cat("Candidate input files:", length(input_files), "\n")

inspect_rank_file <- function(file_path) {
    dt <- tryCatch(
        data.table::fread(
            file_path,
            data.table = FALSE,
            check.names = FALSE,
            showProgress = FALSE
        ),
        error = function(e) NULL
    )

    if (is.null(dt) || nrow(dt) == 0) {
        return(list())
    }

    nms <- colnames(dt)

    gene_col <- pick_col(
        nms,
        c(
            "gene",
            "Gene",
            "symbol",
            "gene_symbol",
            "feature",
            "features"
        )
    )

    if (is.na(gene_col)) {
        return(list())
    }

    model_col <- pick_col(
        nms,
        c("model", "Model", "comparison_model")
    )

    analysis_col <- pick_col(
        nms,
        c("analysis", "Analysis", "dataset_type")
    )

    has_src_column <- any(grepl(
        "srcadj|sourceadj|source_adjust",
        nms,
        ignore.case = TRUE
    ))

    path_is_src <- grepl(
        "source_adjusted|sourceadj",
        file_path,
        ignore.case = TRUE
    )

    model_is_src <- FALSE

    if (!is.na(model_col)) {
        model_is_src <- any(grepl(
            "source",
            as.character(dt[[model_col]]),
            ignore.case = TRUE
        ))

        if (model_is_src) {
            dt <- dt[
                grepl(
                    "source",
                    as.character(dt[[model_col]]),
                    ignore.case = TRUE
                ),
                ,
                drop = FALSE
            ]
        }
    }

    source_qualified <- (
        has_src_column ||
        path_is_src ||
        model_is_src
    )

    if (!source_qualified || nrow(dt) == 0) {
        return(list())
    }

    src_stat_col <- pick_col(
        colnames(dt),
        c(
            "statistic_SrcAdj",
            "stat_SrcAdj",
            "t_SrcAdj",
            "QLF_SrcAdj",
            "F_SrcAdj",
            "LR_SrcAdj"
        )
    )

    effect_col <- pick_col(
        colnames(dt),
        c(
            "logFC_SrcAdj",
            "log2FC_SrcAdj",
            "estimate_SrcAdj",
            "logFC",
            "log2FC",
            "log2FoldChange",
            "estimate"
        )
    )

    generic_stat_col <- pick_col(
        colnames(dt),
        c(
            "statistic",
            "stat",
            "t",
            "QLF",
            "F",
            "LR"
        )
    )

    if (is.na(effect_col) && is.na(src_stat_col) &&
        is.na(generic_stat_col)) {
        return(list())
    }

    analysis_values <- character()

    if (!is.na(analysis_col)) {
        av <- unique(as.character(dt[[analysis_col]]))

        if (any(tolower(av) == "primary")) {
            analysis_values <- c(analysis_values, "Primary")
        }

        if (any(tolower(av) == "core")) {
            analysis_values <- c(analysis_values, "Core")
        }
    }

    if (length(analysis_values) == 0) {
        if (grepl("primary", basename(file_path), ignore.case = TRUE)) {
            analysis_values <- "Primary"
        }

        if (grepl("core", basename(file_path), ignore.case = TRUE)) {
            analysis_values <- c(analysis_values, "Core")
        }
    }

    analysis_values <- unique(analysis_values)

    if (length(analysis_values) == 0) {
        return(list())
    }

    output <- list()

    for (analysis_name in analysis_values) {
        sub <- dt

        if (!is.na(analysis_col)) {
            sub <- sub[
                tolower(as.character(sub[[analysis_col]])) ==
                    tolower(analysis_name),
                ,
                drop = FALSE
            ]
        }

        if (nrow(sub) == 0) {
            next
        }

        genes <- toupper(trimws(as.character(sub[[gene_col]])))
        effect <- if (!is.na(effect_col)) {
            safe_numeric(sub[[effect_col]])
        } else {
            rep(NA_real_, nrow(sub))
        }

        if (!is.na(src_stat_col)) {
            raw_stat <- safe_numeric(sub[[src_stat_col]])

            if (
                all(raw_stat[is.finite(raw_stat)] >= 0) &&
                any(is.finite(effect))
            ) {
                rank_value <- sign(effect) * sqrt(abs(raw_stat))
                rank_method <- paste0(
                    "signed_sqrt_",
                    src_stat_col
                )
            } else {
                rank_value <- raw_stat
                rank_method <- src_stat_col
            }
        } else if (
            path_is_src &&
            !is.na(generic_stat_col) &&
            any(is.finite(effect))
        ) {
            raw_stat <- safe_numeric(sub[[generic_stat_col]])

            if (all(raw_stat[is.finite(raw_stat)] >= 0)) {
                rank_value <- sign(effect) * sqrt(abs(raw_stat))
                rank_method <- paste0(
                    "signed_sqrt_",
                    generic_stat_col
                )
            } else {
                rank_value <- raw_stat
                rank_method <- generic_stat_col
            }
        } else if (!is.na(effect_col)) {
            rank_value <- effect
            rank_method <- effect_col
        } else {
            rank_value <- safe_numeric(sub[[generic_stat_col]])
            rank_method <- generic_stat_col
        }

        rank_df <- data.frame(
            gene = genes,
            rank = rank_value,
            stringsAsFactors = FALSE
        )

        rank_df <- rank_df[
            !is.na(rank_df$gene) &
            nzchar(rank_df$gene) &
            is.finite(rank_df$rank),
            ,
            drop = FALSE
        ]

        rank_df <- rank_df %>%
            group_by(gene) %>%
            slice_max(
                order_by = abs(rank),
                n = 1,
                with_ties = FALSE
            ) %>%
            ungroup() %>%
            as.data.frame()

        if (nrow(rank_df) < 1000) {
            next
        }

        method_priority <- ifelse(
            grepl("stat|sqrt|qlf|lr", rank_method, ignore.case = TRUE),
            100000,
            0
        )

        src_priority <- ifelse(
            grepl(
                "srcadj|source_adjusted",
                paste(file_path, rank_method),
                ignore.case = TRUE
            ),
            1000000,
            0
        )

        output[[length(output) + 1]] <- list(
            analysis = analysis_name,
            file = file_path,
            gene_col = gene_col,
            effect_col = effect_col,
            rank_method = rank_method,
            n_genes = nrow(rank_df),
            score = nrow(rank_df) +
                method_priority +
                src_priority,
            ranks = rank_df
        )
    }

    output
}

rank_candidates <- unlist(
    lapply(input_files, inspect_rank_file),
    recursive = FALSE
)

if (length(rank_candidates) == 0) {
    stop(
        paste0(
            "No complete source-adjusted ranking table was found. ",
            "Sensitivity analysis was not run."
        )
    )
}

selected_rank <- list()

for (analysis_name in c("Primary", "Core")) {
    idx <- which(vapply(
        rank_candidates,
        function(x) identical(x$analysis, analysis_name),
        logical(1)
    ))

    if (length(idx) == 0) {
        stop(
            "No source-adjusted rank available for ",
            analysis_name
        )
    }

    candidate_scores <- vapply(
        rank_candidates[idx],
        function(x) x$score,
        numeric(1)
    )

    selected_rank[[analysis_name]] <-
        rank_candidates[[idx[which.max(candidate_scores)]]]
}

rank_manifest <- bind_rows(lapply(
    selected_rank,
    function(x) {
        data.frame(
            analysis = x$analysis,
            input_file = x$file,
            gene_column = x$gene_col,
            effect_column = ifelse(
                is.na(x$effect_col),
                "",
                x$effect_col
            ),
            rank_method = x$rank_method,
            retained_genes = x$n_genes,
            stringsAsFactors = FALSE
        )
    }
))

data.table::fwrite(
    rank_manifest,
    file.path(table_dir, "rank_input_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

cat("===== SELECTED SOURCE-ADJUSTED RANKS =====\n")
print(rank_manifest)

is_cyto_rp <- function(gene) {
    grepl("^RPL", gene) |
        (
            grepl("^RPS", gene) &
            !grepl("^RPS6K", gene)
        )
}

is_mito_sensitive <- function(gene) {
    grepl("^MRPL", gene) |
        grepl("^MRPS", gene) |
        grepl("^MT-", gene)
}

scenario_order <- c(
    "Full",
    "NoCytoRP",
    "NoRP+MT"
)

scenario_ranks <- list()
exclusion_summary <- list()

for (analysis_name in names(selected_rank)) {
    rank_df <- selected_rank[[analysis_name]]$ranks

    cyto_flag <- is_cyto_rp(rank_df$gene)
    mito_flag <- is_mito_sensitive(rank_df$gene)

    scenario_ranks[[analysis_name]] <- list(
        Full = rank_df,
        NoCytoRP = rank_df[!cyto_flag, , drop = FALSE],
        `NoRP+MT` = rank_df[
            !cyto_flag & !mito_flag,
            ,
            drop = FALSE
        ]
    )

    exclusion_summary[[analysis_name]] <- data.frame(
        analysis = analysis_name,
        scenario = scenario_order,
        input_genes = nrow(rank_df),
        excluded_genes = c(
            0,
            sum(cyto_flag),
            sum(cyto_flag | mito_flag)
        ),
        retained_genes = c(
            nrow(rank_df),
            sum(!cyto_flag),
            sum(!cyto_flag & !mito_flag)
        ),
        stringsAsFactors = FALSE
    )
}

exclusion_summary <- bind_rows(exclusion_summary)

data.table::fwrite(
    exclusion_summary,
    file.path(table_dir, "gene_exclusion_summary.tsv"),
    sep = "\t",
    quote = FALSE
)

cat("===== GENE EXCLUSION SUMMARY =====\n")
print(exclusion_summary)

cat("===== LOAD LOCAL MSIGDB =====\n")

go_raw <- msigdbr::msigdbr(
    db_species = "HS",
    collection = "C5",
    subcollection = "GO:BP"
)

kegg_raw <- msigdbr::msigdbr(
    db_species = "HS",
    collection = "C2",
    subcollection = "CP:KEGG_LEGACY"
)

prepare_gene_sets <- function(x) {
    term2gene <- x %>%
        transmute(
            term = as.character(gs_name),
            gene = toupper(as.character(gene_symbol))
        ) %>%
        filter(
            !is.na(term),
            nzchar(term),
            !is.na(gene),
            nzchar(gene)
        ) %>%
        distinct()

    term2name <- term2gene %>%
        distinct(term) %>%
        mutate(
            name = term %>%
                str_replace("^GOBP_", "") %>%
                str_replace("^KEGG_", "") %>%
                str_replace_all("_", " ") %>%
                str_to_title()
        ) %>%
        select(term, name)

    list(
        term2gene = term2gene,
        term2name = term2name
    )
}

gene_sets <- list(
    GO_BP = prepare_gene_sets(go_raw),
    KEGG = prepare_gene_sets(kegg_raw)
)

cat("GO pathways:", n_distinct(go_raw$gs_name), "\n")
cat("KEGG pathways:", n_distinct(kegg_raw$gs_name), "\n")

run_one_gsea <- function(
    rank_df,
    database_name,
    gene_set_object,
    analysis_name,
    scenario_name
) {
    rank_vector <- rank_df$rank
    names(rank_vector) <- rank_df$gene

    rank_vector <- sort(
        rank_vector,
        decreasing = TRUE
    )

    rank_vector <- rank_vector[
        is.finite(rank_vector) &
        !duplicated(names(rank_vector))
    ]

    set.seed(20260813)

    gsea_object <- tryCatch(
        suppressWarnings(
            clusterProfiler::GSEA(
                geneList = rank_vector,
                TERM2GENE = gene_set_object$term2gene,
                TERM2NAME = gene_set_object$term2name,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                eps = 0,
                seed = TRUE,
                verbose = FALSE,
                by = "fgsea"
            )
        ),
        error = function(e) {
            message(
                "[GSEA_ERROR] ",
                analysis_name,
                " | ",
                scenario_name,
                " | ",
                database_name,
                " | ",
                conditionMessage(e)
            )
            NULL
        }
    )

    if (is.null(gsea_object)) {
        return(data.frame())
    }

    safe_scenario <- gsub(
        "[^A-Za-z0-9]+",
        "_",
        scenario_name
    )

    saveRDS(
        gsea_object,
        file.path(
            result_dir,
            paste0(
                analysis_name,
                "_",
                database_name,
                "_",
                safe_scenario,
                "_GSEA.rds"
            )
        )
    )

    result <- as.data.frame(gsea_object)

    if (nrow(result) == 0) {
        return(data.frame())
    }

    result %>%
        mutate(
            analysis = analysis_name,
            scenario = scenario_name,
            database = database_name,
            .before = 1
        )
}

cat("===== RUN SENSITIVITY GSEA =====\n")

all_results <- list()
result_index <- 0

for (analysis_name in c("Primary", "Core")) {
    for (scenario_name in scenario_order) {
        rank_df <- scenario_ranks[[analysis_name]][[scenario_name]]

        for (database_name in names(gene_sets)) {
            cat(
                "[RUN]",
                analysis_name,
                "|",
                scenario_name,
                "|",
                database_name,
                "| genes=",
                nrow(rank_df),
                "\n"
            )

            result_index <- result_index + 1

            all_results[[result_index]] <- run_one_gsea(
                rank_df = rank_df,
                database_name = database_name,
                gene_set_object = gene_sets[[database_name]],
                analysis_name = analysis_name,
                scenario_name = scenario_name
            )
        }
    }
}

all_gsea <- bind_rows(all_results)

if (nrow(all_gsea) == 0) {
    stop("All sensitivity GSEA analyses returned zero rows")
}

all_gsea$scenario <- factor(
    all_gsea$scenario,
    levels = scenario_order
)

data.table::fwrite(
    as.data.frame(all_gsea),
    file.path(table_dir, "all_sensitivity_GSEA.tsv"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

gsea_counts <- all_gsea %>%
    group_by(
        analysis,
        scenario,
        database
    ) %>%
    summarise(
        tested_pathways = n(),
        fdr_lt_005 = sum(p.adjust < 0.05, na.rm = TRUE),
        positive_fdr = sum(
            p.adjust < 0.05 & NES > 0,
            na.rm = TRUE
        ),
        negative_fdr = sum(
            p.adjust < 0.05 & NES < 0,
            na.rm = TRUE
        ),
        .groups = "drop"
    )

data.table::fwrite(
    gsea_counts,
    file.path(table_dir, "sensitivity_GSEA_summary.tsv"),
    sep = "\t",
    quote = FALSE
)

cat("===== GSEA COUNTS =====\n")
print(gsea_counts)

is_conservative <- function(pathway_id) {
    base_flag <- grepl(
        paste0(
            "RIBOSOME|RIBOSOMAL|",
            "CYTOPLASMIC_TRANSLATION|",
            "RRNA_PROCESSING|",
            "RIBOSOME_BIOGENESIS|",
            "TRANSLATIONAL_INITIATION|",
            "TRANSLATIONAL_ELONGATION|",
            "TRANSLATIONAL_TERMINATION"
        ),
        pathway_id,
        ignore.case = TRUE
    )

    mito_flag <- grepl(
        "MITOCHONDR",
        pathway_id,
        ignore.case = TRUE
    )

    base_flag & !mito_flag
}

all_gsea <- all_gsea %>%
    mutate(
        conservative_pathway = is_conservative(ID)
    )

conservative_results <- all_gsea %>%
    filter(conservative_pathway) %>%
    arrange(
        analysis,
        scenario,
        database,
        p.adjust
    )

data.table::fwrite(
    as.data.frame(conservative_results),
    file.path(table_dir, "ribosome_translation_audit.tsv"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

mito_patterns <- c(
    OXPHOS_K = "^KEGG_OXIDATIVE_PHOSPHORYLATION$",
    OXPHOS_G = "^GOBP_OXIDATIVE_PHOSPHORYLATION$",
    ETC = "^GOBP_ELECTRON_TRANSPORT_CHAIN$",
    ATP_ET = "^GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT$",
    PMF_ATP = "^GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS$",
    NADH_UQ = "^GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE$",
    PROTON_TM = "^GOBP_PROTON_TRANSMEMBRANE_TRANSPORT$"
)

get_short_name <- function(pathway_id) {
    output <- rep(NA_character_, length(pathway_id))

    for (nm in names(mito_patterns)) {
        hit <- grepl(
            mito_patterns[[nm]],
            pathway_id,
            ignore.case = TRUE
        )
        output[hit] <- nm
    }

    output
}

mito_results <- all_gsea %>%
    mutate(
        pathway_short = get_short_name(ID)
    ) %>%
    filter(!is.na(pathway_short)) %>%
    arrange(
        analysis,
        pathway_short,
        scenario
    )

data.table::fwrite(
    as.data.frame(mito_results),
    file.path(table_dir, "mitochondrial_pathway_sensitivity.tsv"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

mito_stability <- mito_results %>%
    select(
        analysis,
        database,
        ID,
        pathway_short,
        scenario,
        NES,
        p.adjust
    ) %>%
    pivot_wider(
        names_from = scenario,
        values_from = c(NES, p.adjust),
        names_sep = "_"
    ) %>%
    mutate(
        same_direction = (
            sign(NES_Full) == sign(NES_NoCytoRP) &
            sign(NES_Full) == sign(`NES_NoRP+MT`)
        ),
        no_cyto_ratio = abs(NES_NoCytoRP) /
            pmax(abs(NES_Full), 1e-8),
        no_rp_mt_ratio = abs(`NES_NoRP+MT`) /
            pmax(abs(NES_Full), 1e-8),
        all_fdr_lt_005 = (
            p.adjust_Full < 0.05 &
            p.adjust_NoCytoRP < 0.05 &
            `p.adjust_NoRP+MT` < 0.05
        ),
        stability = case_when(
            same_direction &
                all_fdr_lt_005 &
                no_cyto_ratio >= 0.70 &
                no_rp_mt_ratio >= 0.70 ~ "Stable",
            same_direction &
                p.adjust_Full < 0.05 ~ "Partial",
            TRUE ~ "Unstable"
        )
    )

data.table::fwrite(
    as.data.frame(mito_stability),
    file.path(table_dir, "mitochondrial_pathway_stability.tsv"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

cat("===== MITOCHONDRIAL STABILITY =====\n")
print(
    mito_stability %>%
        select(
            analysis,
            pathway_short,
            NES_Full,
            NES_NoCytoRP,
            `NES_NoRP+MT`,
            p.adjust_Full,
            p.adjust_NoCytoRP,
            `p.adjust_NoRP+MT`,
            stability
        )
)

classify_gene <- function(gene) {
    case_when(
        is_cyto_rp(gene) ~ "CytoRP",
        grepl("^MRPL|^MRPS", gene) ~ "MitoRP",
        grepl("^MT-", gene) ~ "MTcode",
        grepl(
            "^NDUF|^ATP5|^UQCR|^SDH|^COX|^TIMM|^TOMM",
            gene
        ) ~ "NuResp",
        TRUE ~ "Other"
    )
}

leading_edge <- mito_results %>%
    filter(
        scenario == "Full",
        !is.na(core_enrichment),
        nzchar(core_enrichment)
    ) %>%
    select(
        analysis,
        database,
        ID,
        pathway_short,
        NES,
        p.adjust,
        core_enrichment
    ) %>%
    separate_rows(
        core_enrichment,
        sep = "/"
    ) %>%
    rename(gene = core_enrichment) %>%
    mutate(
        gene = toupper(trimws(gene)),
        gene_class = classify_gene(gene)
    )

data.table::fwrite(
    as.data.frame(leading_edge),
    file.path(table_dir, "mitochondrial_leading_edge_genes.tsv"),
    sep = "\t",
    quote = FALSE
)

leading_edge_summary <- leading_edge %>%
    count(
        analysis,
        database,
        pathway_short,
        gene_class,
        name = "gene_n"
    ) %>%
    group_by(
        analysis,
        database,
        pathway_short
    ) %>%
    mutate(
        fraction = gene_n / sum(gene_n)
    ) %>%
    ungroup()

data.table::fwrite(
    as.data.frame(leading_edge_summary),
    file.path(table_dir, "mitochondrial_leading_edge_summary.tsv"),
    sep = "\t",
    quote = FALSE
)

if (nrow(mito_results) > 0) {
    p_mito <- ggplot(
        mito_results,
        aes(
            x = scenario,
            y = NES,
            group = 1
        )
    ) +
        geom_hline(
            yintercept = 0,
            linewidth = 0.35,
            color = "grey55"
        ) +
        geom_line(
            linewidth = 0.6,
            color = "grey45"
        ) +
        geom_point(
            aes(
                color = NES,
                shape = p.adjust < 0.05,
                size = -log10(
                    pmax(p.adjust, 1e-300)
                )
            )
        ) +
        facet_grid(
            analysis ~ pathway_short,
            scales = "free_y"
        ) +
        scale_color_gradient2(
            low = "#2166AC",
            mid = "white",
            high = "#B2182B",
            midpoint = 0,
            name = "NES"
        ) +
        scale_shape_manual(
            values = c(
                "TRUE" = 16,
                "FALSE" = 1
            ),
            name = "FDR<0.05"
        ) +
        scale_size_continuous(
            range = c(2, 5),
            name = "-log10FDR"
        ) +
        labs(
            title = "Mito Pathway Sensitivity",
            x = "Scenario",
            y = "NES",
            caption = paste0(
                "NoCytoRP: cytosolic ribosomal protein genes removed; ",
                "NoRP+MT: RPL/RPS, MRPL/MRPS and MT-* genes removed."
            )
        ) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(
                angle = 35,
                hjust = 1
            ),
            strip.text = element_text(
                face = "bold"
            ),
            plot.title = element_text(
                face = "bold"
            ),
            plot.caption = element_text(
                hjust = 0,
                size = 8
            )
        )

    ggsave(
        file.path(
            figure_dir,
            "GSE288887_CD4_mito_sensitivity.png"
        ),
        p_mito,
        width = 16,
        height = 7,
        dpi = 300
    )

    ggsave(
        file.path(
            figure_dir,
            "GSE288887_CD4_mito_sensitivity.pdf"
        ),
        p_mito,
        width = 16,
        height = 7
    )
}

if (nrow(leading_edge_summary) > 0) {
    p_le <- ggplot(
        leading_edge_summary,
        aes(
            x = pathway_short,
            y = fraction,
            fill = gene_class
        )
    ) +
        geom_col(
            width = 0.75
        ) +
        facet_wrap(
            ~analysis,
            ncol = 1
        ) +
        scale_y_continuous(
            labels = scales::percent_format(
                accuracy = 1
            )
        ) +
        scale_fill_manual(
            values = c(
                CytoRP = "#999999",
                MitoRP = "#CC79A7",
                MTcode = "#D55E00",
                NuResp = "#0072B2",
                Other = "#56B4E9"
            )
        ) +
        labs(
            title = "Leading-Edge Composition",
            x = "Pathway",
            y = "Gene Fraction",
            fill = "Gene Class",
            caption = paste0(
                "CytoRP: cytosolic ribosomal proteins; ",
                "MitoRP: mitochondrial ribosomal proteins; ",
                "MTcode: mtDNA-encoded genes; ",
                "NuResp: nuclear-encoded respiratory genes."
            )
        ) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(
                angle = 35,
                hjust = 1
            ),
            plot.title = element_text(
                face = "bold"
            ),
            plot.caption = element_text(
                hjust = 0,
                size = 8
            )
        )

    ggsave(
        file.path(
            figure_dir,
            "GSE288887_CD4_mito_leading_edge.png"
        ),
        p_le,
        width = 11,
        height = 8,
        dpi = 300
    )
}

filtered_full <- all_gsea %>%
    filter(
        scenario == "Full",
        !conservative_pathway,
        p.adjust < 0.05
    ) %>%
    mutate(
        direction = ifelse(
            NES > 0,
            "UpCWP",
            "DownCWP"
        )
    ) %>%
    group_by(
        analysis,
        database,
        direction
    ) %>%
    arrange(
        p.adjust,
        desc(abs(NES)),
        .by_group = TRUE
    ) %>%
    slice_head(n = 6) %>%
    ungroup()

if (nrow(filtered_full) > 0) {
    label_mapping <- filtered_full %>%
        distinct(database, ID, Description) %>%
        arrange(database, ID) %>%
        mutate(
            pathway_code = sprintf(
                "P%03d",
                row_number()
            )
        )

    data.table::fwrite(
        as.data.frame(label_mapping),
        file.path(table_dir, "pathway_label_mapping.tsv"),
        sep = "\t",
        quote = FALSE
    )

    filtered_plot <- filtered_full %>%
        left_join(
            label_mapping,
            by = c(
                "database",
                "ID",
                "Description"
            )
        ) %>%
        mutate(
            pathway_code = factor(
                pathway_code,
                levels = rev(unique(pathway_code))
            )
        )

    p_filtered <- ggplot(
        filtered_plot,
        aes(
            x = NES,
            y = pathway_code
        )
    ) +
        geom_vline(
            xintercept = 0,
            linewidth = 0.35,
            color = "grey55"
        ) +
        geom_segment(
            aes(
                x = 0,
                xend = NES,
                yend = pathway_code,
                color = NES
            ),
            linewidth = 0.8
        ) +
        geom_point(
            aes(
                color = NES,
                size = -log10(
                    pmax(p.adjust, 1e-300)
                )
            )
        ) +
        facet_grid(
            analysis ~ database,
            scales = "free_y",
            space = "free_y"
        ) +
        scale_color_gradient2(
            low = "#2166AC",
            mid = "white",
            high = "#B2182B",
            midpoint = 0,
            name = "NES"
        ) +
        scale_size_continuous(
            range = c(2, 5),
            name = "-log10FDR"
        ) +
        labs(
            title = "Filtered GSEA",
            x = "NES",
            y = "Pathway Code",
            caption = paste0(
                "Cytosolic ribosome and basic translation pathways ",
                "are excluded. Full labels are provided in ",
                "pathway_label_mapping.tsv."
            )
        ) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.minor = element_blank(),
            strip.text = element_text(
                face = "bold"
            ),
            plot.title = element_text(
                face = "bold"
            ),
            plot.caption = element_text(
                hjust = 0,
                size = 8
            )
        )

    ggsave(
        file.path(
            figure_dir,
            "GSE288887_CD4_filtered_GSEA.png"
        ),
        p_filtered,
        width = 12,
        height = 10,
        dpi = 300
    )
}

writeLines(
    c(
        "Pathway sensitivity analysis",
        "",
        "Full: complete source-adjusted ranking.",
        "NoCytoRP: RPL* and RPS* genes removed; RPS6K* retained.",
        "NoRP+MT: RPL/RPS, MRPL/MRPS and MT-* genes removed.",
        "",
        "Stable definition:",
        "1. Same NES direction in all three scenarios.",
        "2. FDR < 0.05 in all three scenarios.",
        "3. Absolute NES remains at least 70% of the Full result.",
        "",
        "Conservative cytosolic ribosome and translation pathways are retained",
        "in the complete result table but excluded from the filtered main figure."
    ),
    file.path(
        table_dir,
        "SENSITIVITY_METHODS.txt"
    )
)

cat("===== FINAL OUTPUT =====\n")
cat("Tables:", table_dir, "\n")
cat("Figures:", figure_dir, "\n")
cat("Results:", result_dir, "\n")
cat("STATUS=PASS\n")

sessionInfo()