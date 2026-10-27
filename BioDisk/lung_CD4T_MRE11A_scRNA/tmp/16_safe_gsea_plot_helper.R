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
