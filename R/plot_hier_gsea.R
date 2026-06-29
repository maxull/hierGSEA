#' Plot a hierarchy-aware GSEA result
#'
#' The plot keeps the hierarchy visible directly in the graphic. The left side
#' is a compact parent-child tree strip, while the right side shows a dot plot
#' for enrichment strength and hierarchy-aware adjusted p-values. When the
#' result was processed with `directional = "both"`, the plot separates negative
#' and positive enrichment into `Down` and `Up` columns so the sign is retained
#' visually without sacrificing the shared hierarchy layout. The same function
#' also accepts a list of compatible `hier_gsea_result` objects and will then
#' draw one shared hierarchy with the comparison blocks arranged from left to
#' right in the order of the list.
#'
#' @param x A `hier_gsea_result` object returned by [hier_gsea()], or a
#'   non-empty list of compatible `hier_gsea_result` objects.
#' @param label_col Label column to display beside the hierarchy.
#' @param size_col Column used for point size.
#' @param colour_col Column used for point fill color.
#' @param show_left_hierarchy Logical. If `TRUE`, draw the hierarchy connector
#'   panel on the left.
#' @param tree_width Width of the left-side hierarchy strip in plot x-units.
#'   Increase this when the parent-child structure needs more horizontal room.
#' @param top_n_parents Optional number of top-ranked starting-level branches to
#'   show. This is applied after the hierarchy-aware ordering has already been
#'   computed, so `1` means "show the highest-ranked visible branch at
#'   `level_top`".
#' @param parent_terms Optional character vector, or list coercible to a
#'   character vector, containing starting-level term labels or IDs to show.
#'   This is useful when you want to focus the plot on specific major branches
#'   from the chosen `level_top`. If supplied, `top_n_parents` is ignored.
#' @param significance_cutoff Adjusted p-value threshold used for the white
#'   midpoint in the fill scale and for the optional outline marker. Defaults
#'   to `0.05`.
#'
#' @return A `ggplot2` object.
#' @export
plot_hier_gsea <- function(
    x,
    label_col = "Description",
    size_col = "abs_NES",
    colour_col = "p_adjust_hier",
    show_left_hierarchy = TRUE,
    tree_width = 0.32,
    top_n_parents = NULL,
    parent_terms = NULL,
    significance_cutoff = 0.05
) {
    if (!is.numeric(significance_cutoff) || length(significance_cutoff) != 1) {
        stop("significance_cutoff must be a single numeric value.", call. = FALSE)
    }

    if (significance_cutoff <= 0 || significance_cutoff >= 1) {
        stop("significance_cutoff must be greater than 0 and less than 1.", call. = FALSE)
    }

    if (!is.numeric(tree_width) || length(tree_width) != 1 || tree_width <= 0) {
        stop("tree_width must be a single positive numeric value.", call. = FALSE)
    }

    if (!is.null(top_n_parents)) {
        if (!is.numeric(top_n_parents) || length(top_n_parents) != 1 || top_n_parents < 1) {
            stop("top_n_parents must be a single positive numeric value.", call. = FALSE)
        }

        top_n_parents <- as.integer(top_n_parents)
    }

    if (is.list(parent_terms)) {
        parent_terms <- unlist(parent_terms, use.names = FALSE)
    }

    if (!is.null(parent_terms) && !is.character(parent_terms)) {
        stop("parent_terms must be NULL or a character vector.", call. = FALSE)
    }

    if (inherits(x, "hier_gsea_result")) {
        return(
            .plot_single_hier_gsea(
                x = x,
                label_col = label_col,
                size_col = size_col,
                colour_col = colour_col,
                show_left_hierarchy = show_left_hierarchy,
                tree_width = tree_width,
                top_n_parents = top_n_parents,
                parent_terms = parent_terms,
                significance_cutoff = significance_cutoff
            )
        )
    }

    if (!is.list(x) || length(x) == 0) {
        stop("x must be a 'hier_gsea_result' object or a non-empty list of them.", call. = FALSE)
    }

    if (!all(vapply(x, inherits, logical(1), what = "hier_gsea_result"))) {
        stop("Every element of x must inherit from 'hier_gsea_result'.", call. = FALSE)
    }

    .plot_multi_hier_gsea(
        x = x,
        label_col = label_col,
        size_col = size_col,
        colour_col = colour_col,
        show_left_hierarchy = show_left_hierarchy,
        tree_width = tree_width,
        top_n_parents = top_n_parents,
        parent_terms = parent_terms,
        significance_cutoff = significance_cutoff
    )
}

.plot_single_hier_gsea <- function(
    x,
    label_col,
    size_col,
    colour_col,
    show_left_hierarchy,
    tree_width,
    top_n_parents,
    parent_terms,
    significance_cutoff
) {
    plot_nodes <- x$plot_tbl$nodes

    if (!label_col %in% names(plot_nodes)) {
        stop("label_col was not found in the plot table.", call. = FALSE)
    }

    if (!size_col %in% names(plot_nodes)) {
        stop("size_col was not found in the plot table.", call. = FALSE)
    }

    if (!colour_col %in% names(plot_nodes)) {
        stop("colour_col was not found in the plot table.", call. = FALSE)
    }

    plot_layout <- .prepare_hier_gsea_plot_layout(
        plot_nodes = plot_nodes,
        directional = x$meta$directional,
        level_top = x$meta$level_top,
        label_col = label_col,
        tree_width = tree_width,
        top_n_parents = top_n_parents,
        parent_terms = parent_terms
    )

    plot_nodes <- plot_layout$nodes
    plot_edges <- plot_layout$edges
    plot_verticals <- plot_layout$verticals

    plot_nodes$plot_label <- plot_nodes[[label_col]]
    plot_nodes$plot_size <- plot_nodes[[size_col]]
    plot_nodes$plot_colour <- plot_nodes[[colour_col]]

    plot_nodes <- plot_nodes %>%
        dplyr::mutate(
            plot_label = stringr::str_wrap(.data$plot_label, width = 50),
            plot_size = dplyr::if_else(is.na(.data$plot_size), 0, .data$plot_size),
            plot_colour_raw = dplyr::if_else(is.na(.data$plot_colour), 1, .data$plot_colour),
            plot_colour_raw = pmin(pmax(.data$plot_colour_raw, 0), 1),
            plot_colour = .transform_hier_gsea_colour_value(
                values = .data$plot_colour_raw,
                significance_cutoff = significance_cutoff
            ),
            plot_label_face = dplyr::if_else(.data$level == x$meta$level_top, "bold", "plain"),
            outline_term = .data$term_in_input & .data$plot_colour_raw < significance_cutoff,
            show_point = .data$term_in_input & !is.na(.data$direction_nes)
        )

    point_breaks_raw <- sort(unique(c(0, significance_cutoff, 0.20, 0.50, 1)))
    point_breaks_raw <- point_breaks_raw[point_breaks_raw >= 0 & point_breaks_raw <= 1]
    point_breaks <- .transform_hier_gsea_colour_value(
        values = point_breaks_raw,
        significance_cutoff = significance_cutoff
    )

    horizontal_guide_tbl <- plot_nodes %>%
        dplyr::transmute(
            y = .data$plot_y,
            x = plot_layout$guide_x_min,
            xend = plot_layout$guide_x_max
        )

    vertical_guide_tbl <- tibble::tibble(
        x = plot_layout$dot_axis_breaks,
        xend = plot_layout$dot_axis_breaks,
        y = min(plot_nodes$plot_y, na.rm = TRUE),
        yend = max(plot_nodes$plot_y, na.rm = TRUE)
    )

    hierarchy_plot <- ggplot2::ggplot()

    if (nrow(horizontal_guide_tbl) > 0) {
        hierarchy_plot <- hierarchy_plot +
            ggplot2::geom_segment(
                data = horizontal_guide_tbl,
                ggplot2::aes(
                    x = .data$x,
                    xend = .data$xend,
                    y = .data$y,
                    yend = .data$y
                ),
                linewidth = 0.35,
                colour = "#ECECEC"
            )
    }

    if (nrow(vertical_guide_tbl) > 0) {
        hierarchy_plot <- hierarchy_plot +
            ggplot2::geom_segment(
                data = vertical_guide_tbl,
                ggplot2::aes(
                    x = .data$x,
                    xend = .data$xend,
                    y = .data$y,
                    yend = .data$yend
                ),
                linewidth = 0.35,
                colour = "#E3E3E3"
            )
    }

    if (show_left_hierarchy && nrow(plot_verticals) > 0) {
        hierarchy_plot <- hierarchy_plot +
            ggplot2::geom_segment(
                data = plot_verticals,
                ggplot2::aes(
                    x = .data$parent_x,
                    xend = .data$parent_x,
                    y = .data$y_start,
                    yend = .data$y_end
                ),
                linewidth = 0.35,
                colour = "grey20"
            )
    }

    if (show_left_hierarchy && nrow(plot_edges) > 0) {
        hierarchy_plot <- hierarchy_plot +
            ggplot2::geom_segment(
                data = plot_edges,
                ggplot2::aes(
                    x = .data$parent_x,
                    xend = .data$child_x,
                    y = .data$child_y,
                    yend = .data$child_y
                ),
                linewidth = 0.35,
                colour = "grey20"
            )
    }

    hierarchy_plot <- hierarchy_plot +
        ggplot2::geom_segment(
            data = plot_nodes,
            ggplot2::aes(
                x = .data$tree_x,
                xend = .data$label_x - 0.015,
                y = .data$plot_y,
                yend = .data$plot_y
            ),
            linewidth = 0.35,
            colour = "grey20"
        ) +
        ggplot2::geom_text(
            data = plot_nodes,
            ggplot2::aes(
                x = .data$label_x,
                y = .data$plot_y,
                label = .data$plot_label,
                fontface = .data$plot_label_face
            ),
            hjust = 0,
            size = 4,
            lineheight = 0.85
        ) +
        ggplot2::geom_point(
            data = plot_nodes %>%
                dplyr::filter(.data$show_point),
            ggplot2::aes(
                x = .data$dot_x,
                y = .data$plot_y,
                size = .data$plot_size,
                fill = .data$plot_colour
            ),
            shape = 21,
            colour = "grey25",
            stroke = 0.35
        ) +
        ggplot2::geom_point(
            data = plot_nodes %>%
                dplyr::filter(.data$outline_term),
            ggplot2::aes(
                x = .data$dot_x,
                y = .data$plot_y,
                size = .data$plot_size
            ),
            shape = 21,
            fill = NA,
            colour = "black",
            stroke = 0.75,
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_gradientn(
            colours = c("#440154FF", "#FFFFFF", "#A8A8A8"),
            values = c(0, 0.5, 1),
            limits = c(0, 1),
            breaks = point_breaks,
            labels = format(point_breaks_raw, digits = 2, nsmall = 2),
            oob = function(values, range = c(0, 1)) {
                pmin(pmax(values, range[[1]]), range[[2]])
            },
            name = "adj. p"
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_colorbar(
                barheight = grid::unit(7, "cm"),
                barwidth = grid::unit(0.7, "cm")
            )
        ) +
        ggplot2::scale_size_continuous(
            range = c(2.2, 10),
            name = "|NES|"
        ) +
        ggplot2::scale_x_continuous(
            breaks = plot_layout$dot_axis_breaks,
            labels = plot_layout$dot_axis_labels,
            expand = ggplot2::expansion(mult = c(0.01, 0.02))
        ) +
        ggplot2::scale_y_continuous(
            breaks = plot_nodes$plot_y,
            labels = rep("", nrow(plot_nodes)),
            expand = ggplot2::expansion(mult = c(0.01, 0.01))
        ) +
        ggplot2::coord_cartesian(
            xlim = c(plot_layout$x_min, plot_layout$x_max),
            clip = "off"
        ) +
        ggplot2::labs(
            title = paste0(
                .format_hier_gsea_db_label(x$meta$db),
                ifelse(is.na(x$meta$ontology), "", paste0(" ", x$meta$ontology)),
                " hierarchy-aware GSEA"
            ),
            subtitle = paste0(
                "Directional filter: ",
                x$meta$directional,
                " | Visible levels: ",
                x$meta$level_top,
                " to ",
                x$meta$level_bottom
            ),
            x = NULL,
            y = NULL
        ) +
        ggplot2::theme_classic(base_size = 18) +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white", colour = NA),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(face = "bold", colour = "black"),
            axis.title.x = ggplot2::element_text(face = "bold", colour = "black"),
            legend.title = ggplot2::element_text(face = "bold"),
            legend.text = ggplot2::element_text(face = "bold", colour = "black"),
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(colour = "black"),
            plot.margin = ggplot2::margin(1, 1.2, 1, 0.8, unit = "cm")
        )

    hierarchy_plot
}

.plot_multi_hier_gsea <- function(
    x,
    label_col,
    size_col,
    colour_col,
    show_left_hierarchy,
    tree_width,
    top_n_parents,
    parent_terms,
    significance_cutoff
) {
    .validate_multi_hier_gsea_results(x)

    reference_result <- x[[1]]
    comparison_labels <- names(x)

    if (is.null(comparison_labels) || any(!nzchar(comparison_labels))) {
        comparison_labels <- paste("Result", seq_along(x))
    }

    comparison_labels <- make.unique(comparison_labels, sep = " ")

    shared_plot_nodes <- .build_shared_plot_nodes(
        x = x,
        label_col = label_col
    )

    if (!size_col %in% names(reference_result$results_tbl)) {
        stop("size_col was not found in the hierarchy result tables.", call. = FALSE)
    }

    if (!colour_col %in% names(reference_result$results_tbl)) {
        stop("colour_col was not found in the hierarchy result tables.", call. = FALSE)
    }

    shared_layout <- .prepare_hier_gsea_plot_layout(
        plot_nodes = shared_plot_nodes,
        directional = reference_result$meta$directional,
        level_top = reference_result$meta$level_top,
        label_col = label_col,
        tree_width = tree_width,
        top_n_parents = top_n_parents,
        parent_terms = parent_terms
    )

    tree_nodes <- shared_layout$nodes
    tree_edges <- shared_layout$edges
    tree_verticals <- shared_layout$verticals

    tree_nodes <- tree_nodes %>%
        dplyr::mutate(
            plot_label = stringr::str_wrap(.data[[label_col]], width = 50),
            plot_label_face = dplyr::if_else(.data$level == reference_result$meta$level_top, "bold", "plain")
        )

    point_tbl <- .build_multi_result_point_table(
        x = x,
        comparison_labels = comparison_labels,
        tree_nodes = tree_nodes,
        size_col = size_col,
        colour_col = colour_col,
        significance_cutoff = significance_cutoff,
        directional = reference_result$meta$directional
    )

    x_axis_layout <- .build_multi_result_x_axis_layout(
        comparison_labels = comparison_labels,
        directional = reference_result$meta$directional
    )

    point_tbl <- point_tbl %>%
        dplyr::left_join(
            x_axis_layout$point_positions,
            by = c("comparison_label", "direction_column")
        )

    horizontal_guide_tbl <- tree_nodes %>%
        dplyr::transmute(
            y = .data$plot_y,
            x = shared_layout$guide_x_min,
            xend = x_axis_layout$guide_x_max
        )

    vertical_guide_tbl <- tibble::tibble(
        x = x_axis_layout$axis_breaks,
        xend = x_axis_layout$axis_breaks,
        y = min(tree_nodes$plot_y, na.rm = TRUE),
        yend = max(tree_nodes$plot_y, na.rm = TRUE)
    )

    comparison_header_tbl <- x_axis_layout$comparison_headers %>%
        dplyr::mutate(
            y = max(tree_nodes$plot_y, na.rm = TRUE) + 1.1
        )

    point_breaks_raw <- sort(unique(c(0, significance_cutoff, 0.20, 0.50, 1)))
    point_breaks_raw <- point_breaks_raw[point_breaks_raw >= 0 & point_breaks_raw <= 1]
    point_breaks <- .transform_hier_gsea_colour_value(
        values = point_breaks_raw,
        significance_cutoff = significance_cutoff
    )

    hierarchy_plot <- ggplot2::ggplot()

    if (nrow(horizontal_guide_tbl) > 0) {
        hierarchy_plot <- hierarchy_plot +
            ggplot2::geom_segment(
                data = horizontal_guide_tbl,
                ggplot2::aes(
                    x = .data$x,
                    xend = .data$xend,
                    y = .data$y,
                    yend = .data$y
                ),
                linewidth = 0.35,
                colour = "#ECECEC"
            )
    }

    if (nrow(vertical_guide_tbl) > 0) {
        hierarchy_plot <- hierarchy_plot +
            ggplot2::geom_segment(
                data = vertical_guide_tbl,
                ggplot2::aes(
                    x = .data$x,
                    xend = .data$xend,
                    y = .data$y,
                    yend = .data$yend
                ),
                linewidth = 0.35,
                colour = "#E3E3E3"
            )
    }

    if (show_left_hierarchy && nrow(tree_verticals) > 0) {
        hierarchy_plot <- hierarchy_plot +
            ggplot2::geom_segment(
                data = tree_verticals,
                ggplot2::aes(
                    x = .data$parent_x,
                    xend = .data$parent_x,
                    y = .data$y_start,
                    yend = .data$y_end
                ),
                linewidth = 0.35,
                colour = "grey20"
            )
    }

    if (show_left_hierarchy && nrow(tree_edges) > 0) {
        hierarchy_plot <- hierarchy_plot +
            ggplot2::geom_segment(
                data = tree_edges,
                ggplot2::aes(
                    x = .data$parent_x,
                    xend = .data$child_x,
                    y = .data$child_y,
                    yend = .data$child_y
                ),
                linewidth = 0.35,
                colour = "grey20"
            )
    }

    hierarchy_plot <- hierarchy_plot +
        ggplot2::geom_segment(
            data = tree_nodes,
            ggplot2::aes(
                x = .data$tree_x,
                xend = .data$label_x - 0.015,
                y = .data$plot_y,
                yend = .data$plot_y
            ),
            linewidth = 0.35,
            colour = "grey20"
        ) +
        ggplot2::geom_text(
            data = tree_nodes,
            ggplot2::aes(
                x = .data$label_x,
                y = .data$plot_y,
                label = .data$plot_label,
                fontface = .data$plot_label_face
            ),
            hjust = 0,
            size = 4,
            lineheight = 0.85
        ) +
        ggplot2::geom_text(
            data = comparison_header_tbl,
            ggplot2::aes(
                x = .data$x_center,
                y = .data$y,
                label = .data$comparison_label
            ),
            fontface = "bold",
            size = 4.4
        ) +
        ggplot2::geom_point(
            data = point_tbl %>%
                dplyr::filter(.data$show_point),
            ggplot2::aes(
                x = .data$dot_x,
                y = .data$plot_y,
                size = .data$plot_size,
                fill = .data$plot_colour
            ),
            shape = 21,
            colour = "grey25",
            stroke = 0.35
        ) +
        ggplot2::geom_point(
            data = point_tbl %>%
                dplyr::filter(.data$outline_term),
            ggplot2::aes(
                x = .data$dot_x,
                y = .data$plot_y,
                size = .data$plot_size
            ),
            shape = 21,
            fill = NA,
            colour = "black",
            stroke = 0.75,
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_gradientn(
            colours = c("#440154FF", "#FFFFFF", "#A8A8A8"),
            values = c(0, 0.5, 1),
            limits = c(0, 1),
            breaks = point_breaks,
            labels = format(point_breaks_raw, digits = 2, nsmall = 2),
            oob = function(values, range = c(0, 1)) {
                pmin(pmax(values, range[[1]]), range[[2]])
            },
            name = "adj. p"
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_colorbar(
                barheight = grid::unit(7, "cm"),
                barwidth = grid::unit(0.7, "cm")
            )
        ) +
        ggplot2::scale_size_continuous(
            range = c(2.2, 10),
            name = "|NES|"
        ) +
        ggplot2::scale_x_continuous(
            breaks = x_axis_layout$axis_breaks,
            labels = x_axis_layout$axis_labels,
            expand = ggplot2::expansion(mult = c(0.01, 0.02))
        ) +
        ggplot2::scale_y_continuous(
            breaks = tree_nodes$plot_y,
            labels = rep("", nrow(tree_nodes)),
            expand = ggplot2::expansion(mult = c(0.01, 0.08))
        ) +
        ggplot2::coord_cartesian(
            xlim = c(shared_layout$x_min, x_axis_layout$x_max),
            clip = "off"
        ) +
        ggplot2::labs(
            title = paste0(
                .format_hier_gsea_db_label(reference_result$meta$db),
                ifelse(is.na(reference_result$meta$ontology), "", paste0(" ", reference_result$meta$ontology)),
                " hierarchy-aware GSEA"
            ),
            subtitle = paste0(
                "Directional filter: ",
                reference_result$meta$directional,
                " | Visible levels: ",
                reference_result$meta$level_top,
                " to ",
                reference_result$meta$level_bottom
            ),
            x = NULL,
            y = NULL
        ) +
        ggplot2::theme_classic(base_size = 18) +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white", colour = NA),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(face = "bold", colour = "black"),
            axis.title.x = ggplot2::element_text(face = "bold", colour = "black"),
            legend.title = ggplot2::element_text(face = "bold"),
            legend.text = ggplot2::element_text(face = "bold", colour = "black"),
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(colour = "black"),
            plot.margin = ggplot2::margin(1.2, 1.2, 1, 0.8, unit = "cm")
        )

    hierarchy_plot
}

.validate_multi_hier_gsea_results <- function(x) {
    reference_meta <- x[[1]]$meta
    required_fields <- c("db", "ontology", "directional", "level_top", "level_bottom")

    for (result_index in seq_along(x)) {
        current_meta <- x[[result_index]]$meta

        for (field_name in required_fields) {
            reference_value <- reference_meta[[field_name]]
            current_value <- current_meta[[field_name]]

            values_match <- identical(reference_value, current_value) ||
                (length(reference_value) == 1 &&
                    length(current_value) == 1 &&
                    is.na(reference_value) &&
                    is.na(current_value))

            if (!values_match) {
                stop(
                    paste0(
                        "All hier_gsea_result objects must share the same ",
                        field_name,
                        " to be plotted together."
                    ),
                    call. = FALSE
                )
            }
        }
    }

    invisible(TRUE)
}

.build_shared_plot_nodes <- function(x, label_col) {
    shared_tbl <- purrr::map_dfr(
        x,
        function(current_result) {
            current_result$results_tbl
        }
    )

    if (!label_col %in% names(shared_tbl)) {
        stop("label_col was not found in the hierarchy result tables.", call. = FALSE)
    }

    summarised_tbl <- shared_tbl %>%
        dplyr::group_by(.data$term_id) %>%
        dplyr::summarise(
            Description = {
                valid_values <- stats::na.omit(.data$Description)

                if (length(valid_values) == 0) NA_character_ else as.character(valid_values[[1]])
            },
            level = {
                valid_values <- stats::na.omit(.data$level)

                if (length(valid_values) == 0) NA_integer_ else as.integer(valid_values[[1]])
            },
            canonical_parent_id = {
                valid_values <- stats::na.omit(.data$canonical_parent_id)

                if (length(valid_values) == 0) NA_character_ else as.character(valid_values[[1]])
            },
            canonical_path = list(dplyr::first(.data$canonical_path)),
            canonical_path_string = {
                valid_values <- stats::na.omit(.data$canonical_path_string)

                if (length(valid_values) == 0) NA_character_ else as.character(valid_values[[1]])
            },
            term_in_input = any(.data$term_in_input, na.rm = TRUE),
            p_adjust_hier = {
                valid_values <- .data$p_adjust_hier[!is.na(.data$p_adjust_hier)]

                if (length(valid_values) == 0) {
                    NA_real_
                } else {
                    min(valid_values)
                }
            },
            abs_NES = {
                valid_values <- .data$abs_NES[!is.na(.data$abs_NES)]

                if (length(valid_values) == 0) {
                    NA_real_
                } else {
                    max(valid_values)
                }
            },
            NES = NA_real_,
            .groups = "drop"
        )

    if (identical(label_col, "Description")) {
        shared_tbl <- summarised_tbl
    } else {
        label_tbl <- shared_tbl %>%
            dplyr::group_by(.data$term_id) %>%
            dplyr::summarise(
                plot_label_value = {
                    valid_values <- stats::na.omit(.data[[label_col]])

                    if (length(valid_values) == 0) NA_character_ else as.character(valid_values[[1]])
                },
                .groups = "drop"
            )

        shared_tbl <- summarised_tbl %>%
            dplyr::left_join(label_tbl, by = "term_id")

        names(shared_tbl)[names(shared_tbl) == "plot_label_value"] <- label_col
    }

    shared_tbl <- .compute_branch_metrics(result_tbl = shared_tbl)

    ordered_ids <- .order_branch_terms(
        result_tbl = shared_tbl,
        level_top = x[[1]]$meta$level_top
    )

    shared_tbl %>%
        dplyr::mutate(order_index = match(.data$term_id, ordered_ids)) %>%
        dplyr::arrange(.data$order_index)
}

.build_multi_result_point_table <- function(
    x,
    comparison_labels,
    tree_nodes,
    size_col,
    colour_col,
    significance_cutoff,
    directional
) {
    tree_term_tbl <- tree_nodes %>%
        dplyr::select(
            term_id = .data$term_id,
            plot_y = .data$plot_y
        )

    point_tbl <- purrr::map2_dfr(
        x,
        comparison_labels,
        function(current_result, comparison_label) {
            current_scope_tbl <- current_result$meta$testing_scope_tbl

            if (!size_col %in% names(current_scope_tbl)) {
                stop("size_col was not found in the hierarchy testing scope tables.", call. = FALSE)
            }

            if (!colour_col %in% names(current_scope_tbl)) {
                stop("colour_col was not found in the hierarchy testing scope tables.", call. = FALSE)
            }

            tree_term_tbl %>%
                dplyr::left_join(current_scope_tbl, by = "term_id") %>%
                dplyr::mutate(
                    comparison_label = comparison_label,
                    term_in_input = dplyr::if_else(is.na(.data$term_in_input), FALSE, .data$term_in_input),
                    plot_size = .data[[size_col]],
                    plot_colour_raw = .data[[colour_col]],
                    plot_size = dplyr::if_else(is.na(.data$plot_size), 0, .data$plot_size),
                    plot_colour_raw = dplyr::if_else(is.na(.data$plot_colour_raw), 1, .data$plot_colour_raw),
                    plot_colour_raw = pmin(pmax(.data$plot_colour_raw, 0), 1),
                    plot_colour = .transform_hier_gsea_colour_value(
                        values = .data$plot_colour_raw,
                        significance_cutoff = significance_cutoff
                    ),
                    direction_nes = .data$NES,
                    direction_column = dplyr::case_when(
                        identical(directional, "both") &
                            !is.na(.data$direction_nes) &
                            .data$direction_nes < 0 ~ "Down",
                        identical(directional, "both") &
                            !is.na(.data$direction_nes) &
                            .data$direction_nes >= 0 ~ "Up",
                        identical(directional, "down") &
                            !is.na(.data$direction_nes) ~ "Down",
                        identical(directional, "up") &
                            !is.na(.data$direction_nes) ~ "Up",
                        TRUE ~ NA_character_
                    ),
                    outline_term = .data$term_in_input & .data$plot_colour_raw < significance_cutoff,
                    show_point = .data$term_in_input & !is.na(.data$direction_column)
                ) %>%
                dplyr::select(
                    .data$comparison_label,
                    .data$term_id,
                    .data$plot_y,
                    .data$plot_size,
                    .data$plot_colour_raw,
                    .data$plot_colour,
                    .data$direction_nes,
                    .data$direction_column,
                    .data$outline_term,
                    .data$show_point
                )
        }
    )

    point_tbl
}

.build_multi_result_x_axis_layout <- function(comparison_labels, directional) {
    if (identical(directional, "both")) {
        direction_levels <- c("Down", "Up")
        direction_offsets <- c("Down" = 0, "Up" = 0.50)
        block_width <- 0.50
    } else if (identical(directional, "down")) {
        direction_levels <- "Down"
        direction_offsets <- c("Down" = 0)
        block_width <- 0
    } else {
        direction_levels <- "Up"
        direction_offsets <- c("Up" = 0)
        block_width <- 0
    }

    first_block_start <- 1.85
    comparison_gap <- 0.55
    comparison_step <- block_width + comparison_gap

    comparison_tbl <- tibble::tibble(
        comparison_label = comparison_labels,
        comparison_index = seq_along(comparison_labels),
        block_start = first_block_start + (comparison_index - 1) * comparison_step
    )

    point_positions <- tibble::tibble(
        comparison_label = rep(comparison_tbl$comparison_label, each = length(direction_levels)),
        block_start = rep(comparison_tbl$block_start, each = length(direction_levels)),
        direction_column = rep(direction_levels, times = nrow(comparison_tbl))
    ) %>%
        dplyr::mutate(
            dot_x = .data$block_start + unname(direction_offsets[.data$direction_column])
        ) %>%
        dplyr::select(
            .data$comparison_label,
            .data$direction_column,
            .data$dot_x
        )

    comparison_headers <- comparison_tbl %>%
        dplyr::mutate(
            x_center = .data$block_start + block_width / 2
        ) %>%
        dplyr::select(
            .data$comparison_label,
            .data$x_center
        )

    axis_tbl <- point_positions %>%
        dplyr::arrange(match(.data$comparison_label, comparison_labels), .data$dot_x)

    list(
        point_positions = point_positions,
        comparison_headers = comparison_headers,
        axis_breaks = axis_tbl$dot_x,
        axis_labels = axis_tbl$direction_column,
        x_max = max(point_positions$dot_x) + 0.16,
        guide_x_max = max(point_positions$dot_x) + 0.02
    )
}

# The plot layout is rebuilt here rather than reusing the raw plotting
# coordinates directly, because the visual requirements are more opinionated
# than the analysis object itself. This keeps the tree strip compact and makes
# branch-level spacing easy to tune without touching the statistical output.
.prepare_hier_gsea_plot_layout <- function(
    plot_nodes,
    directional,
    level_top,
    label_col,
    tree_width,
    top_n_parents = NULL,
    parent_terms = NULL
) {
    plot_nodes <- plot_nodes %>%
        dplyr::arrange(.data$order_index)

    root_ids <- plot_nodes %>%
        dplyr::filter(.data$level == level_top) %>%
        dplyr::arrange(.data$order_index) %>%
        dplyr::pull(.data$term_id)

    plot_nodes <- plot_nodes %>%
        dplyr::mutate(
            visible_root_id = purrr::map_chr(
                .data$canonical_path,
                ~ {
                    root_match <- .x[.x %in% root_ids]

                    if (length(root_match) == 0) {
                        return(NA_character_)
                    }

                    root_match[[1]]
                }
            )
        )

    selected_root_ids <- root_ids

    if (!is.null(parent_terms)) {
        root_tbl <- plot_nodes %>%
            dplyr::filter(.data$level == level_top) %>%
            dplyr::mutate(
                root_term_id_lower = tolower(.data$term_id),
                root_label_lower = tolower(as.character(.data[[label_col]]))
            )

        requested_terms <- tolower(parent_terms)

        selected_root_tbl <- root_tbl %>%
            dplyr::filter(
                .data$root_term_id_lower %in% requested_terms |
                    .data$root_label_lower %in% requested_terms
            )

        if (nrow(selected_root_tbl) == 0) {
            stop(
                "None of the requested parent_terms matched a visible starting-level term label or ID.",
                call. = FALSE
            )
        }

        # Preserve the user-supplied branch order whenever manual parent terms
        # are requested. This keeps the final plot compact and predictable even
        # when the selected branches were far apart in the full hierarchy.
        selected_root_ids <- character()

        for (requested_term in requested_terms) {
            matched_ids <- selected_root_tbl %>%
                dplyr::filter(
                    .data$root_term_id_lower == requested_term |
                        .data$root_label_lower == requested_term
                ) %>%
                dplyr::arrange(.data$order_index) %>%
                dplyr::pull(.data$term_id)

            selected_root_ids <- c(selected_root_ids, matched_ids)
        }

        selected_root_ids <- unique(selected_root_ids)
    } else if (!is.null(top_n_parents)) {
        selected_root_ids <- utils::head(root_ids, top_n_parents)
    }

    plot_nodes <- plot_nodes %>%
        dplyr::filter(.data$visible_root_id %in% selected_root_ids)

    if (nrow(plot_nodes) == 0) {
        stop("No rows remained after applying the parent branch selection.", call. = FALSE)
    }

    root_ids <- plot_nodes %>%
        dplyr::filter(.data$level == level_top) %>%
        dplyr::arrange(.data$order_index) %>%
        dplyr::pull(.data$term_id)

    root_order_tbl <- tibble::tibble(
        visible_root_id = selected_root_ids,
        selected_root_rank = seq_along(selected_root_ids)
    )

    plot_nodes <- plot_nodes %>%
        dplyr::left_join(root_order_tbl, by = "visible_root_id") %>%
        dplyr::arrange(.data$selected_root_rank, .data$order_index)

    visible_node_ids <- plot_nodes$term_id

    branch_gap <- 0.8
    branch_position_tbl <- plot_nodes %>%
        dplyr::group_by(.data$selected_root_rank, .data$visible_root_id) %>%
        dplyr::mutate(branch_row_index = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::count(.data$selected_root_rank, .data$visible_root_id, name = "branch_n") %>%
        dplyr::arrange(.data$selected_root_rank) %>%
        dplyr::mutate(
            branch_offset = cumsum(dplyr::lag(.data$branch_n + branch_gap, default = 0))
        )

    total_plot_height <- sum(branch_position_tbl$branch_n) +
        branch_gap * max(nrow(branch_position_tbl) - 1, 0)

    plot_nodes <- plot_nodes %>%
        dplyr::group_by(.data$selected_root_rank, .data$visible_root_id) %>%
        dplyr::mutate(branch_row_index = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(
            branch_position_tbl,
            by = c("selected_root_rank", "visible_root_id")
        ) %>%
        dplyr::mutate(
            plot_y = total_plot_height - .data$branch_offset - .data$branch_row_index + 1
        )

    max_depth <- max(plot_nodes$level - level_top, na.rm = TRUE)
    max_depth <- ifelse(is.finite(max_depth), max_depth, 0)

    label_x <- 0.22
    tree_strip_width <- tree_width
    tree_step <- if (max_depth == 0) {
        0
    } else {
        tree_strip_width / max_depth
    }

    plot_nodes <- plot_nodes %>%
        dplyr::mutate(
            label_x = label_x,
            tree_x = label_x - tree_strip_width + (.data$level - level_top) * tree_step,
            visible_parent_id = dplyr::if_else(
                .data$canonical_parent_id %in% visible_node_ids,
                .data$canonical_parent_id,
                NA_character_
            )
        )

    if (identical(directional, "both")) {
        dot_positions <- c("Down" = 1.85, "Up" = 2.35)

        plot_nodes <- plot_nodes %>%
            dplyr::mutate(
                direction_nes = .data$NES,
                direction_column = dplyr::case_when(
                    !is.na(.data$direction_nes) & .data$direction_nes < 0 ~ "Down",
                    !is.na(.data$direction_nes) & .data$direction_nes >= 0 ~ "Up",
                    TRUE ~ NA_character_
                ),
                dot_x = dplyr::case_when(
                    .data$direction_column == "Down" ~ dot_positions[["Down"]],
                    .data$direction_column == "Up" ~ dot_positions[["Up"]],
                    TRUE ~ NA_real_
                )
            )
    } else if (identical(directional, "down")) {
        dot_positions <- c("Down" = 2.10)

        plot_nodes <- plot_nodes %>%
            dplyr::mutate(
                direction_nes = .data$NES,
                direction_column = dplyr::case_when(
                    !is.na(.data$direction_nes) ~ "Down",
                    TRUE ~ NA_character_
                ),
                dot_x = dplyr::case_when(
                    .data$direction_column == "Down" ~ dot_positions[["Down"]],
                    TRUE ~ NA_real_
                )
            )
    } else {
        dot_positions <- c("Up" = 2.10)

        plot_nodes <- plot_nodes %>%
            dplyr::mutate(
                direction_nes = .data$NES,
                direction_column = dplyr::case_when(
                    !is.na(.data$direction_nes) ~ "Up",
                    TRUE ~ NA_character_
                ),
                dot_x = dplyr::case_when(
                    .data$direction_column == "Up" ~ dot_positions[["Up"]],
                    TRUE ~ NA_real_
                )
            )
    }

    edge_tbl <- plot_nodes %>%
        dplyr::select(
            child_id = term_id,
            child_y = plot_y,
            child_x = tree_x,
            parent_id = visible_parent_id
        ) %>%
        dplyr::filter(!is.na(.data$parent_id)) %>%
        dplyr::left_join(
            plot_nodes %>%
                dplyr::select(parent_id = term_id, parent_y = plot_y, parent_x = tree_x),
            by = "parent_id"
        )

    vertical_tbl <- edge_tbl %>%
        dplyr::group_by(.data$parent_id, .data$parent_x, .data$parent_y) %>%
        dplyr::summarise(
            child_y_min = min(.data$child_y),
            child_y_max = max(.data$child_y),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            y_start = pmin(.data$parent_y, .data$child_y_min),
            y_end = pmax(.data$parent_y, .data$child_y_max)
        )

    list(
        nodes = plot_nodes,
        edges = edge_tbl,
        verticals = vertical_tbl,
        dot_axis_breaks = unname(dot_positions),
        dot_axis_labels = names(dot_positions),
        x_min = label_x - tree_strip_width - 0.14,
        x_max = max(unname(dot_positions)) + 0.16,
        guide_x_min = label_x - 0.02,
        guide_x_max = max(unname(dot_positions)) + 0.02
    )
}

.transform_hier_gsea_colour_value <- function(values, significance_cutoff) {
    dplyr::case_when(
        is.na(values) ~ NA_real_,
        values <= significance_cutoff ~ 0.5 * (values / significance_cutoff),
        TRUE ~ 0.5 + 0.5 * ((values - significance_cutoff) / (1 - significance_cutoff))
    )
}

.format_hier_gsea_db_label <- function(db) {
    if (identical(db, "mitocarta")) {
        return("MitoCarta")
    }

    stringr::str_to_title(db)
}
