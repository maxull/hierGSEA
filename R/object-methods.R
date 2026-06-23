#' Print a hierarchy-aware GSEA result
#'
#' This print method keeps the output short while still surfacing the most
#' important context needed during exploratory analysis.
#'
#' @param x A `hier_gsea_result` object.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.hier_gsea_result <- function(x, ...) {
    cat("hier_gsea_result\n")
    cat("  database:    ", x$meta$db, "\n", sep = "")

    if (!is.null(x$meta$ontology) && !is.na(x$meta$ontology)) {
        cat("  ontology:    ", x$meta$ontology, "\n", sep = "")
    }

    cat("  directional: ", x$meta$directional, "\n", sep = "")
    cat(
        "  levels:      ",
        x$meta$level_top,
        " to ",
        x$meta$level_bottom,
        "\n",
        sep = ""
    )
    cat("  alpha:       ", x$meta$alpha, "\n", sep = "")
    cat("  rows shown:  ", nrow(x$results_tbl), "\n", sep = "")
    cat(
        "  significant: ",
        sum(x$results_tbl$is_significant_hier, na.rm = TRUE),
        "\n",
        sep = ""
    )

    invisible(x)
}

#' Summary for a hierarchy-aware GSEA result
#'
#' @param object A `hier_gsea_result` object.
#' @param ... Unused.
#'
#' @return A named list with compact summary metrics.
#' @export
summary.hier_gsea_result <- function(object, ...) {
    list(
        db = object$meta$db,
        ontology = object$meta$ontology,
        directional = object$meta$directional,
        level_top = object$meta$level_top,
        level_bottom = object$meta$level_bottom,
        n_rows = nrow(object$results_tbl),
        n_significant = sum(object$results_tbl$is_significant_hier, na.rm = TRUE)
    )
}
