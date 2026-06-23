#' Get the MitoCarta TERM2GENE mapping for custom GSEA
#'
#' Returns the human MitoPathways3.0 gene-set mapping in the two-column format
#' expected by `clusterProfiler::GSEA()` via the `TERM2GENE` argument.
#'
#' @return A tibble with columns `term_id` and `gene_symbol`.
#' @export
mitocarta_term2gene <- function() {
    tibble::as_tibble(mitocarta_term2gene_tbl)
}

#' Get the MitoCarta TERM2NAME mapping for custom GSEA
#'
#' Returns the human MitoPathways3.0 term-name mapping in the two-column format
#' expected by `clusterProfiler::GSEA()` via the `TERM2NAME` argument.
#'
#' @return A tibble with columns `term_id` and `term_name`.
#' @export
mitocarta_term2name <- function() {
    tibble::as_tibble(mitocarta_term2name_tbl)
}
