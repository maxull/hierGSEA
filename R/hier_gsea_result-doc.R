#' Result object returned by `hier_gsea()`
#'
#' `hier_gsea()` returns a structured S3 object of class `hier_gsea_result`.
#' The object is designed to keep the original upstream `gseaResult`, the
#' hierarchy-aware result table, the full hierarchy metadata, and the plotting
#' coordinates together in one place so downstream analysis stays reproducible.
#'
#' @section Top-level elements:
#' \describe{
#'   \item{`result_raw`}{The original `DOSE::gseaResult` object supplied to
#'   [hier_gsea()].}
#'   \item{`results_tbl`}{The main hierarchy-aware output table after direction
#'   filtering, visible-level trimming, family-wise BH correction, ancestor
#'   retention, and branch ordering.}
#'   \item{`hierarchy_tbl`}{The backend parent-child edge table for the selected
#'   database. For GO this retains all available parent-child links, even when a
#'   single canonical plotting parent is chosen later.}
#'   \item{`paths_tbl`}{The backend term-path table used to recover root-to-node
#'   lineage, canonical paths, and term levels.}
#'   \item{`plot_tbl`}{A list containing plot-ready node and edge coordinates:
#'   `nodes`, `edges`, and `verticals`.}
#'   \item{`meta`}{A named list with analysis settings, recommended upstream
#'   GSEA settings, backend version metadata, and the pre-pruning testing scope
#'   table used for hierarchy-aware multiple-testing correction.}
#' }
#'
#' @section Key columns in `results_tbl`:
#' \describe{
#'   \item{`term_id`}{Database term identifier used to match the upstream GSEA
#'   result onto the selected hierarchy backend.}
#'   \item{`Description`}{Display label used for tables and plotting.}
#'   \item{`NES`}{Normalized enrichment score from the upstream GSEA result.}
#'   \item{`abs_NES`}{Absolute normalized enrichment score used by default for
#'   point size.}
#'   \item{`pvalue`}{Raw upstream GSEA p-value. This is the value used for the
#'   hierarchy-aware BH recalculation.}
#'   \item{`p.adjust`}{Incoming upstream global adjusted p-value, retained for
#'   reference only.}
#'   \item{`p_adjust_hier`}{Hierarchy-aware adjusted p-value recalculated within
#'   visible hierarchy families.}
#'   \item{`is_significant_hier`}{Logical flag indicating
#'   `p_adjust_hier < alpha`.}
#'   \item{`term_in_input`}{Logical flag marking whether the row came directly
#'   from the upstream GSEA result (`TRUE`) or was added back as a retained
#'   ancestor placeholder (`FALSE`).}
#'   \item{`level`}{Visible hierarchy level used by `hierGSEA`, after applying
#'   backend-specific ontology handling such as GO container-node removal.}
#'   \item{`parent_id`}{Canonical parent used for family assignment and plotting.}
#'   \item{`canonical_path`}{Canonical root-to-node path as a list column.}
#'   \item{`family_id`}{Family label used for the sibling-wise BH correction at
#'   the chosen visible level window.}
#'   \item{`family_n_tested`}{Number of tested terms in that family before
#'   branch pruning.}
#'   \item{`branch_best_p`}{Best descendant `p_adjust_hier` within the retained
#'   branch, used for ordering.}
#'   \item{`branch_best_nes`}{Largest descendant absolute `NES` within the
#'   retained branch, used as a tie-breaker for ordering.}
#'   \item{`order_index`}{Final branch-preserving display order where parents
#'   always appear before descendants.}
#' }
#'
#' @section Notes on interpretation:
#' The `results_tbl` may contain non-significant ancestor rows with
#' `term_in_input = FALSE` and `NES = NA`. These rows are intentionally retained
#' so that significant descendants can be displayed inside an intact hierarchy
#' rather than as disconnected leaf terms.
#'
#' @name hier_gsea_result
NULL
