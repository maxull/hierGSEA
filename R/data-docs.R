#' Single-fiber bulk ranking example data
#'
#' A small bundled example list derived from the single-fiber exercise bulk
#' differential expression outputs. Each element is a tibble with the columns
#' needed to build ranked vectors for `gseGO()` and `gsePathway()` examples.
#'
#' The same source data are also written to `inst/extdata/` as CSV files so the
#' chronological example/debug script can run from a package installation.
#'
#' @format A named list with two tibbles:
#' \describe{
#'   \item{post_vs_pre}{Differential expression rankings for the POST vs PRE comparison.}
#'   \item{rec_vs_pre}{Differential expression rankings for the REC vs PRE comparison.}
#' }
"single_fiber_bulk_rankings"
