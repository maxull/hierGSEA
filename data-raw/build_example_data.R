# The bundled example data are intentionally simple. We keep the chronological
# script focused on the hierarchy-aware GSEA workflow rather than on a large
# amount of unrelated preprocessing.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tibble)
})

build_single_fiber_bulk_rankings <- function() {
    single_fiber_bulk_dir <- "/Users/maxullrich/Documents/GitHub/single-fiber-exercise/data/bulk"

    post_path <- file.path(single_fiber_bulk_dir, "DE_post_vs_pre.csv")
    rec_path <- file.path(single_fiber_bulk_dir, "DE_rec_vs_pre.csv")

    if (!file.exists(post_path) || !file.exists(rec_path)) {
        stop(
            "single-fiber-exercise bulk differential expression files were not found in the expected location.",
            call. = FALSE
        )
    }

    read_bulk_table <- function(path) {
        readr::read_csv(path, show_col_types = FALSE) %>%
            dplyr::rename(gene = 1) %>%
            dplyr::select("gene", "log2FoldChange", "stat", "pvalue", "padj", "xiao") %>%
            dplyr::filter(!is.na(.data$gene), .data$gene != "", !is.na(.data$log2FoldChange)) %>%
            dplyr::distinct(.data$gene, .keep_all = TRUE)
    }

    single_fiber_bulk_rankings <- list(
        post_vs_pre = read_bulk_table(post_path),
        rec_vs_pre = read_bulk_table(rec_path)
    )

    single_fiber_bulk_rankings
}
