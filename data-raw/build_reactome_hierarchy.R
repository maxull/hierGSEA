# This script builds the internal Reactome hierarchy backend that the package
# uses at runtime. The source files come from the existing CrosSys repository so
# the package can ship a stable local snapshot rather than relying on live
# downloads for every analysis run.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tibble)
})

source(file.path("R", "hierarchy-backend.R"))

build_reactome_hierarchy_data <- function() {
    build_reactome_hierarchy_data_from_files(
        reactome_relation_path = "/Users/maxullrich/Documents/GitHub/Sharples-Lab/CrosSys/data/ReactomePathwaysRelation.txt",
        reactome_pathway_path = "/Users/maxullrich/Documents/GitHub/Sharples-Lab/CrosSys/data/ReactomePathways.txt"
    )
}

build_reactome_hierarchy_data_from_files <- function(
    reactome_relation_path,
    reactome_pathway_path,
    reactome_release = NA_character_,
    source_label = "local_files"
) {

    if (!file.exists(reactome_relation_path) || !file.exists(reactome_pathway_path)) {
        stop(
            "Reactome source files were not found in the expected CrosSys data directory.",
            call. = FALSE
        )
    }

    reactome_edges <- readr::read_tsv(
        reactome_relation_path,
        col_names = c("parent_id", "child_id"),
        show_col_types = FALSE
    ) %>%
        dplyr::filter(
            grepl("^R-HSA-", .data$parent_id),
            grepl("^R-HSA-", .data$child_id)
        ) %>%
        dplyr::mutate(relation = "reactome_parent_child")

    reactome_terms <- readr::read_tsv(
        reactome_pathway_path,
        col_names = c("term_id", "term_name", "species"),
        show_col_types = FALSE
    ) %>%
        dplyr::filter(.data$species == "Homo sapiens") %>%
        dplyr::select(.data$term_id, .data$term_name)

    reactome_node_ids <- unique(c(reactome_edges$parent_id, reactome_edges$child_id))

    reactome_terms <- reactome_terms %>%
        dplyr::filter(.data$term_id %in% reactome_node_ids) %>%
        dplyr::arrange(.data$term_id)

    reactome_hierarchy_data <- .prepare_hierarchy_backend(
        edges_tbl = reactome_edges,
        term_tbl = reactome_terms,
        db = "reactome",
        ontology = NA_character_,
        metadata = list(
            reactome_release = reactome_release,
            source_label = source_label,
            source_relation_path = reactome_relation_path,
            source_pathway_path = reactome_pathway_path,
            built_on = as.character(Sys.time()),
            relation_file_mtime = as.character(file.info(reactome_relation_path)$mtime),
            pathway_file_mtime = as.character(file.info(reactome_pathway_path)$mtime)
        )
    )

    reactome_hierarchy_data
}
