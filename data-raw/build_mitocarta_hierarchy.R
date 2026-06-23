# This script builds a hierarchical MitoCarta backend using the Broad
# Institute's MitoPathways3.0 resource. Unlike GO and Reactome, MitoCarta is
# not built into clusterProfiler, so we also prepare TERM2GENE and TERM2NAME
# mappings that analysts can feed into clusterProfiler::GSEA().

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
})

source(file.path("R", "hierarchy-backend.R"))

build_mitocarta_hierarchy_data <- function(
    source_html_path = NULL,
    mitocarta_release = "MitoCarta3.0 (2020)"
) {
    if (is.null(source_html_path)) {
        stop(
            "build_mitocarta_hierarchy_data() requires a local MitoPathways3.0 HTML file.",
            call. = FALSE
        )
    }

    if (!file.exists(source_html_path)) {
        stop("The supplied MitoPathways3.0 HTML file was not found.", call. = FALSE)
    }

    mitopathways_tbl <- .parse_mitopathways_html(source_html_path = source_html_path)

    path_id_lookup <- .build_mitopathways_id_lookup(mitopathways_tbl$hierarchy_path)

    mitopathways_terms <- mitopathways_tbl %>%
        dplyr::transmute(
            hierarchy_path = .data$hierarchy_path,
            term_id = unname(path_id_lookup[.data$hierarchy_path]),
            term_name = .data$term_name,
            genes = .data$genes
        )

    mitopathways_edges <- mitopathways_terms %>%
        dplyr::mutate(
            parent_path = vapply(
                strsplit(.data$hierarchy_path, " > ", fixed = TRUE),
                function(path_parts) {
                    if (length(path_parts) == 1) {
                        return(NA_character_)
                    }

                    paste(path_parts[-length(path_parts)], collapse = " > ")
                },
                FUN.VALUE = character(1)
            ),
            parent_id = unname(path_id_lookup[.data$parent_path]),
            child_id = .data$term_id,
            relation = "mitopathways_parent_child"
        ) %>%
        dplyr::filter(!is.na(.data$parent_id)) %>%
        dplyr::select(.data$parent_id, .data$child_id, .data$relation)

    mitopathways_term_tbl <- mitopathways_terms %>%
        dplyr::select(.data$term_id, .data$term_name) %>%
        dplyr::distinct(.data$term_id, .keep_all = TRUE) %>%
        dplyr::arrange(.data$term_id)

    mitocarta_hierarchy_data <- .prepare_hierarchy_backend(
        edges_tbl = mitopathways_edges,
        term_tbl = mitopathways_term_tbl,
        db = "mitocarta",
        ontology = NA_character_,
        metadata = list(
            mitocarta_release = mitocarta_release,
            source_html_path = source_html_path,
            built_on = as.character(Sys.time()),
            source_label = "MitoPathways3.0_human_html"
        )
    )

    mitocarta_term2gene_tbl <- mitopathways_terms %>%
        dplyr::transmute(
            term_id = .data$term_id,
            gene_symbol = .data$genes
        ) %>%
        tidyr::unnest_longer(.data$gene_symbol) %>%
        dplyr::filter(!is.na(.data$gene_symbol), .data$gene_symbol != "") %>%
        dplyr::distinct(.data$term_id, .data$gene_symbol)

    mitocarta_term2name_tbl <- mitopathways_terms %>%
        dplyr::transmute(
            term_id = .data$term_id,
            term_name = .data$term_name
        ) %>%
        dplyr::distinct(.data$term_id, .keep_all = TRUE)

    list(
        hierarchy_data = mitocarta_hierarchy_data,
        term2gene_tbl = mitocarta_term2gene_tbl,
        term2name_tbl = mitocarta_term2name_tbl
    )
}

.parse_mitopathways_html <- function(source_html_path) {
    html_lines <- readLines(source_html_path, warn = FALSE, encoding = "UTF-8")
    html_text <- paste(html_lines, collapse = "\n")
    row_matches <- gregexpr(
        "<tr><td>.*?</td><td>.*?</td><td>.*?</td></tr>",
        html_text,
        perl = TRUE
    )

    row_strings <- regmatches(html_text, row_matches)[[1]]

    if (length(row_strings) == 0) {
        stop("No MitoPathways3.0 table rows could be parsed from the HTML file.", call. = FALSE)
    }

    row_tbl_list <- lapply(row_strings, function(row_string) {
        cell_matches <- gregexpr("<td>.*?</td>", row_string, perl = TRUE)
        cell_strings <- regmatches(row_string, cell_matches)[[1]]

        if (length(cell_strings) != 3) {
            return(NULL)
        }

        cleaned_cells <- vapply(cell_strings, .clean_mitopathways_html_cell, FUN.VALUE = character(1))
        gene_symbols <- trimws(unlist(strsplit(cleaned_cells[[3]], ",", fixed = TRUE)))

        tibble::tibble(
            term_name = cleaned_cells[[1]],
            hierarchy_path = cleaned_cells[[2]],
            genes = list(gene_symbols)
        )
    })

    dplyr::bind_rows(row_tbl_list) %>%
        dplyr::filter(!is.na(.data$hierarchy_path), .data$hierarchy_path != "") %>%
        dplyr::distinct(.data$hierarchy_path, .keep_all = TRUE)
}

.clean_mitopathways_html_cell <- function(cell_string) {
    cleaned_string <- cell_string
    cleaned_string <- gsub("<br\\s*/?>", " ", cleaned_string, perl = TRUE)
    cleaned_string <- gsub("&nbsp;", " ", cleaned_string, fixed = TRUE)
    cleaned_string <- gsub("&amp;", "&", cleaned_string, fixed = TRUE)
    cleaned_string <- gsub("</?td>", "", cleaned_string, perl = TRUE)
    cleaned_string <- gsub("<[^>]+>", "", cleaned_string, perl = TRUE)
    cleaned_string <- gsub("\\s+", " ", cleaned_string, perl = TRUE)
    trimws(cleaned_string)
}

.build_mitopathways_id_lookup <- function(hierarchy_paths) {
    unique_paths <- unique(hierarchy_paths)
    path_ids <- vapply(unique_paths, .make_mitopathways_term_id, FUN.VALUE = character(1))
    names(path_ids) <- unique_paths
    path_ids
}

.make_mitopathways_term_id <- function(hierarchy_path) {
    path_slug <- tolower(hierarchy_path)
    path_slug <- gsub("[^a-z0-9]+", "_", path_slug, perl = TRUE)
    path_slug <- gsub("^_+|_+$", "", path_slug, perl = TRUE)
    paste0("MITOCARTA_", path_slug)
}
