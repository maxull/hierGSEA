# This script builds ontology-specific GO hierarchy backends from GO.db. GO is
# a DAG rather than a simple tree, so the runtime package stores all parent
# links internally and also records one canonical path per term for plotting and
# family-wise testing.

suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(dplyr)
    library(tibble)
    library(GO.db)
})

source(file.path("R", "hierarchy-backend.R"))

build_go_hierarchy_data <- function(ontology = c("BP", "MF", "CC")) {
    ontology <- match.arg(ontology)

    parent_map <- switch(
        ontology,
        BP = as.list(GOBPPARENTS),
        MF = as.list(GOMFPARENTS),
        CC = as.list(GOCCPARENTS)
    )

    ontology_container_ids <- switch(
        ontology,
        BP = c("all", "GO:0008150"),
        MF = c("all", "GO:0003674"),
        CC = c("all", "GO:0005575")
    )

    edge_list <- list()

    for (term_id in names(parent_map)) {
        parent_ids <- unname(parent_map[[term_id]])

        if (length(parent_ids) == 0 || all(is.na(parent_ids))) {
            next
        }

        relation_types <- names(parent_map[[term_id]])

        edge_list[[term_id]] <- tibble::tibble(
            parent_id = parent_ids,
            child_id = term_id,
            relation = relation_types
        )
    }

    go_edges <- dplyr::bind_rows(edge_list) %>%
        dplyr::filter(!is.na(.data$parent_id), !is.na(.data$child_id)) %>%
        dplyr::filter(
            !.data$parent_id %in% ontology_container_ids,
            !.data$child_id %in% ontology_container_ids
        ) %>%
        dplyr::distinct(.data$parent_id, .data$child_id, .keep_all = TRUE)

    go_node_ids <- unique(c(go_edges$parent_id, go_edges$child_id))
    go_term_names <- AnnotationDbi::Term(GOTERM[go_node_ids])

    go_terms <- tibble::tibble(
        term_id = go_node_ids,
        term_name = unname(go_term_names)
    ) %>%
        dplyr::mutate(term_name = ifelse(is.na(.data$term_name), .data$term_id, .data$term_name)) %>%
        dplyr::arrange(.data$term_id)

    go_hierarchy_data <- .prepare_hierarchy_backend(
        edges_tbl = go_edges,
        term_tbl = go_terms,
        db = "go",
        ontology = ontology,
        metadata = list(
            ontology = ontology,
            dropped_container_ids = ontology_container_ids,
            built_on = as.character(Sys.time()),
            go_db_version = as.character(utils::packageVersion("GO.db"))
        )
    )

    go_hierarchy_data
}

build_go_hierarchy_data_from_obo <- function(
    obo_path,
    ontology = c("BP", "MF", "CC"),
    go_release = NA_character_
) {
    ontology <- match.arg(ontology)

    namespace_lookup <- c(
        BP = "biological_process",
        MF = "molecular_function",
        CC = "cellular_component"
    )

    target_namespace <- unname(namespace_lookup[[ontology]])
    obo_lines <- readLines(obo_path, warn = FALSE, encoding = "UTF-8")
    term_start_indices <- which(obo_lines == "[Term]")

    if (length(term_start_indices) == 0) {
        stop("No [Term] blocks were found in the GO OBO file.", call. = FALSE)
    }

    term_end_indices <- c(term_start_indices[-1] - 1, length(obo_lines))
    term_records <- vector(mode = "list", length = length(term_start_indices))

    for (block_index in seq_along(term_start_indices)) {
        block_lines <- obo_lines[term_start_indices[[block_index]]:term_end_indices[[block_index]]]
        id_line <- block_lines[grepl("^id: GO:", block_lines)][1]

        if (length(id_line) == 0 || is.na(id_line)) {
            next
        }

        namespace_line <- block_lines[grepl("^namespace: ", block_lines)][1]

        if (length(namespace_line) == 0 || is.na(namespace_line)) {
            next
        }

        namespace_value <- sub("^namespace: ", "", namespace_line)

        if (!identical(namespace_value, target_namespace)) {
            next
        }

        if (any(block_lines == "is_obsolete: true")) {
            next
        }

        term_id <- sub("^id: ", "", id_line)
        term_name_line <- block_lines[grepl("^name: ", block_lines)][1]
        term_name <- if (length(term_name_line) == 0 || is.na(term_name_line)) {
            term_id
        } else {
            sub("^name: ", "", term_name_line)
        }

        is_a_parents <- sub(" !.*$", "", sub("^is_a: ", "", block_lines[grepl("^is_a: GO:", block_lines)]))
        relationship_lines <- block_lines[grepl("^relationship: ", block_lines)]

        relationship_tbl <- tibble::tibble(
            relation = character(),
            parent_id = character()
        )

        if (length(relationship_lines) > 0) {
            relationship_tbl <- tibble::tibble(
                relation = sub("^relationship: ([^ ]+) .*", "\\1", relationship_lines),
                parent_id = sub(" !.*$", "", sub("^relationship: [^ ]+ ", "", relationship_lines))
            ) %>%
                dplyr::filter(grepl("^GO:", .data$parent_id))
        }

        parent_tbl <- tibble::tibble(
            parent_id = c(is_a_parents, relationship_tbl$parent_id),
            relation = c(rep("is_a", length(is_a_parents)), relationship_tbl$relation)
        ) %>%
            dplyr::filter(!is.na(.data$parent_id), .data$parent_id != "") %>%
            dplyr::distinct(.data$parent_id, .keep_all = TRUE)

        term_records[[block_index]] <- list(
            term_id = term_id,
            term_name = term_name,
            parent_tbl = parent_tbl
        )
    }

    term_records <- Filter(Negate(is.null), term_records)

    go_terms <- tibble::tibble(
        term_id = purrr::map_chr(term_records, "term_id"),
        term_name = purrr::map_chr(term_records, "term_name")
    ) %>%
        dplyr::arrange(.data$term_id)

    go_term_ids <- go_terms$term_id

    go_edges <- purrr::map2_dfr(
        term_records,
        seq_along(term_records),
        function(term_record, term_index) {
            term_record$parent_tbl %>%
                dplyr::filter(.data$parent_id %in% go_term_ids) %>%
                dplyr::transmute(
                    parent_id = .data$parent_id,
                    child_id = term_record$term_id,
                    relation = .data$relation
                )
        }
    ) %>%
        dplyr::distinct(.data$parent_id, .data$child_id, .keep_all = TRUE)

    go_hierarchy_data <- .prepare_hierarchy_backend(
        edges_tbl = go_edges,
        term_tbl = go_terms,
        db = "go",
        ontology = ontology,
        metadata = list(
            ontology = ontology,
            source_label = "go_basic_obo",
            go_release = go_release,
            obo_path = obo_path,
            built_on = as.character(Sys.time())
        )
    )

    go_hierarchy_data
}
