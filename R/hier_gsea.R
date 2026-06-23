#' Process a GSEA result while respecting ontology hierarchy
#'
#' `hier_gsea()` takes a `gseaResult` object generated upstream by
#' `clusterProfiler` or `ReactomePA`, maps the tested terms onto a precomputed
#' hierarchy, performs family-wise Benjamini-Hochberg correction within sibling
#' families, retains the significant branches and their visible ancestors, and
#' returns a structured object ready for plotting or downstream inspection.
#'
#' The function is intentionally opinionated about the order of operations:
#' hierarchy mapping happens first, then directional filtering, then level
#' windowing, then family-wise correction. This keeps the statistical families
#' aligned with exactly what the analyst is choosing to inspect.
#'
#' To retain broad parent pathways such as top-level Reactome branches, users
#' should usually generate the upstream `gseaResult` with permissive size
#' filters, for example `minGSSize = 5`, `maxGSSize = 10000`,
#' `pvalueCutoff = 1`, and `pAdjustMethod = "none"`. `hierGSEA` recomputes the
#' hierarchy-aware multiple-testing adjustment from the raw `pvalue` column, so
#' the upstream global adjusted p-values are not used directly.
#'
#' For MitoCarta analyses, the upstream result should typically come from
#' `clusterProfiler::GSEA()` using `mitocarta_term2gene()` and
#' `mitocarta_term2name()` as the custom term mapping inputs.
#'
#' @param result A `DOSE::gseaResult` object.
#' @param db Database name. Must be `"reactome"`, `"go"`, or `"mitocarta"`.
#' @param ontology Required when `db = "go"`. Must be one of `"BP"`, `"MF"`,
#'   or `"CC"`. Ignored for `"reactome"` and `"mitocarta"`.
#' @param directional Direction filter. `"both"` keeps all terms, `"up"` keeps
#'   `NES > 0`, and `"down"` keeps `NES < 0`.
#' @param level_top Top hierarchy level to display.
#' @param level_bottom Bottom hierarchy level to display. Defaults to the
#'   deepest level available after hierarchy mapping and directional filtering.
#' @param alpha Significance threshold used for branch retention and outline
#'   labeling.
#' @param correction_method Currently only `"BH_family"` is supported.
#'
#' @return A `hier_gsea_result` object.
#' @export
hier_gsea <- function(
    result,
    db,
    ontology = NULL,
    directional = "both",
    level_top = 1,
    level_bottom = NULL,
    alpha = 0.05,
    correction_method = "BH_family"
) {
    if (!inherits(result, "gseaResult")) {
        stop("result must inherit from 'gseaResult'.", call. = FALSE)
    }

    db <- tolower(db)
    .validate_directional(directional = directional)

    if (identical(db, "mitopathways")) {
        db <- "mitocarta"
    }

    if (!db %in% c("reactome", "go", "mitocarta")) {
        stop("db must be one of 'reactome', 'go', or 'mitocarta'.", call. = FALSE)
    }

    if (identical(db, "go")) {
        if (is.null(ontology)) {
            stop("ontology must be supplied when db = 'go'.", call. = FALSE)
        }

        ontology <- toupper(ontology)

        if (!ontology %in% c("BP", "MF", "CC")) {
            stop("ontology must be one of 'BP', 'MF', or 'CC'.", call. = FALSE)
        }
    } else {
        ontology <- NA_character_
    }

    if (!identical(correction_method, "BH_family")) {
        stop("Only correction_method = 'BH_family' is supported in v1.", call. = FALSE)
    }

    backend <- .get_hierarchy_backend(db = db, ontology = ontology)

    result_tbl <- tibble::as_tibble(result@result)
    names(result_tbl)[names(result_tbl) == "ID"] <- "term_id"

    result_tbl <- result_tbl %>%
        dplyr::left_join(
            backend$terms %>% dplyr::transmute(
                term_id = .data$term_id,
                backend_term_name = .data$term_name,
                level = .data$level,
                level_label = .data$level_label,
                canonical_parent_id = .data$canonical_parent_id,
                canonical_path = .data$canonical_path,
                canonical_path_string = .data$canonical_path_string
            ),
            by = "term_id"
        ) %>%
        dplyr::mutate(
            Description = .coalesce_term_name(.data$Description, .data$backend_term_name),
            parent_id = .data$canonical_parent_id,
            abs_NES = abs(.data$NES),
            term_in_input = TRUE
        ) %>%
        dplyr::select(-"backend_term_name")

    result_tbl <- result_tbl %>%
        dplyr::filter(!is.na(.data$level))

    if (nrow(result_tbl) == 0) {
        stop(
            "No terms from the input GSEA result could be mapped to the selected hierarchy.",
            call. = FALSE
        )
    }

    if (identical(directional, "up")) {
        result_tbl <- result_tbl %>%
            dplyr::filter(.data$NES > 0)
    }

    if (identical(directional, "down")) {
        result_tbl <- result_tbl %>%
            dplyr::filter(.data$NES < 0)
    }

    if (nrow(result_tbl) == 0) {
        stop(
            "No terms remained after applying the directional filter.",
            call. = FALSE
        )
    }

    level_window <- .validate_level_window(
        level_top = level_top,
        level_bottom = level_bottom,
        available_levels = sort(unique(stats::na.omit(backend$terms$level)))
    )

    level_top <- level_window$level_top
    level_bottom <- level_window$level_bottom

    result_tbl <- result_tbl %>%
        dplyr::filter(.data$level >= level_top, .data$level <= level_bottom)

    should_warn_missing_starting_level <- !any(
        result_tbl$level == level_top & result_tbl$term_in_input,
        na.rm = TRUE
    )

    if (identical(db, "go") && identical(level_top, 1L)) {
        should_warn_missing_starting_level <- FALSE
    }

    if (should_warn_missing_starting_level) {
        warning(
            paste0(
                "No tested terms were present at level_top = ",
                level_top,
                ". This often means the upstream clusterProfiler/ReactomePA ",
                "run filtered out large or small parent terms, or that the ",
                "directional filter removed those starting-level pathways. ",
                "Rerun GSEA with minGSSize = 5, maxGSSize = 10000, ",
                "pvalueCutoff = 1, and pAdjustMethod = 'none' to preserve ",
                "broad hierarchy levels."
            ),
            call. = FALSE
        )
    }

    result_tbl <- result_tbl %>%
        dplyr::mutate(
            family_id = dplyr::if_else(
                .data$level == level_top,
                "__VISIBLE_ROOT__",
                .data$parent_id
            )
        ) %>%
        dplyr::group_by(.data$family_id) %>%
        dplyr::mutate(
            family_n_tested = dplyr::n(),
            p_adjust_hier = stats::p.adjust(.data$pvalue, method = "BH")
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            is_significant_hier = .data$p_adjust_hier < alpha
        )

    testing_scope_tbl <- result_tbl

    retained_ids <- result_tbl %>%
        dplyr::filter(.data$is_significant_hier) %>%
        dplyr::pull(.data$term_id) %>%
        unique()

    if (length(retained_ids) == 0) {
        warning(
            "No hierarchy-aware significant terms were detected at the requested alpha level. Returning the visible level window without branch pruning.",
            call. = FALSE
        )

        retained_ids <- result_tbl$term_id
    } else {
        parent_lookup <- backend$terms$canonical_parent_id
        names(parent_lookup) <- backend$terms$term_id

        ancestor_ids <- retained_ids

        for (term_id in retained_ids) {
            current_id <- term_id

            repeat {
                parent_id <- parent_lookup[[current_id]]

                if (is.null(parent_id) || is.na(parent_id) || !nzchar(parent_id)) {
                    break
                }

                parent_level <- backend$terms %>%
                    dplyr::filter(.data$term_id == parent_id) %>%
                    dplyr::pull(.data$level)

                if (length(parent_level) == 0 || is.na(parent_level) || parent_level < level_top) {
                    break
                }

                ancestor_ids <- c(ancestor_ids, parent_id)
                current_id <- parent_id
            }
        }

        retained_ids <- unique(ancestor_ids)
    }

    missing_ancestor_ids <- setdiff(retained_ids, result_tbl$term_id)

    if (length(missing_ancestor_ids) > 0) {
        ancestor_placeholder_tbl <- backend$terms %>%
            dplyr::filter(.data$term_id %in% missing_ancestor_ids) %>%
            dplyr::mutate(
                Description = .data$term_name,
                setSize = NA_integer_,
                enrichmentScore = NA_real_,
                NES = NA_real_,
                pvalue = NA_real_,
                p.adjust = NA_real_,
                qvalues = NA_real_,
                rank = NA_real_,
                leading_edge = NA_character_,
                core_enrichment = NA_character_,
                canonical_parent_id = .data$canonical_parent_id,
                parent_id = .data$canonical_parent_id,
                abs_NES = NA_real_,
                term_in_input = FALSE,
                family_id = dplyr::if_else(
                    .data$level == level_top,
                    "__VISIBLE_ROOT__",
                    .data$canonical_parent_id
                ),
                family_n_tested = NA_integer_,
                p_adjust_hier = NA_real_,
                is_significant_hier = FALSE
            ) %>%
            dplyr::select(
                term_id,
                Description,
                setSize,
                enrichmentScore,
                NES,
                pvalue,
                p.adjust,
                qvalues,
                rank,
                leading_edge,
                core_enrichment,
                level,
                level_label,
                canonical_parent_id,
                parent_id,
                canonical_path,
                canonical_path_string,
                abs_NES,
                term_in_input,
                family_id,
                family_n_tested,
                p_adjust_hier,
                is_significant_hier
            )

        result_tbl <- dplyr::bind_rows(result_tbl, ancestor_placeholder_tbl)
    }

    result_tbl <- result_tbl %>%
        dplyr::filter(.data$term_id %in% retained_ids) %>%
        dplyr::left_join(
            backend$terms %>% dplyr::transmute(
                term_id = .data$term_id,
                all_parent_ids = .data$all_parent_ids
            ),
            by = "term_id"
        )

    result_tbl <- .compute_branch_metrics(result_tbl = result_tbl)

    ordered_ids <- .order_branch_terms(result_tbl = result_tbl, level_top = level_top)

    result_tbl <- result_tbl %>%
        dplyr::mutate(order_index = match(.data$term_id, ordered_ids)) %>%
        dplyr::arrange(.data$order_index)

    plot_components <- .build_plot_table(
        result_tbl = result_tbl,
        level_top = level_top,
        level_bottom = level_bottom
    )

    output <- list(
        result_raw = result,
        results_tbl = result_tbl,
        hierarchy_tbl = backend$edges,
        paths_tbl = backend$paths,
        plot_tbl = plot_components,
        meta = list(
            db = db,
            ontology = ontology,
            directional = directional,
            level_top = level_top,
            level_bottom = level_bottom,
            alpha = alpha,
            correction_method = correction_method,
            recommended_upstream_gsea_args = list(
                minGSSize = 5,
                maxGSSize = 10000,
                pvalueCutoff = 1,
                pAdjustMethod = "none"
            ),
            testing_scope_tbl = testing_scope_tbl,
            backend_metadata = backend$metadata
        )
    )

    class(output) <- "hier_gsea_result"

    output
}
