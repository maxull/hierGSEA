test_that("hier_gsea validates level windows and direction choices", {
    example_results <- get_example_gsea_results()

    expect_error(
        hier_gsea(
            result = example_results$reactome,
            db = "reactome",
            directional = "sideways"
        ),
        "directional"
    )

    expect_error(
        hier_gsea(
            result = example_results$reactome,
            db = "reactome",
            level_top = 3,
            level_bottom = 3
        ),
        "level_bottom"
    )
})

test_that("hier_gsea respects directional filtering", {
    example_results <- get_example_gsea_results()

    reactome_up <- hier_gsea(
        result = example_results$reactome,
        db = "reactome",
        directional = "up",
        level_top = 1,
        level_bottom = 4,
        alpha = 0.01
    )

    reactome_down <- hier_gsea(
        result = example_results$reactome,
        db = "reactome",
        directional = "down",
        level_top = 1,
        level_bottom = 4,
        alpha = 0.01
    )

    expect_true(all(reactome_up$results_tbl$NES[reactome_up$results_tbl$term_in_input] > 0, na.rm = TRUE))
    expect_true(all(reactome_down$results_tbl$NES[reactome_down$results_tbl$term_in_input] < 0, na.rm = TRUE))
})

test_that("family-wise BH correction matches within-family recalculation", {
    example_results <- get_example_gsea_results()

    reactome_processed <- hier_gsea(
        result = example_results$reactome,
        db = "reactome",
        directional = "both",
        level_top = 1,
        level_bottom = 4,
        alpha = 0.01
    )

    family_check_tbl <- reactome_processed$meta$testing_scope_tbl %>%
        dplyr::filter(.data$term_in_input, !is.na(.data$pvalue)) %>%
        dplyr::group_by(.data$family_id) %>%
        dplyr::mutate(expected_p_adjust_hier = stats::p.adjust(.data$pvalue, method = "BH")) %>%
        dplyr::ungroup()

    expect_equal(
        family_check_tbl$p_adjust_hier,
        family_check_tbl$expected_p_adjust_hier,
        tolerance = 1e-12
    )
})

test_that("ancestor retention keeps branch structure above significant descendants", {
    example_results <- get_example_gsea_results()

    reactome_processed <- hier_gsea(
        result = example_results$reactome,
        db = "reactome",
        directional = "both",
        level_top = 1,
        level_bottom = 4,
        alpha = 0.01
    )

    non_significant_rows <- reactome_processed$results_tbl %>%
        dplyr::filter(!.data$is_significant_hier)

    expect_true(nrow(non_significant_rows) > 0)

    for (row_index in seq_len(nrow(non_significant_rows))) {
        current_path <- non_significant_rows$canonical_path[[row_index]]

        descendant_tbl <- reactome_processed$results_tbl %>%
            dplyr::filter(
                purrr::map_lgl(
                    .data$canonical_path,
                    ~ length(current_path) <= length(.x) &&
                        identical(.x[seq_along(current_path)], current_path)
                )
            )

        expect_true(any(descendant_tbl$is_significant_hier, na.rm = TRUE))
    }
})

test_that("ordered output places parents before descendants", {
    example_results <- get_example_gsea_results()

    go_processed <- hier_gsea(
        result = example_results$go_bp,
        db = "go",
        ontology = "BP",
        directional = "up",
        level_top = 1,
        level_bottom = 4,
        alpha = 0.01
    )

    parent_order_lookup <- go_processed$results_tbl$order_index
    names(parent_order_lookup) <- go_processed$results_tbl$term_id

    child_rows <- go_processed$results_tbl %>%
        dplyr::filter(!is.na(.data$parent_id))

    expect_true(all(parent_order_lookup[child_rows$parent_id] < child_rows$order_index))
})

test_that("plot_hier_gsea supports non-root level_top windows and starting-level term selection", {
    example_results <- get_example_gsea_results()

    reactome_processed <- hier_gsea(
        result = example_results$reactome,
        db = "reactome",
        directional = "both",
        level_top = 2,
        level_bottom = 4,
        alpha = 0.01
    )

    expect_no_warning({
        level_top_plot <- plot_hier_gsea(
            reactome_processed,
            tree_width = 0.32,
            top_n_parents = 3
        )

        ggplot2::ggplot_build(level_top_plot)
    })

    expect_no_warning({
        selected_plot <- plot_hier_gsea(
            reactome_processed,
            tree_width = 0.32,
            parent_terms = list("Innate Immune System", "Cytokine Signaling in Immune system")
        )

        ggplot2::ggplot_build(selected_plot)
    })
})

test_that("plot_hier_gsea supports list input with a shared tree", {
    mock_result <- get_example_mitocarta_mock_result()

    mitocarta_a <- hier_gsea(
        result = mock_result,
        db = "mitocarta",
        directional = "both",
        level_top = 1,
        level_bottom = 3,
        alpha = 0.05
    )

    mitocarta_b <- hier_gsea(
        result = mock_result,
        db = "mitocarta",
        directional = "both",
        level_top = 1,
        level_bottom = 3,
        alpha = 0.10
    )

    expect_no_error({
        multi_plot <- plot_hier_gsea(
            x = list(Placebo = mitocarta_a, Antihistamine = mitocarta_b),
            tree_width = 0.4,
            top_n_parents = 2
        )

        multi_build <- ggplot2::ggplot_build(multi_plot)

        expect_equal(
            multi_build$layout$panel_params[[1]]$x$get_labels(),
            c("Down", "Up", "Down", "Up")
        )
    })
})

test_that("plot_hier_gsea list input validates shared metadata", {
    mock_result <- get_example_mitocarta_mock_result()

    mitocarta_both <- hier_gsea(
        result = mock_result,
        db = "mitocarta",
        directional = "both",
        level_top = 1,
        level_bottom = 3,
        alpha = 0.05
    )

    mitocarta_up <- hier_gsea(
        result = mock_result,
        db = "mitocarta",
        directional = "up",
        level_top = 1,
        level_bottom = 3,
        alpha = 0.05
    )

    expect_error(
        plot_hier_gsea(
            x = list(mitocarta_both, mitocarta_up)
        ),
        "same directional"
    )
})
