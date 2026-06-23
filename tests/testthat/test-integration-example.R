test_that("Reactome and GO example workflows produce non-empty hierarchy outputs", {
    example_results <- get_example_gsea_results()

    reactome_processed <- hier_gsea(
        result = example_results$reactome,
        db = "reactome",
        directional = "both",
        level_top = 1,
        level_bottom = 4,
        alpha = 0.01
    )

    go_processed <- hier_gsea(
        result = example_results$go_bp,
        db = "go",
        ontology = "BP",
        directional = "up",
        level_top = 1,
        level_bottom = 4,
        alpha = 0.01
    )

    mitocarta_processed <- hier_gsea(
        result = get_example_mitocarta_mock_result(),
        db = "mitocarta",
        directional = "both",
        level_top = 1,
        level_bottom = 3,
        alpha = 0.05
    )

    expect_s3_class(reactome_processed, "hier_gsea_result")
    expect_s3_class(go_processed, "hier_gsea_result")
    expect_s3_class(mitocarta_processed, "hier_gsea_result")

    expect_true(nrow(reactome_processed$results_tbl) > 0)
    expect_true(nrow(go_processed$results_tbl) > 0)
    expect_true(nrow(mitocarta_processed$results_tbl) > 0)

    expect_true(nrow(reactome_processed$plot_tbl$nodes) > 0)
    expect_true(nrow(go_processed$plot_tbl$nodes) > 0)
    expect_true(nrow(mitocarta_processed$plot_tbl$nodes) > 0)
})

test_that("Hierarchy plot returns tree geometry and a ggplot object", {
    example_results <- get_example_gsea_results()

    reactome_processed <- hier_gsea(
        result = example_results$reactome,
        db = "reactome",
        directional = "both",
        level_top = 1,
        level_bottom = 4,
        alpha = 0.01
    )

    hierarchy_plot <- plot_hier_gsea(reactome_processed)

    expect_s3_class(hierarchy_plot, "ggplot")
    expect_true(nrow(reactome_processed$plot_tbl$edges) > 0)
})

test_that("Chronological example script runs end to end", {
    testthat::skip_if_not_installed("pkgload")

    script_path <- normalizePath(
        testthat::test_path("..", "..", "inst", "scripts", "run_single_fiber_bulk_example.R"),
        winslash = "/"
    )
    output_dir <- file.path(tempdir(), "hierGSEA_example_output")
    old_script_path <- Sys.getenv("HIERGSEA_SCRIPT_PATH", unset = NA_character_)
    old_run_mitocarta <- Sys.getenv("HIERGSEA_RUN_MITOCARTA_EXAMPLE", unset = NA_character_)

    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    on.exit(
        {
            if (is.na(old_script_path)) {
                Sys.unsetenv("HIERGSEA_SCRIPT_PATH")
            } else {
                Sys.setenv(HIERGSEA_SCRIPT_PATH = old_script_path)
            }

            if (is.na(old_run_mitocarta)) {
                Sys.unsetenv("HIERGSEA_RUN_MITOCARTA_EXAMPLE")
            } else {
                Sys.setenv(HIERGSEA_RUN_MITOCARTA_EXAMPLE = old_run_mitocarta)
            }
        },
        add = TRUE
    )

    setwd(tempdir())
    Sys.setenv(HIERGSEA_SCRIPT_PATH = script_path)
    Sys.setenv(HIERGSEA_RUN_MITOCARTA_EXAMPLE = "false")

    source(script_path, local = new.env(parent = globalenv()))

    expect_true(file.exists(file.path(output_dir, "single_fiber_post_vs_pre_reactome_hierarchy.pdf")))
    expect_true(file.exists(file.path(output_dir, "single_fiber_post_vs_pre_go_bp_hierarchy.pdf")))
})
