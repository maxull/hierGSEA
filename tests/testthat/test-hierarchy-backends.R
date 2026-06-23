test_that("Reactome backend keeps expected root and level structure", {
    expect_equal(reactome_hierarchy_data$db, "reactome")
    expect_equal(length(reactome_hierarchy_data$roots), 29)
    expect_true(max(reactome_hierarchy_data$terms$level, na.rm = TRUE) >= 10)
    expect_true(nrow(reactome_hierarchy_data$edges) > 0)
    expect_true(all(!is.na(reactome_hierarchy_data$paths$canonical_path_string)))
})

test_that("GO backends are ontology-specific and non-empty", {
    expect_equal(go_bp_hierarchy_data$db, "go")
    expect_equal(go_bp_hierarchy_data$ontology, "BP")
    expect_equal(go_mf_hierarchy_data$ontology, "MF")
    expect_equal(go_cc_hierarchy_data$ontology, "CC")

    expect_true(nrow(go_bp_hierarchy_data$edges) > 0)
    expect_true(nrow(go_mf_hierarchy_data$edges) > 0)
    expect_true(nrow(go_cc_hierarchy_data$edges) > 0)
    expect_false(any(go_bp_hierarchy_data$terms$term_id %in% c("all", "GO:0008150")))
    expect_false(any(go_mf_hierarchy_data$terms$term_id %in% c("all", "GO:0003674")))
    expect_false(any(go_cc_hierarchy_data$terms$term_id %in% c("all", "GO:0005575")))
    expect_true(length(go_bp_hierarchy_data$roots) > 1)
})

test_that("MitoCarta backend and custom GSEA mappings are available", {
    expect_equal(mitocarta_hierarchy_data$db, "mitocarta")
    expect_true(length(mitocarta_hierarchy_data$roots) > 0)
    expect_true(nrow(mitocarta_hierarchy_data$edges) > 0)
    expect_true(max(mitocarta_hierarchy_data$terms$level, na.rm = TRUE) >= 3)

    term2gene_tbl <- hierGSEA::mitocarta_term2gene()
    term2name_tbl <- hierGSEA::mitocarta_term2name()

    expect_true(nrow(term2gene_tbl) > 0)
    expect_true(nrow(term2name_tbl) > 0)
    expect_true(all(c("term_id", "gene_symbol") %in% names(term2gene_tbl)))
    expect_true(all(c("term_id", "term_name") %in% names(term2name_tbl)))
    expect_true(any(term2name_tbl$term_name == "Metabolism"))
    expect_true(all(grepl("^MITOCARTA_", term2name_tbl$term_id)))
})

test_that("multi-parent terms retain all edges and one canonical parent", {
    reactome_multi_parent_ids <- reactome_hierarchy_data$edges %>%
        dplyr::count(.data$child_id, name = "n_parents") %>%
        dplyr::filter(.data$n_parents > 1) %>%
        dplyr::pull(.data$child_id)

    expect_true(length(reactome_multi_parent_ids) > 0)

    canonical_parent_tbl <- reactome_hierarchy_data$terms %>%
        dplyr::filter(.data$term_id %in% reactome_multi_parent_ids)

    expect_true(all(lengths(canonical_parent_tbl$all_parent_ids) > 1))
    expect_true(all(!is.na(canonical_parent_tbl$canonical_parent_id)))
})
