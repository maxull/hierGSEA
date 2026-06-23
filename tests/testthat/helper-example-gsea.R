suppressPackageStartupMessages({
    library(dplyr)
})

example_cache <- new.env(parent = emptyenv())

.get_example_extdata_path <- function(file_name) {
    installed_path <- system.file("extdata", file_name, package = "hierGSEA")

    if (nzchar(installed_path)) {
        return(installed_path)
    }

    testthat::test_path("..", "..", "inst", "extdata", file_name)
}

get_example_rankings <- function() {
    if (!exists("rankings", envir = example_cache, inherits = FALSE)) {
        post_tbl <- readr::read_csv(
            .get_example_extdata_path("single_fiber_bulk_post_vs_pre.csv"),
            show_col_types = FALSE
        )

        rec_tbl <- readr::read_csv(
            .get_example_extdata_path("single_fiber_bulk_rec_vs_pre.csv"),
            show_col_types = FALSE
        )

        example_cache$rankings <- list(
            post_vs_pre = post_tbl,
            rec_vs_pre = rec_tbl
        )
    }

    example_cache$rankings
}

convert_example_symbols_to_entrez <- function(ranking_df) {
    testthat::skip_if_not_installed("AnnotationDbi")
    testthat::skip_if_not_installed("org.Hs.eg.db")

    id_lookup <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = ranking_df$gene,
        keytype = "SYMBOL",
        columns = c("SYMBOL", "ENTREZID")
    ) %>%
        tibble::as_tibble() %>%
        dplyr::filter(!is.na(.data$ENTREZID))

    ranking_df %>%
        dplyr::inner_join(id_lookup, by = c("gene" = "SYMBOL")) %>%
        dplyr::mutate(abs_log2FoldChange = abs(.data$log2FoldChange)) %>%
        dplyr::arrange(dplyr::desc(.data$abs_log2FoldChange), .data$ENTREZID) %>%
        dplyr::distinct(.data$ENTREZID, .keep_all = TRUE) %>%
        dplyr::select(-abs_log2FoldChange)
}

build_ranked_vector <- function(ranking_df) {
    entrez_tbl <- convert_example_symbols_to_entrez(ranking_df)
    ranked_vector <- entrez_tbl$log2FoldChange
    names(ranked_vector) <- entrez_tbl$ENTREZID
    sort(ranked_vector, decreasing = TRUE)
}

build_symbol_ranked_vector <- function(ranking_df) {
    symbol_tbl <- ranking_df %>%
        dplyr::mutate(abs_log2FoldChange = abs(.data$log2FoldChange)) %>%
        dplyr::arrange(dplyr::desc(.data$abs_log2FoldChange), .data$gene) %>%
        dplyr::distinct(.data$gene, .keep_all = TRUE) %>%
        dplyr::select(-abs_log2FoldChange)

    ranked_vector <- symbol_tbl$log2FoldChange
    names(ranked_vector) <- symbol_tbl$gene
    sort(ranked_vector, decreasing = TRUE)
}

get_example_gsea_results <- function() {
    if (!exists("gsea_results", envir = example_cache, inherits = FALSE)) {
        testthat::skip_if_not_installed("clusterProfiler")
        testthat::skip_if_not_installed("ReactomePA")
        testthat::skip_if_not_installed("org.Hs.eg.db")

        rankings <- get_example_rankings()
        ranked_vector <- build_ranked_vector(rankings$post_vs_pre)

        example_cache$gsea_results <- list(
            reactome = ReactomePA::gsePathway(
                geneList = ranked_vector,
                organism = "human",
                minGSSize = 5,
                maxGSSize = 10000,
                pvalueCutoff = 1,
                pAdjustMethod = "none",
                verbose = FALSE
            ),
            go_bp = clusterProfiler::gseGO(
                geneList = ranked_vector,
                OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                minGSSize = 5,
                maxGSSize = 10000,
                pvalueCutoff = 1,
                pAdjustMethod = "none",
                verbose = FALSE
            )
        )
    }

    example_cache$gsea_results
}

get_example_mitocarta_mock_result <- function() {
    if (!exists("mitocarta_mock_result", envir = example_cache, inherits = FALSE)) {
        testthat::skip_if_not_installed("DOSE")
        loadNamespace("DOSE")

        rankings <- get_example_rankings()
        gene_list <- build_symbol_ranked_vector(rankings$post_vs_pre)
        term2gene_tbl <- hierGSEA::mitocarta_term2gene()
        term2name_tbl <- hierGSEA::mitocarta_term2name()
        backend_terms <- get("mitocarta_hierarchy_data", envir = asNamespace("hierGSEA"))$terms

        selected_terms <- backend_terms %>%
            dplyr::filter(.data$level >= 1, .data$level <= 4) %>%
            dplyr::arrange(.data$level, .data$term_name) %>%
            dplyr::slice_head(n = 24) %>%
            dplyr::left_join(term2name_tbl, by = "term_id") %>%
            dplyr::mutate(
                Description = dplyr::coalesce(.data$term_name.y, .data$term_name.x),
                NES = c(-2.4, -2.1, -1.8, -1.4, -1.1, -0.9, 0.8, 1.0, 1.2, 1.4, 1.6, 1.9,
                        -2.2, -1.7, -1.3, -1.0, 0.9, 1.1, 1.3, 1.5, 1.8, 2.0, 2.2, 2.5),
                enrichmentScore = .data$NES / 3,
                pvalue = c(0.0008, 0.0015, 0.004, 0.009, 0.018, 0.04, 0.06, 0.08,
                           0.012, 0.02, 0.03, 0.045, 0.0005, 0.003, 0.007, 0.015,
                           0.025, 0.035, 0.05, 0.07, 0.011, 0.017, 0.028, 0.042)
            )

        gene_sets <- term2gene_tbl %>%
            dplyr::filter(.data$term_id %in% selected_terms$term_id) %>%
            dplyr::group_by(.data$term_id) %>%
            dplyr::summarise(genes = list(.data$gene_symbol), .groups = "drop")

        result_tbl <- selected_terms %>%
            dplyr::left_join(gene_sets, by = "term_id") %>%
            dplyr::mutate(
                setSize = lengths(.data$genes),
                p.adjust = stats::p.adjust(.data$pvalue, method = "BH"),
                qvalues = .data$p.adjust,
                rank = seq_len(dplyr::n()),
                leading_edge = "tags=0%, list=0%, signal=0%",
                core_enrichment = vapply(
                    .data$genes,
                    function(gene_vector) {
                        paste(utils::head(gene_vector, 10), collapse = "/")
                    },
                    FUN.VALUE = character(1)
                )
            ) %>%
            dplyr::transmute(
                ID = .data$term_id,
                Description = .data$Description,
                setSize = .data$setSize,
                enrichmentScore = .data$enrichmentScore,
                NES = .data$NES,
                pvalue = .data$pvalue,
                p.adjust = .data$p.adjust,
                qvalues = .data$qvalues,
                rank = .data$rank,
                leading_edge = .data$leading_edge,
                core_enrichment = .data$core_enrichment
            )

        term2gene_subset <- term2gene_tbl %>%
            dplyr::filter(.data$term_id %in% result_tbl$ID)

        gene_sets_list <- split(term2gene_subset$gene_symbol, term2gene_subset$term_id)

        example_cache$mitocarta_mock_result <- methods::new(
            "gseaResult",
            result = as.data.frame(result_tbl),
            organism = "human",
            setType = "MitoCarta",
            geneSets = gene_sets_list,
            geneList = gene_list,
            keytype = "SYMBOL",
            permScores = matrix(numeric(0), nrow = 0, ncol = 0),
            params = list(),
            gene2Symbol = character(0),
            readable = FALSE,
            termsim = matrix(numeric(0), nrow = 0, ncol = 0),
            method = "mock_fgsea",
            dr = list()
        )
    }

    example_cache$mitocarta_mock_result
}
