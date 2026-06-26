suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(clusterProfiler)
    library(dplyr)
    library(org.Hs.eg.db)
    library(readr)
    library(ReactomePA)
})

################################################################################################################################################
######################################################      LOAD hierGSEA LOCALLY OR FROM AN INSTALL      ######################################
################################################################################################################################################

# This script is meant to work in two common situations:
# 1. You are developing inside the package repository and want a quick debug run.
# 2. You have installed the package and want a reproducible worked example.
#
# The helper below first tries the installed package. If that is not available,
# it falls back to loading the source package from the repository that contains
# this script. This keeps the example practical during active development.

get_script_path <- function() {
    script_path_env <- Sys.getenv("HIERGSEA_SCRIPT_PATH", unset = "")

    if (nzchar(script_path_env) && file.exists(script_path_env)) {
        return(normalizePath(script_path_env, winslash = "/"))
    }

    argument_matches <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

    if (length(argument_matches) > 0) {
        return(normalizePath(sub("^--file=", "", argument_matches[[1]]), winslash = "/"))
    }

    if (!is.null(sys.frames()[[1]]$ofile)) {
        return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/"))
    }

    source_guess <- file.path(getwd(), "inst", "scripts", "run_single_fiber_bulk_example.R")

    if (file.exists(source_guess)) {
        return(normalizePath(source_guess, winslash = "/"))
    }

    stop("Could not determine the script path for the example workflow.", call. = FALSE)
}

load_hiergsea_package <- function() {
    if (requireNamespace("hierGSEA", quietly = TRUE)) {
        library(hierGSEA)
        return(invisible(TRUE))
    }

    script_path <- get_script_path()
    package_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")

    if (!requireNamespace("pkgload", quietly = TRUE)) {
        stop(
            "pkgload is required to run the example from the source repository when hierGSEA is not installed.",
            call. = FALSE
        )
    }

    pkgload::load_all(package_root, quiet = TRUE, export_all = FALSE)
    invisible(TRUE)
}

load_hiergsea_package()

################################################################################################################################################
######################################################      DEFINE EXAMPLE PATHS AND OUTPUT LOCATION      ######################################
################################################################################################################################################

# We write plots into the current working directory so the example behaves the
# same way from an installed package and from a development checkout.

script_path <- get_script_path()
package_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")

installed_extdata_dir <- system.file("extdata", package = "hierGSEA")
source_extdata_dir <- file.path(package_root, "inst", "extdata")

example_data_dir <- if (nzchar(installed_extdata_dir)) {
    installed_extdata_dir
} else {
    source_extdata_dir
}

example_output_dir <- file.path(getwd(), "hierGSEA_example_output")
dir.create(example_output_dir, showWarnings = FALSE, recursive = TRUE)

# The bundled example runs the MitoCarta custom GSEA workflow by default
# because that is useful for real package demonstrations. Test environments can
# switch it off to keep automated checks faster while still covering MitoCarta
# with a dedicated mocked integration test elsewhere in the package.

run_mitocarta_example <- identical(
    tolower(Sys.getenv("HIERGSEA_RUN_MITOCARTA_EXAMPLE", unset = "true")),
    "true"
)

################################################################################################################################################
######################################################      LOAD BUNDLED RANKING TABLES      ####################################################
################################################################################################################################################

# The package ships lightweight ranking tables rather than a full upstream
# project snapshot. This keeps the example reproducible while still using real
# biological data with meaningful hierarchy structure.

bulk_post_vs_pre <- readr::read_csv(
    file.path(example_data_dir, "single_fiber_bulk_post_vs_pre.csv"),
    show_col_types = FALSE
)

bulk_rec_vs_pre <- readr::read_csv(
    file.path(example_data_dir, "single_fiber_bulk_rec_vs_pre.csv"),
    show_col_types = FALSE
)

################################################################################################################################################
######################################################      CONVERT SYMBOLS TO ENTREZ IDS      ##################################################
################################################################################################################################################

# Reactome and GO GSEA functions work most reliably with Entrez identifiers.
# This helper keeps the conversion logic local to the example script, as
# requested, so the public package API stays focused on hierarchy-aware
# post-processing rather than on identifier wrangling.

convert_symbols_to_entrez <- function(ranking_df) {
    id_lookup <- AnnotationDbi::select(
        org.Hs.eg.db,
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

bulk_post_vs_pre_entrez <- convert_symbols_to_entrez(bulk_post_vs_pre)
bulk_rec_vs_pre_entrez <- convert_symbols_to_entrez(bulk_rec_vs_pre)

# MitoCarta uses a custom TERM2GENE backend rather than an organism annotation
# package, so we also keep a symbol-based ranked vector in parallel. Using gene
# symbols here keeps the example close to the official MitoPathways resource.

bulk_post_vs_pre_symbols <- bulk_post_vs_pre %>%
    dplyr::mutate(abs_log2FoldChange = abs(.data$log2FoldChange)) %>%
    dplyr::arrange(dplyr::desc(.data$abs_log2FoldChange), .data$gene) %>%
    dplyr::distinct(.data$gene, .keep_all = TRUE) %>%
    dplyr::select(-abs_log2FoldChange)

################################################################################################################################################
######################################################      BUILD LOG2FC-RANKED VECTORS      ####################################################
################################################################################################################################################

# The example intentionally uses log2 fold change because that was the agreed
# demonstration ranking metric for this first package version. The vectors are
# sorted decreasingly so positively enriched terms appear naturally on the upper
# end of the ranked list and negative enrichment is still captured by the GSEA
# statistic.

post_vs_pre_ranked_vector <- bulk_post_vs_pre_entrez$log2FoldChange
names(post_vs_pre_ranked_vector) <- bulk_post_vs_pre_entrez$ENTREZID
post_vs_pre_ranked_vector <- sort(post_vs_pre_ranked_vector, decreasing = TRUE)

rec_vs_pre_ranked_vector <- bulk_rec_vs_pre_entrez$log2FoldChange
names(rec_vs_pre_ranked_vector) <- bulk_rec_vs_pre_entrez$ENTREZID
rec_vs_pre_ranked_vector <- sort(rec_vs_pre_ranked_vector, decreasing = TRUE)

post_vs_pre_symbol_ranked_vector <- bulk_post_vs_pre_symbols$log2FoldChange
names(post_vs_pre_symbol_ranked_vector) <- bulk_post_vs_pre_symbols$gene
post_vs_pre_symbol_ranked_vector <- sort(post_vs_pre_symbol_ranked_vector, decreasing = TRUE)

################################################################################################################################################
######################################################      RUN REACTOME, GO, AND MITOCARTA GSEA      ##########################################
################################################################################################################################################

# We set pvalueCutoff to 1 so the GSEA objects retain the broad tested term
# universe. That is important for hierarchy-aware post-processing because
# non-significant parent terms may still need to be displayed to connect a
# significant descendant branch cleanly.
#
# We also use minGSSize = 5 and maxGSSize = 10000 because broad parent terms
# such as top-level Reactome pathways can otherwise be excluded before
# hierGSEA ever sees them. The package expects users to keep the upstream GSEA
# result permissive and then let hierGSEA handle the hierarchy-aware filtering
# and multiple-testing correction afterwards.

reactome_post_vs_pre <- ReactomePA::gsePathway(
    geneList = post_vs_pre_ranked_vector,
    organism = "human",
    maxGSSize = 10000,
    minGSSize = 5,
    pvalueCutoff = 1,
    pAdjustMethod = "none",
    verbose = FALSE
)

go_bp_post_vs_pre <- clusterProfiler::gseGO(
    geneList = post_vs_pre_ranked_vector,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    maxGSSize = 10000,
    minGSSize = 5,
    pvalueCutoff = 1,
    pAdjustMethod = "none",
    verbose = FALSE
)

if (isTRUE(run_mitocarta_example)) {
    mitocarta_post_vs_pre <- clusterProfiler::GSEA(
        geneList = post_vs_pre_symbol_ranked_vector,
        TERM2GENE = hierGSEA::mitocarta_term2gene(),
        TERM2NAME = hierGSEA::mitocarta_term2name(),
        maxGSSize = 10000,
        minGSSize = 5,
        pvalueCutoff = 1,
        pAdjustMethod = "none",
        verbose = FALSE
    )
}

################################################################################################################################################
######################################################      APPLY HIERARCHY-AWARE POST-PROCESSING      ########################################
################################################################################################################################################

# The hierarchy-aware wrapper is the core step that analysts will usually care
# about. The settings below keep the example focused on biologically readable
# top-to-mid-level branches rather than the full deepest ontology depth.

reactome_post_hier <- hierGSEA::hier_gsea(
    result = reactome_post_vs_pre,
    db = "reactome",
    directional = "both",
    level_top = 2,
    level_bottom = 5,
    alpha = 0.05
)

go_bp_post_hier <- hierGSEA::hier_gsea(
    result = go_bp_post_vs_pre,
    db = "go",
    ontology = "BP",
    directional = "both",
    level_top = 1,
    level_bottom = 6,
    alpha = 0.05
)

if (isTRUE(run_mitocarta_example)) {
    mitocarta_post_hier <- hierGSEA::hier_gsea(
        result = mitocarta_post_vs_pre,
        db = "mitocarta",
        directional = "both",
        level_top = 1,
        level_bottom = 3,
        alpha = 0.05
    )
}

################################################################################################################################################
######################################################      PLOT REACTOME HIERARCHY-AWARE RESULTS      #########################################
################################################################################################################################################

# The Reactome plot is a good first diagnostic because the pathway tree is
# relatively compact and easy to inspect. Saving immediately after the plot
# block keeps the script chronological and makes failed rendering easier to
# trace during debugging.
#
# Reactome level-1 parent term options in this bundled example are:
# "Signal Transduction"
# "Hemostasis"
# "Extracellular matrix organization"
# "Immune System"
# "Disease"
# "Cell-Cell communication"
# "Gene expression (Transcription)"
# "Developmental Biology"
# "Autophagy"
# "Vesicle-mediated transport"
# "Cellular responses to stimuli"
#
# If you want to focus the plot on a subset of major branches, you can either:
# 1. use top_n_parents = 3
# 2. use parent_terms = c("Immune System", "Signal Transduction")

reactome_post_plot <- hierGSEA::plot_hier_gsea(
    x = reactome_post_hier,
    label_col = "Description",
    size_col = "abs_NES",
    colour_col = "p_adjust_hier",
    show_left_hierarchy = TRUE,
    tree_width = 0.4,
    top_n_parents = 5
    #parent_terms = c("Immune System", "DNA Repair")
)

ggplot2::ggsave(
    reactome_post_plot,
    filename = file.path(example_output_dir, "single_fiber_post_vs_pre_reactome_hierarchy_top_5.pdf"),
    width = 20,
    height = 45,
    units = "cm"
)

################################################################################################################################################
######################################################      PLOT GO BP HIERARCHY-AWARE RESULTS      ############################################
################################################################################################################################################

# GO biological process has a much denser DAG structure than Reactome, so the
# example keeps only a few top-level branches by default. That makes the saved
# figure much easier to inspect during a first test run and avoids producing a
# crowded plot that needs immediate manual pruning.
#
# For GO, hierGSEA intentionally removes the artificial ontology container
# nodes such as "all" and "biological_process". That means the first visible
# GO level now corresponds to real biological-process branches rather than to
# ontology headers.

go_bp_post_plot <- hierGSEA::plot_hier_gsea(
    x = go_bp_post_hier,
    label_col = "Description",
    size_col = "abs_NES",
    colour_col = "p_adjust_hier",
    show_left_hierarchy = TRUE,
    tree_width = 0.32, 
    top_n_parents = 3, 
    significance_cutoff = 0.05
)

ggplot2::ggsave(
    go_bp_post_plot,
    filename = file.path(example_output_dir, "single_fiber_post_vs_pre_go_bp_hierarchy.pdf"),
    width = 36,
    height = 34,
    units = "cm"
)

################################################################################################################################################
######################################################      PLOT MITOCARTA HIERARCHY-AWARE RESULTS      ########################################
################################################################################################################################################

# MitoCarta is useful here as a proof that hierGSEA can also handle a custom
# hierarchical gene-set universe rather than only databases that ship directly
# inside clusterProfiler. The bundled backend uses the Broad MitoPathways3.0
# hierarchy, while the enrichment itself still comes from clusterProfiler::GSEA.
# We again keep the default branch count low so the example output remains
# publication-style readable on a normal page.

if (isTRUE(run_mitocarta_example)) {
    mitocarta_post_plot <- hierGSEA::plot_hier_gsea(
        x = mitocarta_post_hier,
        label_col = "Description",
        size_col = "abs_NES",
        colour_col = "p_adjust_hier",
        show_left_hierarchy = TRUE,
        tree_width = 0.34,
        top_n_parents = 3,
        significance_cutoff = 0.05
    )

    ggplot2::ggsave(
        mitocarta_post_plot,
        filename = file.path(example_output_dir, "single_fiber_post_vs_pre_mitocarta_hierarchy.pdf"),
        width = 36,
        height = 28,
        units = "cm"
    )
}

################################################################################################################################################
######################################################      SAVE TABLES FOR QUICK INSPECTION      ################################################
################################################################################################################################################

# Writing the processed tables is useful during debugging because it lets you
# inspect the retained branch order and the hierarchy-aware adjusted p-values
# without needing to interrogate objects interactively.

readr::write_csv(
    reactome_post_hier$results_tbl,
    file.path(example_output_dir, "single_fiber_post_vs_pre_reactome_hierarchy_results.csv")
)

readr::write_csv(
    go_bp_post_hier$results_tbl,
    file.path(example_output_dir, "single_fiber_post_vs_pre_go_bp_hierarchy_results.csv")
)

if (isTRUE(run_mitocarta_example)) {
    readr::write_csv(
        mitocarta_post_hier$results_tbl,
        file.path(example_output_dir, "single_fiber_post_vs_pre_mitocarta_hierarchy_results.csv")
    )
}

message("hierGSEA example completed. Output files were written to: ", example_output_dir)
