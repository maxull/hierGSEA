suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(arrow)
    library(clusterProfiler)
    library(dplyr)
    library(org.Hs.eg.db)
    library(readr)
    library(ReactomePA)
    library(stringr)
    library(tibble)
})

################################################################################################################################################
######################################################      LOAD hierGSEA LOCALLY OR FROM AN INSTALL      ######################################
################################################################################################################################################

# This script is a real-world stress test for the multi-result plotting mode.
# It uses an external MetaMex checkout, builds matched ranked vectors for
# several exercise settings, runs fresh Reactome GSEA analyses, and then passes
# named lists of hier_gsea_result objects into plot_hier_gsea() so the shared
# tree logic can be inspected visually.

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

    source_guess <- file.path(getwd(), "inst", "scripts", "run_metamex_multi_result_example.R")

    if (file.exists(source_guess)) {
        return(normalizePath(source_guess, winslash = "/"))
    }

    stop("Could not determine the script path for the MetaMex example workflow.", call. = FALSE)
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
######################################################      DEFINE INPUT AND OUTPUT PATHS      ##################################################
################################################################################################################################################

# This workflow depends on a local MetaMex checkout because the MetaMex data
# are not bundled inside hierGSEA. You can override the default path with the
# HIERGSEA_METAMEX_ROOT environment variable if your local folder lives
# elsewhere.

script_path <- get_script_path()
package_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")

metamex_root <- Sys.getenv(
    "HIERGSEA_METAMEX_ROOT",
    unset = "/Users/maxullrich/Documents/GitHub/Side-Projects/MetaMex"
)

metamex_root <- normalizePath(metamex_root, winslash = "/", mustWork = FALSE)

if (!dir.exists(metamex_root)) {
    stop(
        paste0(
            "The MetaMex repository was not found at: ",
            metamex_root,
            ". Set HIERGSEA_METAMEX_ROOT to your local MetaMex checkout."
        ),
        call. = FALSE
    )
}

metamex_data_dir <- file.path(metamex_root, "data")

if (!dir.exists(metamex_data_dir)) {
    stop("The MetaMex data/ directory could not be found.", call. = FALSE)
}

example_output_dir <- file.path(getwd(), "hierGSEA_example_output")
dir.create(example_output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################################################################################
######################################################      DEFINE SELECTION SETTINGS      #####################################################
################################################################################################################################################

# The acute comparison uses an exact shared cohort signature across aerobic,
# resistance, and HIIT so the three ranked vectors represent genuinely matched
# MetaMex subgroup definitions.
#
# The chronic comparison is handled as one four-group figure. After checking
# the available MetaMex signatures, there are no exact young-or-middle-aged
# healthy lean-versus-overweight overlaps where every other metadata field is
# identical within each training modality. The script therefore uses the
# closest practical pooled comparison: healthy YNG + MDL contrasts only,
# separated by modality and bodyweight class.
#
# We also restrict the chronic biopsies to explicit post-exercise timepoints of
# at least 24 hours. This avoids mixing likely near-bout samples into the
# chronic adaptation view. HNA is excluded here because it does not guarantee a
# sufficiently delayed post-exercise biopsy.

acute_selection <- list(
    age = "YNG",
    training_status = "ACT",
    weight_class = "LEA",
    health_status = "HLY",
    timepoint = "H03"
)

chronic_training_healthy_selection <- list(
    exercise_type = "MIX",
    age = c("YNG", "MDL"),
    health_status = "HLY",
    timepoint = c("H24", "H48", "H72", "H96")
)

################################################################################################################################################
######################################################      LOAD METAMEX STAT TABLES      ######################################################
################################################################################################################################################

# We only load the tables needed for the requested comparisons. Each table
# contains one gene column and many contrast columns, where the contrast name
# encodes the subgroup metadata that we will parse in the next step.

stats_human_acute_aerobic <- arrow::read_feather(
    file.path(metamex_data_dir, "stats_human_acute_aerobic.feather")
)

stats_human_acute_hit <- arrow::read_feather(
    file.path(metamex_data_dir, "stats_human_acute_hit.feather")
)

stats_human_acute_resistance <- arrow::read_feather(
    file.path(metamex_data_dir, "stats_human_acute_resistance.feather")
)

stats_human_training_aerobic <- arrow::read_feather(
    file.path(metamex_data_dir, "stats_human_training_aerobic.feather")
)

stats_human_training_resistance <- arrow::read_feather(
    file.path(metamex_data_dir, "stats_human_training_resistance.feather")
)

################################################################################################################################################
######################################################      PARSE METAMEX CONTRAST METADATA      ################################################
################################################################################################################################################

# MetaMex stores subgroup information directly in the contrast column names.
# Parsing those names up front lets us select matched exercise comparisons
# explicitly, which is much safer than hard-coding raw column names later.

parse_metamex_metadata <- function(stats_tbl) {
    contrast_columns <- setdiff(colnames(stats_tbl), "SYMBOL")
    parsed <- stringr::str_split_fixed(contrast_columns, pattern = "_", n = 12)

    tibble::tibble(
        column = contrast_columns,
        stat = parsed[, 1],
        contrast_id = sub("^[^_]+_", "", contrast_columns),
        gse = dplyr::na_if(parsed[, 2], ""),
        condition = dplyr::na_if(parsed[, 3], ""),
        exercise_type = dplyr::na_if(parsed[, 4], ""),
        muscle = dplyr::na_if(parsed[, 5], ""),
        sex = dplyr::na_if(parsed[, 6], ""),
        age = dplyr::na_if(parsed[, 7], ""),
        training_status = dplyr::na_if(parsed[, 8], ""),
        weight_class = dplyr::na_if(parsed[, 9], ""),
        health_status = dplyr::na_if(parsed[, 10], ""),
        timepoint = dplyr::na_if(parsed[, 11], ""),
        intervention_length = dplyr::na_if(parsed[, 12], "")
    )
}

acute_aerobic_meta <- parse_metamex_metadata(stats_human_acute_aerobic)
acute_hit_meta <- parse_metamex_metadata(stats_human_acute_hit)
acute_resistance_meta <- parse_metamex_metadata(stats_human_acute_resistance)
training_aerobic_meta <- parse_metamex_metadata(stats_human_training_aerobic)
training_resistance_meta <- parse_metamex_metadata(stats_human_training_resistance)

################################################################################################################################################
######################################################      HELPER TO SELECT MATCHED CONTRASTS      ############################################
################################################################################################################################################

# We keep the filters explicit and human-readable so that analysts can quickly
# change the cohort definition when they want to probe another MetaMex subgroup.

filter_metamex_metadata <- function(metadata_tbl, selection_filters) {
    selected_tbl <- metadata_tbl

    for (filter_name in names(selection_filters)) {
        filter_value <- selection_filters[[filter_name]]

        selected_tbl <- selected_tbl %>%
            dplyr::filter(.data[[filter_name]] %in% filter_value)
    }

    selected_tbl
}

################################################################################################################################################
######################################################      HELPER TO AGGREGATE STUDY-LEVEL CONTRASTS      #####################################
################################################################################################################################################

# For each chosen MetaMex subgroup we average the available logFC contrasts
# gene-wise, weighting by study size when the size columns are available.
# This gives us one stable ranked vector per exercise setting while still
# respecting that the subgroup may be represented by multiple studies.

aggregate_metamex_group <- function(stats_tbl, metadata_tbl, selection_filters, group_label) {
    selected_logfc <- filter_metamex_metadata(metadata_tbl, selection_filters) %>%
        dplyr::filter(.data$stat == "logFC")

    if (nrow(selected_logfc) == 0) {
        stop(
            paste0(
                "No MetaMex contrasts matched the requested filters for group: ",
                group_label
            ),
            call. = FALSE
        )
    }

    selected_contrast_order <- selected_logfc$contrast_id

    selected_logfc <- selected_logfc %>%
        dplyr::mutate(contrast_id = factor(.data$contrast_id, levels = selected_contrast_order)) %>%
        dplyr::arrange(.data$contrast_id)

    selected_size <- metadata_tbl %>%
        dplyr::filter(.data$stat == "size", .data$contrast_id %in% selected_contrast_order) %>%
        dplyr::mutate(contrast_id = factor(.data$contrast_id, levels = selected_contrast_order)) %>%
        dplyr::arrange(.data$contrast_id)

    logfc_mat <- as.matrix(stats_tbl[, selected_logfc$column, drop = FALSE])
    mode(logfc_mat) <- "numeric"

    if (nrow(selected_size) == nrow(selected_logfc)) {
        size_mat <- as.matrix(stats_tbl[, selected_size$column, drop = FALSE])
        mode(size_mat) <- "numeric"
        weighted_numerator <- rowSums(logfc_mat * size_mat, na.rm = TRUE)
        weighted_denominator <- rowSums(size_mat, na.rm = TRUE)
        aggregated_logfc <- weighted_numerator / ifelse(weighted_denominator == 0, NA_real_, weighted_denominator)
        mean_size <- rowMeans(size_mat, na.rm = TRUE)
    } else {
        aggregated_logfc <- rowMeans(logfc_mat, na.rm = TRUE)
        mean_size <- rep(NA_real_, nrow(stats_tbl))
    }

    tibble::tibble(
        SYMBOL = stats_tbl$SYMBOL,
        logFC = aggregated_logfc,
        mean_size = mean_size,
        n_contrasts = rowSums(!is.na(logfc_mat)),
        group_label = group_label
    ) %>%
        dplyr::filter(!is.na(.data$SYMBOL), nzchar(.data$SYMBOL), is.finite(.data$logFC)) %>%
        dplyr::distinct(.data$SYMBOL, .keep_all = TRUE) %>%
        dplyr::arrange(dplyr::desc(abs(.data$logFC)), .data$SYMBOL)
}

################################################################################################################################################
######################################################      HELPER TO CONVERT SYMBOLS TO ENTREZ IDS      #######################################
################################################################################################################################################

# Reactome GSEA is most stable with Entrez IDs. We keep this conversion local
# to the script so the package API remains focused on hierarchy-aware
# post-processing rather than on MetaMex-specific identifier preparation.

convert_symbols_to_entrez <- function(ranking_df) {
    id_lookup <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = ranking_df$SYMBOL,
        keytype = "SYMBOL",
        columns = c("SYMBOL", "ENTREZID")
    ) %>%
        tibble::as_tibble() %>%
        dplyr::filter(!is.na(.data$ENTREZID))

    ranking_df %>%
        dplyr::inner_join(id_lookup, by = "SYMBOL") %>%
        dplyr::arrange(dplyr::desc(.data$logFC), .data$ENTREZID) %>%
        dplyr::distinct(.data$ENTREZID, .keep_all = TRUE) %>%
        dplyr::arrange(dplyr::desc(.data$logFC), .data$ENTREZID)
}

build_ranked_vector <- function(ranking_df) {
    entrez_tbl <- convert_symbols_to_entrez(ranking_df)
    ranked_vector <- entrez_tbl$logFC
    names(ranked_vector) <- entrez_tbl$ENTREZID
    sort(ranked_vector, decreasing = TRUE)
}

################################################################################################################################################
######################################################      BUILD MATCHED ACUTE META ANALYSIS GROUPS      ######################################
################################################################################################################################################

# We test the three-way multi-result plot on an exact early acute signature
# that exists in aerobic, resistance, and HIIT: young, active, lean, healthy,
# at H03. The resistance arm intentionally pools the available matched
# resistance modes at that timepoint because MetaMex stores those cohorts as
# ECC and CON rather than as one generic MIX label.

acute_aerobic_rank_tbl <- aggregate_metamex_group(
    stats_tbl = stats_human_acute_aerobic,
    metadata_tbl = acute_aerobic_meta,
    selection_filters = acute_selection,
    group_label = "Acute Aerobic"
)

acute_resistance_rank_tbl <- aggregate_metamex_group(
    stats_tbl = stats_human_acute_resistance,
    metadata_tbl = acute_resistance_meta,
    selection_filters = acute_selection,
    group_label = "Acute Resistance"
)

acute_hiit_rank_tbl <- aggregate_metamex_group(
    stats_tbl = stats_human_acute_hit,
    metadata_tbl = acute_hit_meta,
    selection_filters = acute_selection,
    group_label = "Acute HIIT"
)

################################################################################################################################################
######################################################      BUILD MATCHED CHRONIC META ANALYSIS GROUPS      ####################################
################################################################################################################################################

# The chronic comparison is a pooled healthy lean-versus-overweight view for
# aerobic and resistance training, restricted to young and middle-aged cohorts
# with biopsies taken at least 24 hours after the last exercise bout.

chronic_healthy_lean_aerobic_rank_tbl <- aggregate_metamex_group(
    stats_tbl = stats_human_training_aerobic,
    metadata_tbl = training_aerobic_meta,
    selection_filters = c(chronic_training_healthy_selection, list(weight_class = "LEA")),
    group_label = "Training Aerobic Healthy"
)

chronic_healthy_overweight_aerobic_rank_tbl <- aggregate_metamex_group(
    stats_tbl = stats_human_training_aerobic,
    metadata_tbl = training_aerobic_meta,
    selection_filters = c(chronic_training_healthy_selection, list(weight_class = "OWE")),
    group_label = "Training Aerobic Overweight"
)

chronic_healthy_lean_resistance_rank_tbl <- aggregate_metamex_group(
    stats_tbl = stats_human_training_resistance,
    metadata_tbl = training_resistance_meta,
    selection_filters = c(chronic_training_healthy_selection, list(weight_class = "LEA")),
    group_label = "Training Resistance Healthy"
)

chronic_healthy_overweight_resistance_rank_tbl <- aggregate_metamex_group(
    stats_tbl = stats_human_training_resistance,
    metadata_tbl = training_resistance_meta,
    selection_filters = c(chronic_training_healthy_selection, list(weight_class = "OWE")),
    group_label = "Training Resistance Overweight"
)

################################################################################################################################################
######################################################      SAVE MATCHED CONTRAST TABLES FOR DEBUGGING      ####################################
################################################################################################################################################

# These tables make it easy to confirm exactly which MetaMex contrasts fed each
# ranked vector, which is especially useful when you want to adjust the cohort
# filter logic later.

readr::write_csv(
    filter_metamex_metadata(acute_aerobic_meta, acute_selection),
    file.path(example_output_dir, "metamex_acute_aerobic_selected_contrasts.csv")
)

readr::write_csv(
    filter_metamex_metadata(acute_resistance_meta, acute_selection),
    file.path(example_output_dir, "metamex_acute_resistance_selected_contrasts.csv")
)

readr::write_csv(
    filter_metamex_metadata(acute_hit_meta, acute_selection),
    file.path(example_output_dir, "metamex_acute_hiit_selected_contrasts.csv")
)

readr::write_csv(
    filter_metamex_metadata(
        training_aerobic_meta,
        c(chronic_training_healthy_selection, list(weight_class = "LEA"))
    ),
    file.path(example_output_dir, "metamex_training_aerobic_healthy_selected_contrasts.csv")
)

readr::write_csv(
    filter_metamex_metadata(
        training_aerobic_meta,
        c(chronic_training_healthy_selection, list(weight_class = "OWE"))
    ),
    file.path(example_output_dir, "metamex_training_aerobic_overweight_selected_contrasts.csv")
)

readr::write_csv(
    filter_metamex_metadata(
        training_resistance_meta,
        c(chronic_training_healthy_selection, list(weight_class = "LEA"))
    ),
    file.path(example_output_dir, "metamex_training_resistance_healthy_selected_contrasts.csv")
)

readr::write_csv(
    filter_metamex_metadata(
        training_resistance_meta,
        c(chronic_training_healthy_selection, list(weight_class = "OWE"))
    ),
    file.path(example_output_dir, "metamex_training_resistance_overweight_selected_contrasts.csv")
)

################################################################################################################################################
######################################################      BUILD RANKED VECTORS      ##########################################################
################################################################################################################################################

# The ranked vectors are built directly from the aggregated MetaMex logFC
# values, which keeps the workflow close to the single-fiber example while now
# operating on study-combined contrasts instead of on one single DE table.

acute_aerobic_ranked_vector <- build_ranked_vector(acute_aerobic_rank_tbl)
acute_resistance_ranked_vector <- build_ranked_vector(acute_resistance_rank_tbl)
acute_hiit_ranked_vector <- build_ranked_vector(acute_hiit_rank_tbl)

chronic_healthy_lean_aerobic_ranked_vector <- build_ranked_vector(chronic_healthy_lean_aerobic_rank_tbl)
chronic_healthy_overweight_aerobic_ranked_vector <- build_ranked_vector(chronic_healthy_overweight_aerobic_rank_tbl)
chronic_healthy_lean_resistance_ranked_vector <- build_ranked_vector(chronic_healthy_lean_resistance_rank_tbl)
chronic_healthy_overweight_resistance_ranked_vector <- build_ranked_vector(chronic_healthy_overweight_resistance_rank_tbl)

################################################################################################################################################
######################################################      RUN FRESH REACTOME GSEA      #######################################################
################################################################################################################################################

# We keep the upstream Reactome run permissive for the same reason as in the
# package examples: hierGSEA needs the broad tested term universe so it can do
# the hierarchy-aware retention and correction afterwards.

run_reactome_gsea <- function(ranked_vector) {
    ReactomePA::gsePathway(
        geneList = ranked_vector,
        organism = "human",
        maxGSSize = 10000,
        minGSSize = 5,
        pvalueCutoff = 1,
        pAdjustMethod = "none",
        verbose = FALSE
    )
}

acute_aerobic_reactome <- run_reactome_gsea(acute_aerobic_ranked_vector)
acute_resistance_reactome <- run_reactome_gsea(acute_resistance_ranked_vector)
acute_hiit_reactome <- run_reactome_gsea(acute_hiit_ranked_vector)

chronic_healthy_lean_aerobic_reactome <- run_reactome_gsea(chronic_healthy_lean_aerobic_ranked_vector)
chronic_healthy_overweight_aerobic_reactome <- run_reactome_gsea(chronic_healthy_overweight_aerobic_ranked_vector)
chronic_healthy_lean_resistance_reactome <- run_reactome_gsea(chronic_healthy_lean_resistance_ranked_vector)
chronic_healthy_overweight_resistance_reactome <- run_reactome_gsea(chronic_healthy_overweight_resistance_ranked_vector)

################################################################################################################################################
######################################################      APPLY HIERARCHY-AWARE POST-PROCESSING      ########################################
################################################################################################################################################

# All results within one multi-result plot must share the same hierarchy
# settings. We therefore process the matched MetaMex groups with the same
# directional and level-window choices before combining them in the plot call.

acute_aerobic_hier <- hierGSEA::hier_gsea(
    result = acute_aerobic_reactome,
    db = "reactome",
    directional = "both",
    level_top = 1,
    level_bottom = 4,
    alpha = 0.05
)

acute_resistance_hier <- hierGSEA::hier_gsea(
    result = acute_resistance_reactome,
    db = "reactome",
    directional = "both",
    level_top = 1,
    level_bottom = 4,
    alpha = 0.05
)

acute_hiit_hier <- hierGSEA::hier_gsea(
    result = acute_hiit_reactome,
    db = "reactome",
    directional = "both",
    level_top = 1,
    level_bottom = 4,
    alpha = 0.05
)

chronic_healthy_lean_aerobic_hier <- hierGSEA::hier_gsea(
    result = chronic_healthy_lean_aerobic_reactome,
    db = "reactome",
    directional = "both",
    level_top = 2,
    level_bottom = 10,
    alpha = 0.1
)

chronic_healthy_overweight_aerobic_hier <- hierGSEA::hier_gsea(
    result = chronic_healthy_overweight_aerobic_reactome,
    db = "reactome",
    directional = "both",
    level_top = 2,
    level_bottom = 10,
    alpha = 0.1
)

chronic_healthy_lean_resistance_hier <- hierGSEA::hier_gsea(
    result = chronic_healthy_lean_resistance_reactome,
    db = "reactome",
    directional = "both",
    level_top = 2,
    level_bottom = 10,
    alpha = 0.1
)

chronic_healthy_overweight_resistance_hier <- hierGSEA::hier_gsea(
    result = chronic_healthy_overweight_resistance_reactome,
    db = "reactome",
    directional = "both",
    level_top = 2,
    level_bottom = 10,
    alpha = 0.1
)

################################################################################################################################################
######################################################      PLOT ACUTE MULTI-RESULT SHARED TREE      ###########################################
################################################################################################################################################

# This is the main test of the new functionality: three independent
# hierarchy-aware GSEA results are supplied as a named list, and the plot uses
# one shared Reactome tree with the comparison blocks arranged left to right in
# the order of that list.

acute_reactome_multiplot <- hierGSEA::plot_hier_gsea(
    x = list(
        "Acute Aerobic" = acute_aerobic_hier,
        "Acute Resistance" = acute_resistance_hier,
        "Acute HIIT" = acute_hiit_hier
    ),
    label_col = "Description",
    size_col = "abs_NES",
    colour_col = "p_adjust_hier",
    show_left_hierarchy = TRUE,
    tree_width = 0.42,
    top_n_parents = 2,
    significance_cutoff = 0.05
)

ggplot2::ggsave(
    acute_reactome_multiplot,
    filename = file.path(example_output_dir, "metamex_acute_aerobic_vs_resistance_vs_hiit_reactome_hierarchy.pdf"),
    width = 48,
    height = 30,
    units = "cm"
)

################################################################################################################################################
######################################################      PLOT CHRONIC HEALTHY VS OVERWEIGHT SHARED TREE      ################################
################################################################################################################################################

# This figure compares aerobic and resistance training while keeping the
# chronic selection restricted to healthy YNG + MDL cohorts and splitting the
# pooled contrasts by bodyweight class. Because exact lean-versus-overweight
# matched signatures were not available in MetaMex for this age range, this
# broad pooled comparison is the closest stable test case for the new
# multi-result shared-tree layout.

chronic_weight_reactome_multiplot <- hierGSEA::plot_hier_gsea(
    x = list(
        "Training Aerobic \nHealthy" = chronic_healthy_lean_aerobic_hier,
        "Training Aerobic \nOverweight" = chronic_healthy_overweight_aerobic_hier,
        "Training Resistance \nHealthy" = chronic_healthy_lean_resistance_hier,
        "Training Resistance \nOverweight" = chronic_healthy_overweight_resistance_hier
    ),
    label_col = "Description",
    size_col = "abs_NES",
    colour_col = "p_adjust_hier",
    show_left_hierarchy = TRUE,
    tree_width = 0.42,
    parent_terms = c("Cytokine Signaling in Immune System"),
    #top_n_parents = 2,
    significance_cutoff = 0.05
)

ggplot2::ggsave(
    chronic_weight_reactome_multiplot,
    filename = file.path(example_output_dir, "metamex_chronic_healthy_vs_overweight_aerobic_and_resistance_reactome_hierarchy.pdf"),
    width = 56,
    height = 40,
    units = "cm"
)

################################################################################################################################################
######################################################      SAVE HIERARCHY TABLES FOR DEBUGGING      ############################################
################################################################################################################################################

# These processed tables are often the fastest way to debug branch ordering,
# retained ancestors, and comparison-specific significance patterns before
# looking at the figure itself.

readr::write_csv(
    acute_aerobic_hier$results_tbl,
    file.path(example_output_dir, "metamex_acute_aerobic_reactome_hierarchy_results.csv")
)

readr::write_csv(
    acute_resistance_hier$results_tbl,
    file.path(example_output_dir, "metamex_acute_resistance_reactome_hierarchy_results.csv")
)

readr::write_csv(
    acute_hiit_hier$results_tbl,
    file.path(example_output_dir, "metamex_acute_hiit_reactome_hierarchy_results.csv")
)

readr::write_csv(
    chronic_healthy_lean_aerobic_hier$results_tbl,
    file.path(example_output_dir, "metamex_chronic_healthy_aerobic_reactome_hierarchy_results.csv")
)

readr::write_csv(
    chronic_healthy_overweight_aerobic_hier$results_tbl,
    file.path(example_output_dir, "metamex_chronic_overweight_aerobic_reactome_hierarchy_results.csv")
)

readr::write_csv(
    chronic_healthy_lean_resistance_hier$results_tbl,
    file.path(example_output_dir, "metamex_chronic_healthy_resistance_reactome_hierarchy_results.csv")
)

readr::write_csv(
    chronic_healthy_overweight_resistance_hier$results_tbl,
    file.path(example_output_dir, "metamex_chronic_overweight_resistance_reactome_hierarchy_results.csv")
)

message("MetaMex multi-result hierGSEA example completed. Output files were written to: ", example_output_dir)
