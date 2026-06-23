# This convenience script rebuilds all internal package data from the local
# source repositories and installed annotation packages. Running this script is
# the fastest way to refresh the hierarchy snapshots after code changes.

source(file.path("data-raw", "build_reactome_hierarchy.R"))
source(file.path("data-raw", "build_go_hierarchies.R"))
source(file.path("data-raw", "build_mitocarta_hierarchy.R"))
source(file.path("data-raw", "build_example_data.R"))

find_mitocarta_html_source <- function() {
    candidate_paths <- c(
        Sys.getenv("HIERGSEA_MITOCARTA_HTML", unset = ""),
        file.path("data-raw", "human.mitopathways3.0.html"),
        "/tmp/human.mitopathways3.0.html",
        file.path(tempdir(), "human.mitopathways3.0.html")
    )

    candidate_paths <- candidate_paths[nzchar(candidate_paths)]
    existing_path <- candidate_paths[file.exists(candidate_paths)][1]

    if (length(existing_path) == 0 || is.na(existing_path)) {
        stop(
            paste(
                "No local MitoPathways3.0 HTML source was found.",
                "Set HIERGSEA_MITOCARTA_HTML to a downloaded copy of",
                "'human.mitopathways3.0.html', place that file in data-raw/,",
                "or run hierGSEA::update_backend_data(force = TRUE, databases = 'mitocarta')",
                "to refresh the backend from the official Broad source first."
            ),
            call. = FALSE
        )
    }

    normalizePath(existing_path, winslash = "/")
}

reactome_hierarchy_data <- build_reactome_hierarchy_data()
go_bp_hierarchy_data <- build_go_hierarchy_data("BP")
go_mf_hierarchy_data <- build_go_hierarchy_data("MF")
go_cc_hierarchy_data <- build_go_hierarchy_data("CC")
mitocarta_build <- build_mitocarta_hierarchy_data(
    source_html_path = find_mitocarta_html_source()
)
mitocarta_hierarchy_data <- mitocarta_build$hierarchy_data
mitocarta_term2gene_tbl <- mitocarta_build$term2gene_tbl
mitocarta_term2name_tbl <- mitocarta_build$term2name_tbl
single_fiber_bulk_rankings <- build_single_fiber_bulk_rankings()

dir.create("R", showWarnings = FALSE, recursive = TRUE)
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("inst", "extdata"), showWarnings = FALSE, recursive = TRUE)

save(
    reactome_hierarchy_data,
    go_bp_hierarchy_data,
    go_mf_hierarchy_data,
    go_cc_hierarchy_data,
    mitocarta_hierarchy_data,
    mitocarta_term2gene_tbl,
    mitocarta_term2name_tbl,
    file = file.path("R", "sysdata.rda"),
    compress = "xz"
)

save(
    single_fiber_bulk_rankings,
    file = file.path("data", "single_fiber_bulk_rankings.rda"),
    compress = "xz"
)

readr::write_csv(
    single_fiber_bulk_rankings$post_vs_pre,
    file = file.path("inst", "extdata", "single_fiber_bulk_post_vs_pre.csv")
)

readr::write_csv(
    single_fiber_bulk_rankings$rec_vs_pre,
    file = file.path("inst", "extdata", "single_fiber_bulk_rec_vs_pre.csv")
)
