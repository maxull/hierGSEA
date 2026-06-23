#' Update stored Reactome and GO hierarchy backends
#'
#' `update_backend_data()` checks the official Reactome and Gene Ontology
#' releases online, compares them with the versions stored in the package
#' backend metadata, and only rebuilds `R/sysdata.rda` when a newer release is
#' available or when `force = TRUE`.
#'
#' Because the package ships prebuilt internal data, this function is intended
#' to be run from the source repository rather than from an installed package
#' library. After rebuilding the backend snapshot, reinstall or reload the
#' package so the new internal data are used in future sessions.
#'
#' @param package_root Path to the package source directory. Defaults to `"."`.
#'   The directory must contain `DESCRIPTION`, `R/`, and `data-raw/`.
#' @param databases Character vector containing one or more of `"reactome"`,
#'   `"go"`, and `"mitocarta"`.
#' @param force Logical. If `TRUE`, rebuild even when the stored versions match
#'   the online versions.
#' @param verbose Logical. If `TRUE`, print version summaries and status
#'   messages.
#'
#' @return Invisibly returns a named list describing the local versions, the
#'   online versions, whether each backend was rebuilt, and the output
#'   `sysdata.rda` path.
#' @export
update_backend_data <- function(
    package_root = ".",
    databases = c("reactome", "go"),
    force = FALSE,
    verbose = TRUE
) {
    package_root <- .resolve_hiergsea_package_root(package_root = package_root)
    databases <- unique(tolower(databases))

    if (!all(databases %in% c("reactome", "go", "mitocarta"))) {
        stop("databases must contain only 'reactome', 'go', and/or 'mitocarta'.", call. = FALSE)
    }

    backend_snapshot <- .read_local_backend_versions(package_root = package_root)
    online_versions <- .fetch_online_backend_versions(databases = databases)

    rebuild_reactome <- "reactome" %in% databases && (
        isTRUE(force) ||
            !identical(backend_snapshot$reactome$release, online_versions$reactome$release)
    )

    rebuild_go <- "go" %in% databases && (
        isTRUE(force) ||
            !identical(backend_snapshot$go$release, online_versions$go$release)
    )

    rebuild_mitocarta <- "mitocarta" %in% databases && (
        isTRUE(force) ||
            !identical(backend_snapshot$mitocarta$release, online_versions$mitocarta$release)
    )

    if (isTRUE(verbose)) {
        .print_backend_version_summary(
            local_versions = backend_snapshot,
            online_versions = online_versions,
            rebuild_reactome = rebuild_reactome,
            rebuild_go = rebuild_go,
            rebuild_mitocarta = rebuild_mitocarta
        )
    }

    if (!rebuild_reactome && !rebuild_go && !rebuild_mitocarta) {
        if (isTRUE(verbose)) {
            message("hierGSEA backend data are already up to date. No rebuild was needed.")
        }

        return(invisible(list(
            local_versions = backend_snapshot,
            online_versions = online_versions,
            rebuilt = c(reactome = FALSE, go = FALSE, mitocarta = FALSE),
            sysdata_path = file.path(package_root, "R", "sysdata.rda")
        )))
    }

    total_steps <- 3L +
        as.integer(rebuild_reactome) +
        as.integer(rebuild_go) +
        as.integer(rebuild_mitocarta)
    progress_bar <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
    on.exit(close(progress_bar), add = TRUE)

    progress_step <- 0
    tmp_dir <- tempfile(pattern = "hiergsea_backend_update_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    updated_objects <- .load_existing_backend_objects(package_root = package_root)
    progress_step <- progress_step + 1
    utils::setTxtProgressBar(progress_bar, progress_step)

    if (rebuild_reactome) {
        reactome_downloads <- .download_reactome_backend_files(
            out_dir = tmp_dir,
            reactome_release = online_versions$reactome$release
        )

        progress_step <- progress_step + 1
        utils::setTxtProgressBar(progress_bar, progress_step)

        updated_objects$reactome_hierarchy_data <- .build_reactome_backend_from_downloads(
            package_root = package_root,
            relation_path = reactome_downloads$relation_path,
            pathway_path = reactome_downloads$pathway_path,
            reactome_release = online_versions$reactome$release
        )
    }

    progress_step <- progress_step + 1
    utils::setTxtProgressBar(progress_bar, progress_step)

    if (rebuild_go) {
        go_download <- .download_go_obo(
            out_dir = tmp_dir,
            go_release = online_versions$go$release
        )

        progress_step <- progress_step + 1
        utils::setTxtProgressBar(progress_bar, progress_step)

        updated_objects$go_bp_hierarchy_data <- .build_go_backend_from_obo(
            package_root = package_root,
            obo_path = go_download$obo_path,
            ontology = "BP",
            go_release = online_versions$go$release
        )

        updated_objects$go_mf_hierarchy_data <- .build_go_backend_from_obo(
            package_root = package_root,
            obo_path = go_download$obo_path,
            ontology = "MF",
            go_release = online_versions$go$release
        )

        updated_objects$go_cc_hierarchy_data <- .build_go_backend_from_obo(
            package_root = package_root,
            obo_path = go_download$obo_path,
            ontology = "CC",
            go_release = online_versions$go$release
        )
    }

    progress_step <- progress_step + 1
    utils::setTxtProgressBar(progress_bar, progress_step)

    if (rebuild_mitocarta) {
        mitocarta_download <- .download_mitocarta_html(
            out_dir = tmp_dir,
            mitocarta_release = online_versions$mitocarta$release
        )

        progress_step <- progress_step + 1
        utils::setTxtProgressBar(progress_bar, progress_step)

        mitocarta_build <- .build_mitocarta_backend_from_html(
            package_root = package_root,
            html_path = mitocarta_download$html_path,
            mitocarta_release = online_versions$mitocarta$release
        )

        updated_objects$mitocarta_hierarchy_data <- mitocarta_build$hierarchy_data
        updated_objects$mitocarta_term2gene_tbl <- mitocarta_build$term2gene_tbl
        updated_objects$mitocarta_term2name_tbl <- mitocarta_build$term2name_tbl
    }

    progress_step <- progress_step + 1
    utils::setTxtProgressBar(progress_bar, progress_step)

    save(
        updated_objects$reactome_hierarchy_data,
        updated_objects$go_bp_hierarchy_data,
        updated_objects$go_mf_hierarchy_data,
        updated_objects$go_cc_hierarchy_data,
        updated_objects$mitocarta_hierarchy_data,
        updated_objects$mitocarta_term2gene_tbl,
        updated_objects$mitocarta_term2name_tbl,
        file = file.path(package_root, "R", "sysdata.rda"),
        compress = "xz"
    )

    progress_step <- progress_step + 1
    utils::setTxtProgressBar(progress_bar, progress_step)

    if (isTRUE(verbose)) {
        message("hierGSEA backend update completed: ", file.path(package_root, "R", "sysdata.rda"))
    }

    invisible(list(
        local_versions = backend_snapshot,
        online_versions = online_versions,
        rebuilt = c(
            reactome = rebuild_reactome,
            go = rebuild_go,
            mitocarta = rebuild_mitocarta
        ),
        sysdata_path = file.path(package_root, "R", "sysdata.rda")
    ))
}

.resolve_hiergsea_package_root <- function(package_root = ".") {
    package_root <- normalizePath(package_root, winslash = "/", mustWork = FALSE)

    required_paths <- c(
        file.path(package_root, "DESCRIPTION"),
        file.path(package_root, "R"),
        file.path(package_root, "data-raw")
    )

    if (!all(file.exists(required_paths))) {
        stop(
            paste0(
                "package_root must point to the hierGSEA source repository. Missing one or more of: ",
                paste(required_paths, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    package_root
}

.read_local_backend_versions <- function(package_root) {
    sysdata_path <- file.path(package_root, "R", "sysdata.rda")

    if (!file.exists(sysdata_path)) {
        return(list(
            reactome = list(release = NA_character_),
            go = list(release = NA_character_),
            mitocarta = list(release = NA_character_)
        ))
    }

    backend_env <- new.env(parent = emptyenv())
    load(sysdata_path, envir = backend_env)

    list(
        reactome = list(
            release = backend_env$reactome_hierarchy_data$metadata$reactome_release %||%
                backend_env$reactome_hierarchy_data$metadata$relation_file_mtime %||%
                NA_character_
        ),
        go = list(
            release = backend_env$go_bp_hierarchy_data$metadata$go_release %||%
                backend_env$go_bp_hierarchy_data$metadata$go_db_version %||%
                NA_character_
        ),
        mitocarta = list(
            release = backend_env$mitocarta_hierarchy_data$metadata$mitocarta_release %||%
                NA_character_
        )
    )
}

.fetch_online_backend_versions <- function(databases = c("reactome", "go")) {
    versions <- list(
        reactome = list(release = NA_character_),
        go = list(release = NA_character_),
        mitocarta = list(release = NA_character_)
    )

    if ("reactome" %in% databases) {
        reactome_release <- trimws(.read_first_remote_line(
            "https://reactome.org/ContentService/data/database/version"
        ))

        versions$reactome$release <- reactome_release
    }

    if ("go" %in% databases) {
        go_header <- .read_remote_lines(
            "https://current.geneontology.org/ontology/go-basic.obo",
            n = 40
        )

        go_release_line <- go_header[grepl("^data-version: ", go_header)][1]

        if (length(go_release_line) == 0 || is.na(go_release_line)) {
            stop("Could not determine the current GO release from go-basic.obo.", call. = FALSE)
        }

        versions$go$release <- sub("^data-version: ", "", go_release_line)
    }

    if ("mitocarta" %in% databases) {
        mitocarta_page_lines <- .read_remote_lines(
            "https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways",
            n = 220
        )

        release_line <- mitocarta_page_lines[grepl("MitoCarta3\\.0, released 2020", mitocarta_page_lines)][1]

        if (length(release_line) == 0 || is.na(release_line)) {
            stop("Could not determine the current MitoCarta release from the Broad website.", call. = FALSE)
        }

        versions$mitocarta$release <- "MitoCarta3.0 (2020)"
    }

    versions
}

.download_reactome_backend_files <- function(out_dir, reactome_release) {
    relation_path <- file.path(out_dir, "ReactomePathwaysRelation.txt")
    pathway_path <- file.path(out_dir, "ReactomePathways.txt")

    .download_remote_file(
        url = "https://download.reactome.org/current/ReactomePathwaysRelation.txt",
        destfile = relation_path
    )

    .download_remote_file(
        url = "https://download.reactome.org/current/ReactomePathways.txt",
        destfile = pathway_path
    )

    list(
        relation_path = relation_path,
        pathway_path = pathway_path,
        reactome_release = reactome_release
    )
}

.download_go_obo <- function(out_dir, go_release) {
    obo_path <- file.path(out_dir, "go-basic.obo")

    .download_remote_file(
        url = "https://current.geneontology.org/ontology/go-basic.obo",
        destfile = obo_path
    )

    list(
        obo_path = obo_path,
        go_release = go_release
    )
}

.download_mitocarta_html <- function(out_dir, mitocarta_release) {
    html_path <- file.path(out_dir, "human.mitopathways3.0.html")

    .download_remote_file(
        url = "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/human.mitopathways3.0.html",
        destfile = html_path
    )

    list(
        html_path = html_path,
        mitocarta_release = mitocarta_release
    )
}

.build_reactome_backend_from_downloads <- function(package_root, relation_path, pathway_path, reactome_release) {
    script_env <- new.env(parent = globalenv())
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(package_root)
    source(file.path(package_root, "data-raw", "build_reactome_hierarchy.R"), local = script_env)

    script_env$build_reactome_hierarchy_data_from_files(
        reactome_relation_path = relation_path,
        reactome_pathway_path = pathway_path,
        reactome_release = reactome_release,
        source_label = "reactome_download_current"
    )
}

.build_go_backend_from_obo <- function(package_root, obo_path, ontology, go_release) {
    script_env <- new.env(parent = globalenv())
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(package_root)
    source(file.path(package_root, "data-raw", "build_go_hierarchies.R"), local = script_env)

    script_env$build_go_hierarchy_data_from_obo(
        obo_path = obo_path,
        ontology = ontology,
        go_release = go_release
    )
}

.build_mitocarta_backend_from_html <- function(package_root, html_path, mitocarta_release) {
    script_env <- new.env(parent = globalenv())
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(package_root)
    source(file.path(package_root, "data-raw", "build_mitocarta_hierarchy.R"), local = script_env)

    script_env$build_mitocarta_hierarchy_data(
        source_html_path = html_path,
        mitocarta_release = mitocarta_release
    )
}

.load_existing_backend_objects <- function(package_root) {
    backend_env <- new.env(parent = emptyenv())
    load(file.path(package_root, "R", "sysdata.rda"), envir = backend_env)

    list(
        reactome_hierarchy_data = backend_env$reactome_hierarchy_data,
        go_bp_hierarchy_data = backend_env$go_bp_hierarchy_data,
        go_mf_hierarchy_data = backend_env$go_mf_hierarchy_data,
        go_cc_hierarchy_data = backend_env$go_cc_hierarchy_data,
        mitocarta_hierarchy_data = backend_env$mitocarta_hierarchy_data,
        mitocarta_term2gene_tbl = backend_env$mitocarta_term2gene_tbl,
        mitocarta_term2name_tbl = backend_env$mitocarta_term2name_tbl
    )
}

.print_backend_version_summary <- function(local_versions, online_versions, rebuild_reactome, rebuild_go, rebuild_mitocarta) {
    message("Reactome backend version:")
    message("  local:  ", local_versions$reactome$release %||% "unknown")
    message("  online: ", online_versions$reactome$release %||% "unavailable")
    message("  action: ", if (isTRUE(rebuild_reactome)) "rebuild" else "up to date")

    message("GO backend version:")
    message("  local:  ", local_versions$go$release %||% "unknown")
    message("  online: ", online_versions$go$release %||% "unavailable")
    message("  action: ", if (isTRUE(rebuild_go)) "rebuild" else "up to date")

    message("MitoCarta backend version:")
    message("  local:  ", local_versions$mitocarta$release %||% "unknown")
    message("  online: ", online_versions$mitocarta$release %||% "unavailable")
    message("  action: ", if (isTRUE(rebuild_mitocarta)) "rebuild" else "up to date")
}

.read_first_remote_line <- function(url) {
    lines <- .read_remote_lines(url = url, n = 1)

    if (length(lines) == 0) {
        stop("No content was returned from: ", url, call. = FALSE)
    }

    lines[[1]]
}

.read_remote_lines <- function(url, n = -1L) {
    connection <- url(url, open = "rb")
    on.exit(close(connection), add = TRUE)
    readLines(connection, n = n, warn = FALSE, encoding = "UTF-8")
}

.download_remote_file <- function(url, destfile) {
    utils::download.file(
        url = url,
        destfile = destfile,
        mode = "wb",
        quiet = TRUE
    )

    if (!file.exists(destfile) || file.info(destfile)$size <= 0) {
        stop("Download failed for: ", url, call. = FALSE)
    }

    invisible(destfile)
}

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
        return(y)
    }

    x
}
