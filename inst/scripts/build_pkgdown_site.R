suppressPackageStartupMessages({
    library(pkgdown)
})

################################################################################################################################################
######################################################      BUILD THE hierGSEA PKGDOWN SITE      ###############################################
################################################################################################################################################

# This helper script is intentionally simple so it can be run from the package
# root during active development. It assumes that the documentation source
# files, vignettes, and pkgdown configuration already live in the repository.

package_root <- getwd()

message("Building pkgdown site from: ", package_root)

pkgdown::build_site(pkg = package_root, install = FALSE, preview = FALSE, new_process = FALSE)

message("pkgdown site build completed.")
