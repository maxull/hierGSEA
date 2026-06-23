# hierGSEA news

## hierGSEA 0.1.0-beta.1

First public beta release.

### Added

- Reactome hierarchy-aware post-processing of upstream `gseaResult` objects
- GO hierarchy-aware post-processing for separate `BP`, `MF`, and `CC` ontologies
- MitoCarta hierarchical backend using Broad `MitoPathways3.0`
- branch-local Benjamini-Hochberg correction using visible hierarchy families
- branch retention logic that keeps significant descendants and required ancestors
- hierarchy-aware `ggplot2` visualization with left-side tree structure
- bundled example workflows for single-fiber transcriptomics and HIRC proteomics
- backend update helper for Reactome, GO, and MitoCarta
- pkgdown site scaffold and long-form vignettes

### Notes

- This is a beta release and the plotting API may still evolve as more real
  datasets are tested.
- The package is currently GitHub-first and is not yet prepared for CRAN or
  Bioconductor submission.
