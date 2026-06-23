# hierGSEA public beta release checklist

## Repository readiness

- Confirm `README.md` shows the package logo and correct GitHub install command.
- Confirm `NEWS.md` contains the beta release notes.
- Confirm `LICENSE`, `DESCRIPTION`, and `NAMESPACE` are current.
- Confirm `local_preprint/` stays ignored and is not included in the public repo.
- Confirm large temporary output files are not tracked unless intentionally kept as examples.

## Package quality checks

- Run `R CMD check` locally if your environment has Pandoc and the required Suggests packages.
- Confirm example scripts still run end to end on your machine.
- Confirm the vignettes knit successfully.
- Confirm the `hier_gsea_result` documentation matches the real output object.
- Confirm the manual branch-selection plotting bug is fixed in a fresh R session.

## GitHub automation

- Enable GitHub Actions for the repository if prompted.
- Confirm the `R-CMD-check` workflow passes on `main`.
- Confirm the `pkgdown` workflow successfully deploys to GitHub Pages.
- In the GitHub repository settings, enable Pages from GitHub Actions.

## Public presentation

- Add a short repository description on GitHub.
- Add relevant repository topics such as `r`, `bioinformatics`, `gsea`, `reactome`, `gene-ontology`, `omics`, and `bioconductor`.
- Add the pkgdown site URL to the GitHub repository “About” section.
- Decide whether to pin a beta disclaimer at the top of the README or in the first GitHub Release.

## Beta release creation

- Create or confirm the tag `v0.1.0-beta.1`.
- Draft the GitHub Release using `release/v0.1.0-beta.1_release_notes.md`.
- Attach any optional overview figures if you want the release page to double as a project landing page.

## After release

- Test installation from GitHub using `remotes::install_github("maxull/hierGSEA")`.
- Ask one or two collaborators to run the example scripts on their machines.
- Collect feedback on API clarity, plot defaults, and documentation gaps before `0.1.0`.
