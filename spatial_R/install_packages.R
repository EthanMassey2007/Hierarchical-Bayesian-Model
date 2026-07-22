# Install R packages required by the spatial INLA workflow.

cran_packages <- c(
  "arrow",
  "data.table",
  "ggplot2",
  "Matrix",
  "sf"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages(
    "INLA",
    repos = c(
      getOption("repos"),
      INLA = "https://inla.r-inla-download.org/R/stable"
    )
  )
}
