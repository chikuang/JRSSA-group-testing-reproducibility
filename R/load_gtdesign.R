# Load the author's installed package, rather than sourcing its implementation.
# Run scripts from the repository root after scripts/install_dependencies.R.
load_gtdesign <- function(source_dir = "vendor/gtDesign") {
  local_library <- normalizePath(".library", mustWork = FALSE)
  .libPaths(c(local_library, .libPaths()))
  if (!dir.exists(file.path(local_library, "gtDesign"))) {
    stop("Install the archived package first: Rscript scripts/install_dependencies.R")
  }
  expected <- unname(read.dcf(file.path(source_dir, "DESCRIPTION"))[1, "Version"])
  actual <- as.character(utils::packageVersion("gtDesign", lib.loc = local_library))
  if (!identical(actual, expected)) {
    stop("Installed gtDesign differs from the archived version; rerun scripts/install_dependencies.R")
  }
  suppressPackageStartupMessages(library("gtDesign", lib.loc = local_library))
  invisible(actual)
}
