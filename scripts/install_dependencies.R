# Run from the repository root. Existing package versions are not silently changed.
pkgs <- c("CVXR","tibble","MASS","dplyr","tidyr","purrr","ggplot2","forcats","ggrepel","knitr","rmarkdown","numDeriv","jsonlite")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) install.packages(missing,repos="https://cloud.r-project.org")
dir.create(".library",showWarnings=FALSE)
status <- system2(file.path(R.home("bin"),"R"),c("CMD","INSTALL","--library=.library","vendor/gtDesign"))
if(status!=0) stop("Installation of the vendored gtDesign source failed")
cat("Installed gtDesign locally. Recorded validation versions: dependencies.csv\n")
