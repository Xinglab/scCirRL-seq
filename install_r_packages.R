repos <- "http://cran.us.r-project.org"

utils::install.packages('remotes', repos=repos)
base::stopifnot(base::require('remotes'))

utils::install.packages('Seurat', repos=repos)
base::stopifnot(base::require('Seurat'))

remotes::install_github("satijalab/azimuth")
base::stopifnot(base::require('Azimuth'))
