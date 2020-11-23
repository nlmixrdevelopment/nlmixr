
Sys.setenv(PATH=paste0("c:/Rtools/bin;c:/Rtools/mingw_64/bin;", Sys.getenv("PATH")))
install.packages(c("assertthat", "backports", "BH", "brew", "broom.mixed", "checkmate",
                   "cli", "covr", "crayon", "curl", "data.table", "Deriv", "devtools",
                   "digest", "dotwhisker", "dparser", "dplyr", "DT", "expm", "fastGHQuad",
                   "flextable", "generics", "ggplot2", "ggrepel", "ggtext", "gridExtra",
                   "htmltools", "huxtable", "inline", "knitr", "lattice", "lbfgs",
                   "lbfgsb3c", "learnr", "lotri", "madness", "magrittr",
                   "matrixcalc", "memoise", "methods", "microbenchmark", "minqa",
                   "n1qn1", "nlme", "nloptr", "officer", "pillar", "pkgdown",
                   "PreciseSums", "Rcpp",  "RcppEigen", "remotes",
                   "reshape2", "rex", "rlang", "rmarkdown", "roxygen2", "Rvmmin",
                   "scales", "shiny", "sitmo", "stringi", "symengine",
                   "sys", "testthat", "tibble", "tidyr", "tidyverse", "tools", "ucminf",
                   "units", "usethis", "utils", "vdiffr", "vpc", "xgxr", "xpose",
                   "yaml",
                   "RcppArmadillo", "StanHeaders") , type="source",
                 # For R 3.6.1 use the release date of R 3.6.3
                 #repos="https://cran.microsoft.com/snapshot/2020-02-29/"
                 )

#devtools::install_github("nlmixrdevelopment/n1qn1c", dependencies = FALSE)

devtools::install_github("nlmixrdevelopment/RxODE", ref="pruneBranch", dependencies = FALSE)

devtools::install_github("nlmixrdevelopment/nlmixr", ref="foceiDur", dependencies = FALSE)

## devtools::install_github("nlmixrdevelopment/rxModels")
library(RxODE)
devtools::install_github("richardhooijmaijers/R3port")
devtools::install_github("nlmixrdevelopment/xpose.nlmixr")

devtools::install_github("AdeelK93/collapsibleTree")

devtools::install_github("richardhooijmaijers/shinyMixR")




