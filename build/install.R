Sys.setenv(PATH=paste0("c:/Rtools/bin;c:/Rtools/mingw_64/bin;", Sys.getenv("PATH")))
install.packages(c("assertthat", "backports", "BH", "brew", "broom.mixed", "checkmate",
                   "cli", "covr", "crayon", "curl", "data.table", "Deriv", "devtools",
                   "digest", "dotwhisker", "dparser", "dplyr", "DT", "expm", "fastGHQuad",
                   "flextable", "generics", "ggplot2", "ggrepel", "ggtext", "gridExtra",
                   "htmltools", "huxtable", "inline", "knitr", "lattice", "lbfgs",
                   "lbfgsb3c", "learnr", "lotri", "madness", "magrittr",
                   "matrixcalc", "memoise", "methods", "microbenchmark", "minqa",
                   "n1qn1", "nlme", "nloptr", "officer", "parallel", "pillar", "pkgdown",
                   "PreciseSums", "Rcpp", "RcppArmadillo", "RcppEigen", "remotes",
                   "reshape2", "rex", "rlang", "rmarkdown", "roxygen2", "Rvmmin",
                   "scales", "shiny", "sitmo", "StanHeaders", "stringi", "symengine",
                   "sys", "testthat", "tibble", "tidyr", "tidyverse", "tools", "ucminf",
                   "units", "usethis", "utils", "vdiffr", "vpc", "xgxr", "xpose",
                   "yaml") , type="source")

devtools::install_github("nlmixrdevelopment/RxODE", ref="pruneBranch")
devtools::install_github("nlmixrdevelopment/nlmixr", ref="foceiDur")

## devtools::install_github("nlmixrdevelopment/rxModels")
library(RxODE)
devtools::install_github("richardhooijmaijers/R3port")
devtools::install_github("nlmixrdevelopment/xpose.nlmixr")

devtools::install_github("AdeelK93/collapsibleTree")

devtools::install_github("richardhooijmaijers/shinyMixR")




