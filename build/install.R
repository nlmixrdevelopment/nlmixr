Sys.setenv(PATH=paste0("c:/Rtools/bin;c:/Rtools/mingw_64/bin;", Sys.getenv("PATH")))
install.packages(c("tidyverse", "data.table",
                   "devtools", "Rcpp", "brew",
                   "lattice","lbfgs",
                   "inline","dparser","ggplot2",
                   "rex","minqa","Matrix",
                   "numDeriv","R.utils","mvtnorm","n1qn1",
                   "PreciseSums","fastGHQuad","crayon",
                   "cli", "RcppArmadillo","dplyr","ggforce",
                   "purrr","readr","rlang","stringr",
                   "tibble","tidyr","gridExtra","rmarkdown","knitr",
                   "testthat","plotly","webshot","mvtnorm",
                   "shiny",
                   "rmarkdown",
                   "tidyr",
                   "tibble",
                   "curl",
                   "ggplot2",
                   "gridExtra",
                   "microbenchmark",
                   "scales",
                   "stringi", "lbfgsb3",
                   "lbfgsb3c", "madness", "expm", "matrixcalc", "bookdown", "roxygen2", "xpose",
                   ## "reticulate",
                   "nloptr", "ucminf", "vpc", "installr", "DT", "dotwhisker", "broom", "broom.mixed",
                   "Rvmmin", "pkgdown", "xgxr", "RxODE", "nlmixr", "ggPMX"), type="source")

## devtools::install_github("nlmixrdevelopment/rxModels")
library(RxODE)
devtools::install_github("nlmixrdevelopment/SnakeCharmR")
devtools::install_github("richardhooijmaijers/R3port")
devtools::install_github("nlmixrdevelopment/xpose.nlmixr")

devtools::install_github("AdeelK93/collapsibleTree")

devtools::install_github("richardhooijmaijers/shinyMixR")




