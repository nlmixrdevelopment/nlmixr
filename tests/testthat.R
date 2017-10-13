Sys.setenv("R_TESTS" = "")
Sys.setenv("nlmixr_silent"="TRUE")
library(testthat)
library(nlmixr)

test_check("nlmixr")
