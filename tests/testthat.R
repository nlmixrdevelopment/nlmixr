Sys.setenv("R_TESTS" = "")
library(testthat)
library(nlmixr)

test_check("nlmixr")
