library(testthat)
library(RxODE)
library(nlmixr)
verbose_minimization <- FALSE

test_check("nlmixr", stop_on_failure = FALSE, wrap=TRUE,
           reporter = testthat::LocationReporter)


