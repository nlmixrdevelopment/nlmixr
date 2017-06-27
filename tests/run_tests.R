### run tests

Sys.setenv(NLMIXR_VALIDATION_FULL="true")
Sys.setenv(NLMIXR_VALIDATION="true")
testthat::test_dir("E:/Occams/Local/General/nlmixr/nlmixr/tests/testthat/")

