##' Validate nlmixr
##'
##' This allows easy vaildation/qualification of nlmixr by running the
##' testing suite on your system.
##' @param type of test to be run
##' @param check Use devtools::check to run checks
##' @author Matthew L. Fidler
##' @export
nlmixrValidate <- function(type = NULL, check = FALSE) {
  pt <- proc.time()
  .filter <- NULL
  if (is.null(type)) type <- FALSE
  if (is.character(type)) {
    .filter <- type
    type <- TRUE
  }
  if (type == TRUE) {
    .oldCran <- Sys.getenv("NOT_CRAN")
    Sys.setenv("NOT_CRAN"="true")
    on.exit(Sys.setenv("NOT_CRAN"=.oldCran))
  }
  RxODE::.rxWithOptions(list(testthat.progress.max_fails=10000000000), {
    path <- file.path(system.file("tests", package = "nlmixr"), "testthat")
    RxODE::.rxWithWd(path, {
      try(testthat::test_dir(path, filter = .filter))
      message("================================================================================")
      print(proc.time() - pt)
      message("================================================================================")
    })
  })
}

##' @rdname nlmixrValidate
##' @export
nmTest <- nlmixrValidate
