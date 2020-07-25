##' Validate nlmixr
##'
##' This allows easy vaildation/qualification of nlmixr by running the
##' testing suite on your system.
##' @param full Should a full validation be performed?  (By default
##'     \code{FALSE})
##' @author Matthew L. Fidler
##' @export
nlmixrValidate <- function(full = FALSE) {
  ## rxVersion(" Validation", TRUE);
  if (full == "site") {
    message("Needs to be run on installed nlmixr")
    if (!file.exists("DESCRIPTION")) stop("need to be in package root!")
    .opts <- options()
    on.exit(options(.opts))
    options(
      nlmixr.save = TRUE,
      nlmixr.save.dir = system.file(package = "nlmixr")
    )

    sapply(
      c(
        list.files(system.file(package = "nlmixr"), pattern = "\\.rds$", full.names = TRUE),
        list.files(file.path("inst"), pattern = "\\.rds$", full.names = TRUE)
      ),
      function(x) {
        unlink(x)
      }
    )
    ## For some reason when tangling the nlme runs crash.
    ## run them here first.
    eta.ka <- eta.cl <- eta.v <- eta.ka <- eta.cl <- eta.v <- NULL
    ## Cannot see this anywhere...
    `/<-` <- NULL
    dt <- function(...) {}
    depot <- dt <- center <- NULL

    one.cmt <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45 # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.err <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.err)
      })
    }
    one.cmt <- one.cmt()
    one.compartment <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45 # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.err <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.err)
      })
    }
    one.compartment <- one.compartment()
    fit <- nlmixr(one.cmt, nlmixr::theo_sd, est = "nlme")
    fit4 <- nlmixr(one.compartment, nlmixr::theo_sd, est = "nlme", control = nlmeControl(pnlsTol = .5))
    devtools::build_vignettes(quiet = FALSE, install = FALSE)
    sapply(
      list.files(system.file(package = "nlmixr"), pattern = "\\.rds$", full.names = TRUE),
      function(x) {
        message("\t", x)
        file.copy(x, file.path("inst", basename(x)))
      }
    )
    pkgdown::build_site()
  } else {
    old.wd <- getwd()
    on.exit({
      setwd(old.wd)
      Sys.setenv(NLMIXR_VALIDATION = "false", NLMIXR_VALIDATION_FULL = "false", NOT_CRAN = "")
    })
    if (is.character(full)) {
      Sys.setenv(NLMIXR_VALIDATION_FULL = "false", "NOT_CRAN" = "true")
      path <- file.path(system.file("tests", package = "nlmixr"), "testthat")
      setwd(path)
      testthat::test_dir(path, filter = full)
      Sys.setenv(NLMIXR_VALIDATION = "true")
      testthat::test_dir(path, filter = full)
      Sys.setenv(NLMIXR_VALIDATION_FULL = "true")
      testthat::test_dir(path, filter = full)
    } else {
      old.wd <- getwd()
      path <- file.path(system.file("tests", package = "nlmixr"), "testthat")
      setwd(path)
      Sys.setenv("NOT_CRAN" = "funny")
      message("CRAN only tests")
      message("================================================================================")
      pt <- proc.time()
      testthat::test_dir(path)
      message("================================================================================")
      message("Timing of CRAN tests (should be under 60 seconds)")
      message("================================================================================")
      print(proc.time() - pt)
      message("================================================================================")
      message("Normal tests")
      message("================================================================================")
      Sys.setenv("NOT_CRAN" = "true")
      testthat::test_dir(path)
      if (full) {
        message("================================================================================")
        message("Validation tests")
        message("================================================================================")
        Sys.setenv("NLMIXR_VALIDATION" = "true")
        testthat::test_dir(path)
        message("================================================================================")
        message("Full Validation tests")
        message("================================================================================")
        Sys.setenv("NLMIXR_VALIDATION" = "false")
        Sys.setenv("NLMIXR_VALIDATION_FULL" = "true")
        testthat::test_dir(path)
      }
    }
  }
}
##' @rdname nlmixrValidate
##' @export
nmTest <- nlmixrValidate
