##' Validate nlmixr
##'
##' This allows easy vaildation/qualification of nlmixr by running the
##' testing suite on your system.
##' @param full Should a full validation be performed?  (By default
##'     \code{FALSE})
##' @author Matthew L. Fidler
##' @export
nlmixrValidate <- function(full=FALSE){
    ## rxVersion(" Validation", TRUE);
    old.wd <- getwd();
    on.exit({setwd(old.wd); Sys.setenv(NLMIXR_VALIDATION="false", NLMIXR_VALIDATION_FULL="false", NOT_CRAN="")});
    if (is.character(full)){
        Sys.setenv(NLMIXR_VALIDATION_FULL="false", "NOT_CRAN"="true")
        path <- file.path(system.file("tests", package = "nlmixr"),"testthat")
        setwd(path)
        testthat::test_dir(path, filter=full);
        Sys.setenv(NLMIXR_VALIDATION="true")
        testthat::test_dir(path, filter=full);
        Sys.setenv(NLMIXR_VALIDATION_FULL="true")
        testthat::test_dir(path, filter=full);
    } else {
        old.wd <- getwd();
        path <- file.path(system.file("tests", package = "nlmixr"),"testthat")
        setwd(path)
        Sys.setenv("NOT_CRAN"="funny")
        message("CRAN only tests")
        message("================================================================================")
        pt <- proc.time();
        testthat::test_dir(path);
        message("================================================================================")
        message("Timing of CRAN tests (should be under 60 seconds)")
        message("================================================================================")
        print(proc.time() - pt);
        message("================================================================================")
        message("Normal tests")
        message("================================================================================")
        Sys.setenv("NOT_CRAN"="true")
        testthat::test_dir(path);
        if (full){
            message("================================================================================")
            message("Validation tests")
            message("================================================================================")
            Sys.setenv("NLMIXR_VALIDATION"="true")
            testthat::test_dir(path);
            message("================================================================================")
            message("Full Validation tests")
            message("================================================================================")
            Sys.setenv("NLMIXR_VALIDATION"="false")
            Sys.setenv("NLMIXR_VALIDATION_FULL"="true")
            testthat::test_dir(path);
        }
    }
}
