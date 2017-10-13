##' Validate nlmixr
##'
##' This allows easy vaildation/qualification of nlmixr by running the
##' testing suite on your system.
##' @param full Should a full validation be performed?  (By default
##'     \code{FALSE})
##' @author Matthew L. Fidler
##' @export
nlmixrValidate <- function(full=FALSE){
    Sys.setenv("nlmixr_silent"="TRUE")
    on.exit({Sys.setenv("nlmixr_silent"="")})
    if (full){
        Sys.setenv(NLMIXR_VALIDATION_FULL="true")
    }
    Sys.setenv(NLMIXR_VALIDATION="true")
    old.wd <- getwd();
    on.exit({setwd(old.wd)});
    path <- file.path(system.file("tests", package = "nlmixr"),"testthat")
    setwd(path)
    testthat::test_dir(path);
}
