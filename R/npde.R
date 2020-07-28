##' NPDE calculation for nlmixr
##'
##' @param object nlmixr fit object
##' @param nsim Number of simulations.  By default this is 300
##' @param ties When TRUE, the npde distribution can have ties.  When
##'     FALSE, the npde distribution uses uniform random numbers to
##'     prevent ties.
##' @param seed Seed for running nlmixr simulation.  By default 1009
##' @param updateObject Boolean indicating if original object should be updated.  By default this is TRUE.
##' @inheritParams foceiControl
##' @param ... Other ignored parameters.
##' @return New nlmixr fit object
##' @author Matthew L. Fidler
##' @export
addNpde <- function(object, nsim = 300, ties = TRUE, seed = 1009, updateObject = TRUE,
                    cholSEtol = (.Machine$double.eps)^(1 / 3), ...) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  .pt <- proc.time()
  .objName <- substitute(object)
  if (any(names(object) == "NPDE")) {
    warning("Already contains NPDE")
    return(object)
  }
  set.seed(seed)
  .si <- object$simInfo
  .rx <- .si$rx
  .rx <- gsub(rex::rex(capture("ipred"), or("=", "~"), capture(except_any_of("\n;")), any_of("\n;")), "ipred~\\2;\n", .rx)
  .rx <- gsub(rex::rex("d/dt(", capture(except_any_of("\n;)")), ")", or("=", "~")), "d/dt(\\1)~", .rx)
  .rx <- gsub(
    rex::rex("sim", or("=", "~"), "rxTBSi(", capture(except_any_of(",)")), ",", anything, any_of("\n;")),
    "sim=\\1", .rx
  )
  .si$rx <- .rx
  .dat <- nlmixrData(.nmGetData(object))
  .dat <- .dat[.dat$EVID == 0, ]
  .si$object <- object
  .si$returnType <- "data.frame.TBS"
  .si$nsim <- nsim
  .si <- c(.si, list(...))
  .si$modelName <- "NPDE"
  .pt <- proc.time()
  .si$dfObs <- 0
  .si$dfSub <- 0
  .si$thetaMat <- NA
  .sim <- do.call("nlmixrSim", .si)
  .dv <- object$DV
  .dvl <- length(.dv)
  .cls <- class(object)
  .evid <- rep(0L, .dvl)
  .evid[is.na(object$RES) & !is.na(object$PRED)] <- 2L
  .new <- cbind(object, .Call(
    `_nlmixr_npde`, object$ID, .dv, .evid, .sim$sim, .sim$rxLambda, .sim$rxYj, ties,
    cholSEtol
  ))
  class(.new) <- .cls
  if (updateObject) {
    .parent <- parent.frame(2)
    .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
      if (.cur == .objName && identical(.parent[[.cur]], object)) {
        return(.cur)
      }
      return(NULL)
    }))
    if (length(.bound) == 1) {
      if (exists(.bound, envir = .parent)) {
        assign(.bound, .new, envir = .parent)
      }
    }
  }
  .env <- .new$env
  .env$time <- .data.frame(.env$time, npde = (proc.time() - .pt)["elapsed"], check.names = FALSE)
  return(.new)
}
##' Output table/data.frame options
##'
##' @param npde When TRUE, request npde regardless of the algorithm used.
##'
##' @param cwres When TRUE, request CWRES and FOCEi likelihood
##'     regardless of the algorithm used.
##'
##' @param saemNPDE When TRUE and estimating with SAEM, adds NPDE
##'     metrics to fit including EPRED, ERES, and NPDE. (default
##'     TRUE);
##'
##' @param saemCWRES When TRUE and estimating with SAEM, adds CWRES
##'     metrics to the fit including CPRED, CRES and CWRES.  It also
##'     evaluates the function with the FOCEi objective function to
##'     allow comparison between estimation methods. (default FALSE)
##'
##' @param nlmeNPDE When TRUE and estimating with nlme, adds NPDE
##'     metrics to fit including EPRED, ERES, and NPDE. (default
##'     TRUE);
##'
##' @param nlmeCWRES When TRUE and estimating with nlme, adds CWRES
##'     metrics to the fit including CPRED, CRES and CWRES.  It also
##'     evaluates the function with the FOCEi objective function to
##'     allow comparison between estimation methods. (default FALSE)
##'
##' @param foceiNPDE When TRUE and estimating with FOCEi, adds NPDE
##'     metrics to fit including EPRED, ERES, and NPDE. (default
##'     TRUE);
##'
##' @param foceNPDE When TRUE and estimating with FOCEi, adds NPDE
##'     metrics to fit including EPRED, ERES, and NPDE. (default
##'     TRUE);
##'
##' @inheritParams addNpde
##'
##' @details
##'
##' If you ever want to add CWRES/FOCEi objective function you can use the \code{\link{addCwres}}
##'
##' If you ever want to add NPDE/EPRED columns you can use the \code{\link{addNpde}}
##'
##' @return A list of table options for nlmixr
##' @author Matthew L. Fidler
##' @export
tableControl <- function(npde = NULL,
                         cwres = NULL,
                         saemNPDE = FALSE,
                         saemCWRES = FALSE,
                         nlmeNPDE = FALSE,
                         nlmeCWRES = FALSE,
                         foceiNPDE = FALSE,
                         foceNPDE = FALSE,
                         nsim = 300, ties = TRUE, seed = 1009) {
  .ret <- list(
    npde = npde, cwres = cwres,
    saemNPDE = saemNPDE,
    saemCWRES = saemCWRES,
    nlmeNPDE = nlmeNPDE,
    nlmeCWRES = nlmeCWRES,
    foceiNPDE = foceiNPDE,
    foceNPDE = foceNPDE,
    nsim = nsim, ties = ties, seed = seed
  )
  class(.ret) <- "tableControl"
  return(.ret)
}
