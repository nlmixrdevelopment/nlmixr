#' Combines theta/eta parameters needed for residuals/shrinkage calculations
#'
#' @param fit focei style fit
#'
#' @return list with:
#'
#'  - A RxODE `params` dataset `pred` predictions
#'
#'  - A RxODE `params` dataset for `ipred` predictions
#'
#'  - `eta.lst` is a numerical vector for each of the ETAs listed. The first 5 components are the mean, sd, variance, kurtosis and
#'    skewness statistics.  The rest of the components will be filled in later when calculating the shrinkage dataframe
#'
#' @author Matthew Fidler
#' @noRd
.foceiThetaEtaParameters <- function(fit) {
  .etas <- fit$ranef
  .thetas <- fit$fixef
  .Call(`_nlmixr_nlmixrParameters`, .thetas, .etas)
}

#' Solve making sure that ID is dropped
#'
#' This also suppresses the warning for id sorting
#' @param ... All parameters set to `RxODE::rxSolve()`
#' @param fit focei style fit
#' @return solved tataset
#' @author Matthew Fidler
#' @noRd
.foceiSolveWithId <- function(...) {
  .ret <- RxODE::rxSolve(..., warnIdSort = FALSE)
  if (names(.ret)[1] == "time") {
    ## For single subject ID is dropped.
    .ret <- data.frame(ID = 1, .ret)
  }
  return(.ret)
}
#' Solve for pred/ipred types of calculations (including residuals)
#'
#' @param fit focei style fit
#' @param model RxODE model
#' @param pars parameters to solve
#' @param keep vector of columns to keep
#' @param what string of what type of calculation is being performed
#' @inheritParams RxODE::rxSolve
#' @return Solved RxODE data
#' @author Matthew Fidler
#' @noRd
.foceiSolvePars <- function(fit, model, pars=NULL, returnType="data.frame", keep=NULL, what="pred",
                            addDosing=FALSE, subsetNonmem=TRUE) {
  .res <- .foceiSolveWithId(model, pars, fit$dataSav,
                            returnType = returnType,
                            atol = fit$control$atol[1], rtol = fit$control$rtol[1],
                            maxsteps = fit$control$maxstepsOde,
                            hmin = fit$control$hmin, hmax = fit$control$hmax, hini = fit$control$hini,
                            transitAbs = fit$control$transitAbs, maxordn = fit$control$maxordn,
                            maxords = fit$control$maxords, method = fit$control$method,
                            keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem
                            )
  RxODE::rxSolveFree()
  if (any(is.na(.res$rx_pred_)) && fit$control$method == 2L) {
    .res <- .foceiSolveWithId(model, pars, fit$dataSav,
                              returnType = returnType,
                              atol = fit$control$atol[1], rtol = fit$control$rtol[1],
                              maxsteps = fit$control$maxstepsOde * 2,
                              hmin = fit$control$hmin, hmax = fit$control$hmax / 2, hini = fit$control$hini,
                              transitAbs = fit$control$transitAbs, maxordn = fit$control$maxordn,
                              maxords = fit$control$maxords, method = "lsoda",
                              keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem)
    RxODE::rxSolveFree()
    if (any(is.na(.res$rx_pred_))) {
      .res <- .foceiSolveWithId(model, pars, fit$dataSav,
                                returnType = returnType,
                                atol = fit$control$atol[1], rtol = fit$control$rtol[1],
                                maxsteps = fit$control$maxstepsOde * 2,
                                hmin = fit$control$hmin, hmax = fit$control$hmax / 2, hini = fit$control$hini,
                                transitAbs = fit$control$transitAbs, maxordn = fit$control$maxordn,
                                maxords = fit$control$maxords, method = "dop853",
                                keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem)
      RxODE::rxSolveFree()
      if (any(is.na(.res$rx_pred_))) {
        warning("Problems solving ", what, " liblsoda, lsoda and dop853")
      } else {
        warning("Problems solving ", what, " liblsoda and lsoda switched to dop853")
      }
    } else {
      warning("Problems solving ", what, " liblsoda switched to lsoda")
    }
  }
  return(.res)
}

#' Create a ipred/pred list from the focei style model
#'
#' @param fit focei style fit
#' @param thetaEtaParameters Theta/eta parameter list generated from `.foceiThetaEtaParameters()`
#' @param predOnly Pred Only for .ipred model (useful for mean/population models)
#' @inheritParams RxODE::rxSolve
#' @return list with ipred and pred datasets
#' @author Matthew Fidler
#' @noRd
.foceiPredIpredList <- function(fit, thetaEtaParameters=.foceiThetaEtaParameters(fit), keep=NULL, predOnly=!is.null(fit$model$inner),
                                addDosing=FALSE, subsetNonmem=TRUE) {
  .ipredModel <- fit$model$inner
  if (predOnly) .ipredModel <- fit$model$pred.only
  list(ipred = .foceiSolvePars(fit, .ipredModel, thetaEtaParameters$ipred,
                               returnType="data.frame.TBS", keep=keep, what="ipred",
                               addDosing=addDosing, subsetNonmem=subsetNonmem),
       pred = .foceiSolvePars(fit, .ipredModel, thetaEtaParameters$pred,returnType="data.frame", what="pred",
                              addDosing=addDosing, subsetNonmem=subsetNonmem),
       etaLst=thetaEtaParameters$eta.lst)
}

.getRelevantLhs <- function(fit, keep=NULL) {
  .ret <- setdiff(fit$model$pred.only$lhs,fit$ini$name)
  .w <- which(regexpr("^rx", .ret) == -1)
  unique(c(.ret[.w], keep))
}

.calcCwres <- function(fit, data=fit$dataSav, thetaEtaParameters=.foceiThetaEtaParameters(fit),
                       table=tableControl(), dv=NULL, predOnly=FALSE,
                       addDosing=FALSE, subsetNonmem=TRUE, keep=NULL) {
  if (!inherits(table, "tableControl")) table <- do.call(tableControl, table)
  .keep <- keep
  .names <- names(data)
  .lowerNames <- tolower(.names)
  for (.n in c("dv", "cens", "limit")) {
    .w <- which(.lowerNames == .n)
    if (length(.w) == 1L) .keep <- c(.keep, .names[.w])
  }
  .keep <- unique(.keep)
  .prdLst <- .foceiPredIpredList(fit, keep=.keep, thetaEtaParameters=thetaEtaParameters, predOnly=predOnly,
                                 addDosing=addDosing, subsetNonmem=subsetNonmem)
  if (!inherits(dv, "numeric")) {
    dv <- .prdLst$ipred$dv
    table$doSim <- TRUE
  } else {
    table$doSim <- FALSE
  }
  if (predOnly){
    .Call(`_nlmixr_resCalc`, .prdLst, fit$omega,
          fit$eta, .prdLst$ipred$dv, .prdLst$ipred$evid, .prdLst$ipred$cens,
          .prdLst$ipred$limit, .getRelevantLhs(fit, keep), fit$model$pred.only$state, table)
  } else {
    .Call(`_nlmixr_cwresCalc`, .prdLst, fit$omega,
          fit$eta, .prdLst$ipred$dv, .prdLst$ipred$evid, .prdLst$ipred$cens,
          .prdLst$ipred$limit, .getRelevantLhs(fit, keep), fit$model$pred.only$state, table)
  }
}

.calcRes <- function(..., predOnly=TRUE) {
  .calcCwres(..., predOnly=predOnly)
}

.calcIres <- function(fit, data=fit$dataSav, table=tableControl(), dv=NULL,
                      addDosing=FALSE, subsetNonmem=TRUE) {
  if (!inherits(table, "tableControl")) table <- do.call(tableControl, table)
  .keep <- NULL
  .names <- names(data)
  .lowerNames <- tolower(.names)
  for (.n in c("dv", "cens", "limit")) {
    .w <- which(.lowerNames == .n)
    if (length(.w) == 1L) .keep <- c(.keep, .names[.w])
  }
  .thetas <- fit$fixef
  names(.thetas) <- paste0("THETA[", seq_along(.thetas), "]")
  .eta <- fit$eta
  if (inherits(.eta, "data.frame")) {
    .n <- length(.eta) - 1
    .thetas <- c(.thetas, setNames(rep(0, .n), paste0("ETA[", seq_len(.n), "]")))
  }
  .ipred <- .foceiSolvePars(fit, fit$model$pred.only, .thetas,
                            returnType="data.frame.TBS", keep=.keep, what="ipred",
                            addDosing=addDosing, subsetNonmem=subsetNonmem)
  if (!inherits(dv, "numeric")) {
    dv <- .ipred$dv
    table$doSim <- TRUE
  } else {
    table$doSim <- FALSE
  }
  .Call(`_nlmixr_iresCalc`, .ipred, dv, .ipred$evid, .ipred$cens, .ipred$limit, table)
}

.calcShrinkOnly <- function(fit, thetaEtaParameters=.foceiThetaEtaParameters(fit)) {
  .omega <- fit$omega
  .ret <- .Call(`_nlmixr_calcShrinkOnly`, .omega, thetaEtaParameters$eta.lst, length(fit$eta[,1]))
  .ret[, -dim(.omega)[1] - 1]
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
##' @param censMethod Handle censoring method:
##'
##'  - `"truncated-normal"` Simulates from a truncated normal distribution under the assumption of the model and censoring.
##'
##'  - `"cdf"` Use the cdf-method for censoring with npde and use this for any other residuals (`cwres` etc)
##'
##'  - `"omit"` omit the residuals for censoring
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
                         nsim = 300, ties = TRUE,
                         censMethod=c("truncated-normal", "cdf", "omit"),
                         seed = 1009) {
  .ret <- list(
    npde = npde, cwres = cwres, saemNPDE = saemNPDE,
    saemCWRES = saemCWRES, nlmeNPDE = nlmeNPDE,
    nlmeCWRES = nlmeCWRES, foceiNPDE = foceiNPDE,
    foceNPDE = foceNPDE, nsim = nsim, ties = ties, seed = seed,
    censMethod=setNames(c("truncated-normal"=3L, "cdf"=2L, "omit"=1L)[match.arg(censMethod)], NULL)
  )
  class(.ret) <- "tableControl"
  return(.ret)
}
