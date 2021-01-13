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
.foceiSolvePars <- function(fit, model, pars=NULL, returnType="data.frame", keep=NULL, what="pred") {
  .res <- .foceiSolveWithId(model, pars, fit$dataSav,
                            returnType = returnType,
                            atol = fit$control$atol[1], rtol = fit$control$rtol[1],
                            maxsteps = fit$control$maxstepsOde,
                            hmin = fit$control$hmin, hmax = fit$control$hmax, hini = fit$control$hini,
                            transitAbs = fit$control$transitAbs, maxordn = fit$control$maxordn,
                            maxords = fit$control$maxords, method = fit$control$method,
                            keep=keep
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
                              keep=keep)
    RxODE::rxSolveFree()
    if (any(is.na(.res$rx_pred_))) {
      .res <- .foceiSolveWithId(model, pars, fit$dataSav,
                                returnType = returnType,
                                atol = fit$control$atol[1], rtol = fit$control$rtol[1],
                                maxsteps = fit$control$maxstepsOde * 2,
                                hmin = fit$control$hmin, hmax = fit$control$hmax / 2, hini = fit$control$hini,
                                transitAbs = fit$control$transitAbs, maxordn = fit$control$maxordn,
                                maxords = fit$control$maxords, method = "dop853",
                                keep=keep)
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
.foceiPredIpredList <- function(fit, thetaEtaParameters=.foceiThetaEtaParameters(fit), keep=NULL, predOnly=!is.null(fit$model$inner)) {
  .ipredModel <- fit$model$inner
  if (predOnly) .ipredModel <- fit$model$pred.only
  list(ipred = .foceiSolvePars(fit, .ipredModel, thetaEtaParameters$ipred,
                               returnType="data.frame.TBS", keep=keep, what="ipred"),
        pred = .foceiSolvePars(fit, .ipredModel, thetaEtaParameters$pred,returnType="data.frame", what="pred"))
}
