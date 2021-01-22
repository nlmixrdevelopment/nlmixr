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
    .ret <- data.frame(ID = 1L, .ret)
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
                            addDosing=FALSE, subsetNonmem=TRUE, addCov=FALSE) {
  .res <- .foceiSolveWithId(model, pars, fit$dataSav,
                            returnType = returnType,
                            atol = fit$control$atol[1], rtol = fit$control$rtol[1],
                            maxsteps = fit$control$maxstepsOde,
                            hmin = fit$control$hmin, hmax = fit$control$hmax, hini = fit$control$hini,
                            transitAbs = fit$control$transitAbs, maxordn = fit$control$maxordn,
                            maxords = fit$control$maxords, method = fit$control$method,
                            keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=addCov
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
                              keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=addCov)
    RxODE::rxSolveFree()
    if (any(is.na(.res$rx_pred_))) {
      .res <- .foceiSolveWithId(model, pars, fit$dataSav,
                                returnType = returnType,
                                atol = fit$control$atol[1], rtol = fit$control$rtol[1],
                                maxsteps = fit$control$maxstepsOde * 2,
                                hmin = fit$control$hmin, hmax = fit$control$hmax / 2, hini = fit$control$hini,
                                transitAbs = fit$control$transitAbs, maxordn = fit$control$maxordn,
                                maxords = fit$control$maxords, method = "dop853",
                                keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=addCov)
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
.foceiPredIpredList <- function(fit, data=fit$dataSav, thetaEtaParameters=.foceiThetaEtaParameters(fit), keep=NULL, predOnly=!is.null(fit$model$inner),
                                addDosing=FALSE, subsetNonmem=TRUE) {
  if (!predOnly & is.null(fit$model$inner)) {
    # Add inner problem calculation for cwres calculation
    assign("model", fit$uif$inner, envir=fit$env)
    if (is.null(fit$model$inner)) {
      stop("problem calculating focei's inner ODEs", call.=FALSE) # nocov
    }
  }
  .keep <- keep
  .names <- names(data)
  .lowerNames <- tolower(.names)
  for (.n in c("dv", "cens", "limit")) {
    .w <- which(.lowerNames == .n)
    if (length(.w) == 1L) .keep <- c(.keep, .names[.w])
  }
  .keep <- unique(.keep)
  .ipredModel <- fit$model$inner
  if (predOnly) .ipredModel <- fit$model$pred.only
  .ret <- list(ipred = .foceiSolvePars(fit, .ipredModel, thetaEtaParameters$ipred,
                                       returnType="data.frame.TBS", keep=.keep, what="ipred",
                                       addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=predOnly),
               pred = .foceiSolvePars(fit, .ipredModel, thetaEtaParameters$pred,returnType="data.frame", what="pred",
                                      addDosing=addDosing, subsetNonmem=subsetNonmem),
               etaLst=thetaEtaParameters$eta.lst)
  if (!predOnly) {
    .ret <- c(.ret, list(pred.only=.foceiSolvePars(fit, fit$model$pred.only, thetaEtaParameters$ipred,
                                   returnType="data.frame", keep=.keep, what="ebe",
                                   addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=TRUE)))
  }
  .ret
}

.getRelevantLhs <- function(fit, keep=NULL, ipred=NULL) {
  .ret <- setdiff(fit$model$pred.only$lhs,fit$ini$name)
  .w <- which(regexpr("^rx", .ret) == -1)
  .ret <- unique(c(.ret[.w], keep))
  if (any(.ret == "tad")) {
    if (all(is.na(ipred$tad))) {
      .ret <- setdiff(.ret, c("tad", "dosenum"))
    }
  }
  .ret
}


.calcCwres0 <- function(fit, data=fit$dataSav, thetaEtaParameters=.foceiThetaEtaParameters(fit),
                        table=tableControl(), dv=NULL, predOnly=FALSE,
                        addDosing=FALSE, subsetNonmem=TRUE, keep=NULL, npde=FALSE,
                        .prdLst) {
  if (!inherits(dv, "numeric")) {
    dv <- .prdLst$ipred$dv
    table$doSim <- TRUE
  } else {
    table$doSim <- FALSE
  }
  if (npde) {
    .sim <- .npdeSim(fit, nsim = table$nsim, ties = table$ties, seed = table$seed,
                     cholSEtol = table$cholSEtol, addDosing=addDosing, subsetNonmem=subsetNonmem, cores=table$cores)
    .Call(`_nlmixr_npdeCalc`, .sim, .prdLst$ipred$dv, .prdLst$ipred$evid,
          .prdLst$ipred$cens, .prdLst$ipred$limit, table)
  } else {
    if (predOnly){
      .state <- c(fit$model$ipred$state, fit$model$ipred$stateExtra)
      .lhs <- setdiff(unique(.getRelevantLhs(fit, keep, .prdLst$ipred)), .state)
      .params <- setdiff(intersect(names(fit$dataSav),fit$model$ipred$params),c("CMT","cmt","Cmt", .state, .lhs))
      .Call(`_nlmixr_resCalc`, .prdLst, fit$omega,
            fit$eta, .prdLst$ipred$dv, .prdLst$ipred$evid, .prdLst$ipred$cens,
            .prdLst$ipred$limit, .lhs, .state, .params, table)
    } else {
      .state <- c(fit$model$pred.only$state, fit$model$pred.only$stateExtra)
      .lhs <- setdiff(unique(.getRelevantLhs(fit, keep, .prdLst$pred.only)), .state)
      .params <- setdiff(intersect(names(fit$dataSav),fit$model$pred.only$params),c("CMT","cmt","Cmt", .state, .lhs))
      .Call(`_nlmixr_cwresCalc`, .prdLst, fit$omega,
            fit$eta, .prdLst$ipred$dv, .prdLst$ipred$evid, .prdLst$ipred$cens,
            .prdLst$ipred$limit, .lhs, .state, .params, table)
    }
  }
}

.calcCwres <- function(fit, data=fit$dataSav, thetaEtaParameters=.foceiThetaEtaParameters(fit),
                       table=tableControl(), dv=NULL, predOnly=FALSE,
                       addDosing=FALSE, subsetNonmem=TRUE, keep=NULL, npde=FALSE,
                       .prdLst=NULL) {
  if (!inherits(table, "tableControl")) table <- do.call(tableControl, table)
  if (is.null(.prdLst)) {
    .prdLst <- .foceiPredIpredList(fit, data=data, keep=keep, thetaEtaParameters=thetaEtaParameters, predOnly=predOnly,
                                   addDosing=addDosing, subsetNonmem=subsetNonmem)
  }
  ## Split out so that .prdLst can be shared between npde/cwres npde/res
  .ret <- .calcCwres0(fit, data, thetaEtaParameters, table, dv=dv, predOnly,
                      addDosing, subsetNonmem, keep, npde, .prdLst=.prdLst)
  .dups <- which(duplicated(names(.ret)))
  if (length(.dups) > 0) {
    warning("some duplicate columns were dropped", call.=FALSE)
    .ret <- .ret[, -.dups]
  }
  .ret
}

.calcRes <- function(..., predOnly=TRUE) {
  .calcCwres(..., predOnly=predOnly)
}

.calcNpde <- function(..., npde=TRUE, predOnly=TRUE) {
  .calcCwres(..., npde=npde, predOnly=predOnly)
}

.calcIres <- function(fit, data=fit$dataSav, table=tableControl(), dv=NULL,
                      addDosing=FALSE, subsetNonmem=TRUE, keep=NULL) {
  if (!inherits(table, "tableControl")) table <- do.call(tableControl, table)
  .keep <- keep
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
  .state <- c(.ipred$state, .ipred$stateExtra)
  .lhs <- setdiff(unique(.getRelevantLhs(fit, keep, .ipred)), .state)
  .params <- setdiff(intersect(names(fit$dataSav),.ipred$params),c("CMT","cmt","Cmt", .state, .lhs))
  .ret <- .Call(`_nlmixr_iresCalc`, .ipred, dv, .ipred$evid, .ipred$cens, .ipred$limit,
                .lhs, .state, .params, table)
  .dups <- which(duplicated(names(.ret)))
  if (length(.dups) > 0) {
    warning("some duplicate columns were dropped", call.=FALSE)
    .ret <- .ret[, -.dups]
  }
  .ret
}

.calcShrinkOnly <- function(fit, thetaEtaParameters=.foceiThetaEtaParameters(fit)) {
  .omega <- fit$omega
  .ret <- .Call(`_nlmixr_calcShrinkOnly`, .omega, thetaEtaParameters$eta.lst, length(fit$eta[,1]))
  .ret[, -dim(.omega)[1] - 1]
}

.calcTables <- function(fit, data=fit$dataSav, thetaEtaParameters=.foceiThetaEtaParameters(fit),
                        table=tableControl(), keep=NULL) {
  if (!inherits(table, "tableControl")) table <- do.call(tableControl, table)
  if (is.null(table$cwres)) {
    table$cwres <- !is.null(fit$model$inner)
  }
  if (is.null(table$npde)) {
    table$npde <- FALSE
  }
  .predOnly <- !table$cwres
  .censMethod <- table$censMethod
  .ret <- vector("list",2)
  .thetaEtaParameters <- .foceiThetaEtaParameters(fit)
  .prdLst <- .foceiPredIpredList(fit, data=fit$dataSav, keep=keep, thetaEtaParameters=.thetaEtaParameters, predOnly=.predOnly,
                                 addDosing=table$addDosing, subsetNonmem=table$subsetNonmem)
  if (.censMethod %in% c(2L, 6L)) {
    if (!table$npde) {
      warning("censoring method requires npde, adding npde", call.=FALSE)
      table$npde <- TRUE
    }
    .npde1 <- TRUE
    .npde2 <- FALSE
  } else {
    .npde1 <- FALSE
    .npde2 <- TRUE
  }
  if ((.npde1 & table$npde) | !.npde1)
    .ret[[1]] <- .calcCwres(fit, data=fit$dataSav, thetaEtaParameters=.thetaEtaParameters, table=table,
                            predOnly=.predOnly, addDosing=table$addDosing, subsetNonmem=table$subsetNonmem,
                            keep=keep, .prdLst=.prdLst, npde=.npde1)
  if ((.npde2 & table$npde) | !.npde2)
    .ret[[2]] <- .calcCwres(fit, data=fit$dataSav, thetaEtaParameters=.thetaEtaParameters, table=table, dv=.ret[[1]][[1]],
                            predOnly=.predOnly, addDosing=table$addDosing, subsetNonmem=table$subsetNonmem,
                            keep=keep, .prdLst=.prdLst, npde=.npde2)
  .ret <- .Call(`_nlmixr_popResFinal`, .ret)
  .dups <- which(duplicated(names(.ret)))
  if (length(.dups) > 0) {
    warning("some duplicate columns were dropped", call.=FALSE)
    .ret <- .ret[, -.dups]
  }
  .ret
}

.cloneEnv <- function(env) {
  .env <- new.env(parent = emptyenv())
  for (.x in ls(env, all.names=TRUE)) {
    assign(.x, get(.x, env), .env)
  }
  return(.env)
}

#' NPDE calculation for nlmixr
#'
#' @param object nlmixr fit object
#' @param updateObject Boolean indicating if original object should be updated.  By default this is TRUE.
#' @param table `tableControl()` list of options
#' @inheritParams foceiControl
#' @inheritParams addCwres
#' @param ... Other ignored parameters.
#'
#'
#' @return New nlmixr fit object
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'   ini({
#'     ## You may label each parameter with a comment
#'     tka <- 0.45 # Log Ka
#'     tcl <- log(c(0, 2.7, 100)) # Log Cl
#'     ## This works with interactive models
#'     ## You may also label the preceding line with label("label text")
#'     tv <- 3.45; label("log V")
#'     ## the label("Label name") works with all models
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' f <- nlmixr(one.cmt, theo_sd, "saem")
#'
#' # even though you may have forgotten to add the NPDE, you can add it to the data.frame:
#'
#' f <- addNpde(f)
#'
#' }
##' @export
addNpde <- function(object, updateObject = TRUE,
                    table = tableControl(), ...,
                    envir=parent.frame(1)) {
  .pt <- proc.time()
  .objName <- substitute(object)
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  if (any(names(object) == "NPDE")) {
    warning("already contains NPDE")
    return(object)
  }
  message("add npde")
  table$npde <- TRUE
  .npde <- .calcNpde(object, dv=object$DV)
  .cls <- class(object)
  .new <- cbind(object, .npde[[2]])
  class(.new) <- .cls
  .env <- .new$env
  if (inherits(updateObject, "logical")) {
    if (updateObject) {
      .parent <- envir
      .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
        if (.cur == .objName && identical(.parent[[.cur]]$env, .env)) {
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
  }
  .env$time <- .data.frame(.env$time, npde = (proc.time() - .pt)["elapsed"], check.names = FALSE)
  message("done")
  .new
}

#' Add table information to nlmixr fit object without tables
#'
#' @param object nlmixr family of objects
#' @param updateObject Update the object (default FALSE)
#' @param data Saved data from
#' @param thetaEtaParameters Intenral theta/eta parameters
#' @param table a `tableControl()` list of options
#' @param keep Character Vector of items to keep
#' @param drop Character Vector of items to drop or NULL
#' @param envir ENvironment to search for updating
#' @return Fit with table information attached
#' @author Matthew Fidler
#' @export
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'   ini({
#'     ## You may label each parameter with a comment
#'     tka <- 0.45 # Log Ka
#'     tcl <- log(c(0, 2.7, 100)) # Log Cl
#'     ## This works with interactive models
#'     ## You may also label the preceding line with label("label text")
#'     tv <- 3.45; label("log V")
#'     ## the label("Label name") works with all models
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' # run without tables step
#' f <- nlmixr(one.cmt, theo_sd, "saem", control=list(calcTables=FALSE))
#'
#' print(f)
#'
#' # Now add the tables
#'
#' f <- addTable(f)
#'
#' print(f)
#'
#' }

addTable <- function(object, updateObject = FALSE, data=object$dataSav, thetaEtaParameters=.foceiThetaEtaParameters(object),
                     table=tableControl(), keep=NULL, drop=NULL,
                     envir = parent.frame(1)) {
  .pt <- proc.time()
  message("Calculating residuals/tables")
  .objName <- substitute(object)
  if (!inherits(object, "nlmixrFitCore")) {
    stop("requires a nlmixr fit object")
  }
  .fit <- object$env
  if (exists("origControl", .fit)) {
    .control <- .fit$origControl
  } else if (exists("control", .fit)) {
    .control <- .fit$control
  } else {
    .control <- foceiControl()
  }
  if (is.null(.fit$omega)) {
    .df <- .calcIres(.fit, data=data, table=table, dv=NULL,
                     addDosing=table$addDosing, subsetNonmem=table$subsetNonmem, keep=keep)
  } else {
    .tabs <- .calcTables(.fit, data=data, table=table, keep=keep)
    assign("shrink", .tabs$shrink, .fit)
    .df <- .tabs$resid
  }
  drop <- c(drop, "rxLambda", "rxYj")
  .w <- -which(names(.df) %in% drop)
  if (length(.w) > 0) .df <- .df[, .w, drop=FALSE]
  .isDplyr <- requireNamespace("dplyr", quietly = TRUE)
  if (!.isDplyr) {
    .isDataTable <- requireNamespace("data.table", quietly = TRUE)
    if (.isDataTable) {
      .df <- data.table::data.table(.df)
    }
  } else {
    .df <- tibble::as_tibble(.df)
  }
  .cls <- class(.df)
  if (!any(names(.control) == "interaction")) {
    .control$interaction <- FALSE
  }
  if (.control$interaction) {
    .cls <- c(paste0("nlmixr", .fit$method, "i"), "nlmixrFitData", "nlmixrFitCore", .cls)
  } else {
    .cls <- c(paste0("nlmixr", .fit$method), "nlmixrFitData", "nlmixrFitCore", .cls)
  }
  if (inherits(updateObject, "logical")) {
    if (!updateObject){
      .fit <- .cloneEnv(.fit)
    }
  }
  class(.fit) <- "nlmixrFitCoreSilent"
  attr(.cls, ".foceiEnv") <- .fit
  class(.df) <- .cls
  if (inherits(updateObject, "logical")) {
    if (updateObject) {
      .parent <- envir
      .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
        if (.cur == .objName && identical(.parent[[.cur]]$env, .fit$env)) {
          return(.cur)
        }
        return(NULL)
      }))
      if (length(.bound) == 1) {
        if (exists(.bound, envir = .parent)) {
          assign(.bound, .df, envir = .parent)
        }
      }
    }
  }
  message("done")
  .fit$time <- .data.frame(.fit$time, table = (proc.time() - .pt)["elapsed"], check.names = FALSE)
  .df
}

##' Output table/data.frame options
##'
##' @param npde When TRUE, request npde regardless of the algorithm used.
##'
##' @param cwres When TRUE, request CWRES and FOCEi likelihood
##'     regardless of the algorithm used.
##'
##' @param censMethod Handle censoring method:
##'
##'  - `"truncated-normal"` Simulates from a truncated normal distribution under the assumption of the model and censoring.
##'
##'  - `"cdf"` Use the cdf-method for censoring with npde and use this for any other residuals (`cwres` etc)
##'
##'  - `"omit"` omit the residuals for censoring
##'
##' @param ties When `TRUE` jitter prediction-discrepancy points to discourage ties in cdf.
##'
##' @param cholSEtol The tolerance for the `RxODE::choleSE` function
##'
##' @param eta is a Boolean indicating if `eta` values will be included (default `TRUE`)
##'
##' @param state is a Boolean indicating if `state` values will be included (default `TRUE`)
##'
##' @param lhs is a Boolean indicating if remaining `lhs` values will be included (default `TRUE`)
##'
##' @param covariates is a Boolean indicating if covariates will be included (default `TRUE`)
##'
##' @inheritParams addNpde
##' @inheritParams RxODE::rxSolve
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
                         nsim = 300, ties = TRUE,
                         censMethod=c("truncated-normal", "cdf", "ipred", "pred", "epred", "omit"),
                         seed = 1009,
                         cholSEtol=(.Machine$double.eps)^(1/3),
                         state=TRUE,
                         lhs=TRUE,
                         eta=TRUE,
                         covariates=TRUE,
                         addDosing=FALSE, subsetNonmem = TRUE,
                         cores=NULL) {
  if (is.null(cores)) {
    cores = RxODE::rxCores()
  }
  .ret <- list(
    npde = npde, cwres = cwres, nsim = nsim, ties = ties, seed = seed,
    censMethod=setNames(c("truncated-normal"=3L, "cdf"=2L, "omit"=1L, "pred"=5L, "ipred"=4L, "epred"=6L)[match.arg(censMethod)], NULL),
    cholSEtol=cholSEtol, state=state, lhs=lhs, eta=eta, covariates=covariates, addDosing=addDosing, subsetNonmem=subsetNonmem, cores=cores)
  class(.ret) <- "tableControl"
  return(.ret)
}
