.setCov <- function(obj, ...) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  .pt <- proc.time()
  .env <- obj
  if (RxODE::rxIs(obj, "nlmixrFitData")) {
    .env <- obj$env
  }
  if (exists("cov", .env)) {
    .cur <- list(.env$cov)
    names(.cur) <- .env$covMethod
    if (exists("covList", .env)) {
      if (is.null(.env$covList[[.env$covMethod]])) {
        .covList <- c(.env$covList, .cur)
      } else {
        .covList <- .env$covList
      }
    } else {
      .covList <- .cur
    }
    assign("covList", .covList, .env)
  }
  .control <- .env$control
  .control$maxInnerIterations <- 0L;
  .control$maxOuterIterations <- 0L
  .control$boundTol <- 0 # turn off boundary
  .control$calcTables <- FALSE
  .lst <- list(...)
  .env2 <- new.env(parent = emptyenv())
  for (.n in names(.lst)) {
    .control[[.n]] <- .lst[[.n]]
  }
  if (!is.null(.lst$covMethod)) {
    if (RxODE::rxIs(.lst$covMethod, "character")) {
      .lst$covMethod <- match.arg(.lst$covMethod, c("r,s", "r", "s"))
      .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
      .control$covMethod <- .covMethodIdx[.lst$covMethod]
    } else if (inherits(.lst$covMethod, "matrix")){
      .env2$cov <- as.matrix(.lst$covMethod)
      .env2$uif <- obj$uif
      .control$covMethod <- 0L
    } else if (length(.lst$covMethod) == 1) {
      if (.lst$covMethod == "") {
        .control$covMethod <- 0L
      }
    }
  } else if (.control$covMethod == 0L) {
    .control$covMethod <- 1L
  }
  ## covDerivMethod=c("central", "forward"),
  if (!is.null(.lst$hessEps)) {
    .control$hessEps <- .lst$hessEps
    .lst$hessEps <- NULL
  }
  if (!is.null(.lst$gillKcov)) {
    .control$gillKcov <- .lst$gillKcov
    .lst$gillKcov <- NULL
  }
  if (!is.null(.lst$gillStepCov)) {
    .control$gillStepCov <- .lst$gillStepCov
    .lst$gillStepCov <- NULL
  }
  if (!is.null(.lst$gillFtolCov)) {
    .control$gillFtolCov <- .lst$gillFtolCov
    .lst$gillFtolCov <- NULL
  }
  if (!is.null(.lst$rmatNorm)) {
    .control$rmatNorm <- .lst$rmatNorm
    .lst$rmatNorm <- NULL
  }
  if (!is.null(.lst$smatNorm)) {
    .control$smatNorm <- .lst$smatNorm
    .lst$smatNorm <- NULL
  }
  if (!is.null(.lst$covGillF)) {
    .control$covGillF <- .lst$covGillF
    .lst$covGillF <- NULL
  }
  if (!is.null(.lst$covSmall)) {
    .control$covSmall <- .lst$covSmall
    .lst$covSmall <- NULL
  }
  .dat <- getData(obj)
  .uif <- obj$uif
  .mat <- as.matrix(nlme::random.effects(obj)[, -1])
  .skipCov <- obj$skipCov
  .inits <- list(
    THTA = as.vector(nlme::fixed.effects(obj)),
    OMGA = focei.eta.nlmixrFitCore(obj)
  )
  .fit2 <- foceiFit.data.frame0(
    data = .dat,
    inits = .inits,
    PKpars = .uif$theta.pars,
    ## par_trans=fun,
    model = .uif$rxode.pred,
    pred = function() {
      return(nlmixr_pred)
    },
    err = .uif$error,
    lower = .uif$focei.lower,
    upper = .uif$focei.upper,
    thetaNames = .uif$focei.names,
    etaNames = .uif$eta.names,
    etaMat = .mat,
    skipCov = .skipCov,
    control = .control,
    env=.env2
  )
  .env$cov <- .fit2$cov
  .env$popDf <- .fit2$popDf
  .env$popDfSig <- .fit2$popDfSig
  .env$covMethod <- .fit2$covMethod
  .updateParFixed(.env)
  .parent <- parent.frame(2)
  .bound <- do.call("c", lapply(ls(.parent), function(.cur) {
    if (identical(.parent[[.cur]], obj)) {
      return(.cur)
    }
    return(NULL)
  }))
  message(paste0("Updated original fit object ", ifelse(is.null(.bound), "", crayon::yellow(.bound))))
  .env$time$covariance <- (proc.time() - .pt)["elapsed"]
  return(.env$cov)
}

##' Set the covariance type based on prior calculated covariances
##'
##' @param fit nlmixr fit
##'
##' @param method covariance method
##'
##' @return Fit object with covariance updated
##'
##' @author Matt Fidler
##'
##' @export
setCov <- function(fit, method) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  .pt <- proc.time()
  .env <- fit
  if (RxODE::rxIs(fit, "nlmixrFitData")) {
    .env <- fit$env
  }
  if (method == .env$covMethod) {
    stop("no need to switch covariance methods, already set to '",
         method,
         "'", call.=FALSE)
  }
  if (exists("covList", .env)) {
    .covList <- .env$covList
    .cov <- .covList[[method]]
    if (!is.null(.cov)) {
      .setCov(fit, covMethod=.cov)
      .env$covMethod <- method
      return(fit)
    }
  }
  stop("different covariance types have not been calculated",
         call.=FALSE)
}

##' @export
getVarCov.nlmixrFitCore <- function(obj, ...) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  .env <- obj
  if (RxODE::rxIs(obj, "nlmixrFitData")) {
    .env <- obj$env
  }
  .force <- FALSE
  .args <- list(...)
  if (!is.null(.args$force)) {
    .force <- .args$force
  }
  if (exists("cov", envir = .env) && !.force) {
    if (RxODE::rxIs(.env$cov, "matrix")) {
      return(.env$cov)
    }
  }
  .setCov(obj, ...)
}

##' @export
getVarCov.nlmixrFitCoreSilent <- getVarCov.nlmixrFitCore
