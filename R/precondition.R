#' Linearly re-parameterize the model to be less sensitive to rounding errors
#'
#' @param fit A nlmixr fit to be preconditioned
#' @param estType Once the fit has been linearly reparametrized,
#'   should a "full" estimation, "posthoc" estimation or simply a
#'   estimation of the covariance matrix "none" before the fit is
#'   updated
#' @param ntry number of tries before giving up on a pre-conditioned
#'   covariance estimate
#'
#'@export
preconditionFit <- function(fit, estType = c("full", "posthoc", "none"),
                            ntry=10L) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  if (!exists("R", fit$env)) {
    stop("this assumes a covariance matrix with a R matrix",
         call.=FALSE)
  }
  .R <- fit$R
  .covMethod <- ""
  .i <- 1
  while(.i < ntry & .covMethod != "r,s"){
    .i <- .i + 1
    pre <- preCondInv(.R)
    P <- symengine::Matrix(pre)
    d0 <- dimnames(fit$R)[[1]]
    d <- paste0("nlmixrPre_", dimnames(fit$R)[[1]])
    d2 <- symengine::Matrix(d)
    modExtra <- paste(paste0(d0, "=", sapply(as.vector(P %*% d2), as.character)), collapse = "\n")
    preInv <- solve(pre)

    .ini <- as.data.frame(fit$uif$ini)
    newEst <- setNames(as.vector(preInv %*% matrix(fit$theta[d0])), d0)
    for (v in d0) {
      .w <- which(.ini$name == v)
      .ini$lower[.w] <- -Inf
      .ini$upper[.w] <- Inf
      .ini$est[.w] <- newEst[v]
      .ini$name[.w] <- paste0("nlmixrPre_", v)
    }
    .newModel <- eval(parse(text = paste0("function(){", modExtra, "\n", fit$uif$fun.txt, "}"), keep.source = TRUE))
    class(ini) <- c("nlmixrBounds", "data.frame")
    newModel <- .finalizeUiModel(
      nlmixrUIModel(.newModel, .ini, NULL),
      as.list(fit$uif$meta)
    )
    .ctl <- fit$origControl
    estType <- match.arg(estType)
    if (estType == "none") {
      .ctl$maxInnerIterations <- 0
      .ctl$maxOuterIterations <- 0
      .ctl$boundTol <- 0
      .ctl$etaMat <- as.matrix(fit$eta[, -1])
      .ctl$calcTables <- FALSE
    } else if (estType == "posthoc") {
      .ctl$maxOuterIterations <- 0
      .ctl$boundTol <- 0
      .ctl$calcTables <- FALSE
    } else if (estType == "full") {
      .ctl$boundTol <- 0
      .ctl$calcTables <- FALSE
    }
    ## FIXME compare objective functions
    newFit <- suppressWarnings(nlmixr(newModel, getData(fit), est = "focei", control = .ctl))
    .R <- newFit$R
    .covMethod <- newFit$covMethod
  }
  if (.covMethod != "r,s") {
    stop("preconditioning failed after ", ntry, "tries",
         call.=FALSE)
  }
  cov <- pre %*% newFit$cov %*% t(pre)
  dimnames(cov) <- dimnames(pre)
  assign("precondition", cov, env=fit$env)
  .setCov(fit, covMethod=cov)
  assign("covMethod", "precondition", fit$env)
  return(invisible(fit$env$precondition))
}
