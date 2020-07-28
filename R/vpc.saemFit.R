multi2 <- function(mu, vmat, n) {
  eta <- matrix(rnorm(length(mu) * n), ncol = n, nrow = length(mu))
  Q <- chol(vmat, pivot = TRUE)
  pivot <- attr(Q, "pivot")
  oo <- order(pivot)
  para <- t(Q[, oo]) %*% eta
  sweep(para, 1, mu, "+")
}

rmvnorm <- function(n, mu, vmat) multi2(mu, vmat, n)

##' VPC for nlmixr saemFit objects
##'
##' @param fit  saemFit object
##' @param dat  Data to augment the saemFit vpc simulation
##' @param nsim Number of simulations for the VPC
##' @param by Variables to condition the VPC
##' @inheritParams vpc::vpc
##' @param ... Other arguments sent to \code{\link[vpc:vpc]{vpc_vpc}}
##' @return vpc object from the \code{\link[vpc:vpc]{vpc_vpc}} package
##' @author Wenping Wang
##' @export
vpc_saemFit <- function(fit, dat, nsim = 100, by = NULL, ...) {
  if (class(fit) == "nlmixr.ui.saem") fit <- as.saem(fit)
  .env <- attr(fit, "env")
  saem.cfg <- attr(fit, "saem.cfg")
  if (.env$is.ode) {
    .evtM <- saem.cfg$evtM
    .rx <- .env$model
    .pars <- .rx$params
    .pars <- setNames(rep(1.1, length(.pars)), .pars)
    suppressWarnings(do.call(
      RxODE::rxSolve,
      c(
        list(
          object = .rx, params = .pars,
          events = .evtM, .setupOnly = 2L
        ),
        saem.cfg$optM
      )
    ))
    RxODE::rxLock(.rx)
    RxODE::rxAllowUnload(FALSE)
    on.exit({
      RxODE::rxUnlock(.rx)
      RxODE::rxAllowUnload(TRUE)
    })
  }
  dyn.load(.env$saem.dll)
  assignInMyNamespace(".protectSaemDll", .env$saem.dll)
  on.exit(
    {
      assignInMyNamespace(".protectSaemDll", "")
    },
    add = TRUE
  )
  dopred <- attr(fit, "dopred")
  resMat <- fit$resMat
  ares <- resMat[, 1]
  bres <- resMat[, 2]
  ntotal <- saem.cfg$ntotal
  ix_endpnt <- saem.cfg$ix_endpnt[1:ntotal] + 1
  ares <- ares[ix_endpnt]
  bres <- bres[ix_endpnt]

  if (!is.null(by)) {
    if (by %in% names(dat)) {
      dat$grp <- eval(parse(text = paste0("dat$", by)))
    } else {
      msg <- paste0(by, " not found in data")
      stop(msg)
    }
  }
  else {
    dat$grp <- T
  }

  xd <- subset(dat, EVID == 0)
  nsub <- length(unique(xd$ID))
  ntim <- dim(xd)[1]
  ord <- rep(1:ntim, nsim)
  sim <- rep(1:nsim, each = ntim)

  mpost_rand <- t(matrix(fit$Plambda, length(fit$Plambda), nsub))
  s <- lapply(1:nsim, function(k) {
    mpost_rand1 <- t(rmvnorm(nsub, fit$Plambda[saem.cfg$i1 + 1], fit$Gamma2_phi1))
    mpost_rand[, saem.cfg$i1 + 1] <- mpost_rand1
    p <- dopred(mpost_rand, saem.cfg$evt, saem.cfg$opt)
    dim(p) <- NULL
    res <- (ares + p * bres) * rnorm(ntim, 0, 1)
    p + res
  })
  xs <- do.call("cbind", s)
  df <- cbind(xd[ord, c("ID", "TIME", "grp")], DV = as.vector(xs), SIM = sim)

  ns <- loadNamespace("vpc")
  if (exists("vpc_vpc", ns)) {
    vpcn <- "vpc_vpc"
  } else {
    vpcn <- "vpc"
  }
  call <- as.list(match.call(expand.dots = TRUE))[-1]
  if (!is.null(by)) {
    call$stratify <- c("grp")
    call$facet <- "wrap"
  }
  call <- call[names(call) %in% methods::formalArgs(getFromNamespace(vpcn, "vpc"))]
  p <- do.call(getFromNamespace(vpcn, "vpc"), c(list(sim = df, obs = dat), call), envir = parent.frame(1))
  print(p)

  invisible(df)
}

#' @rdname vpc_saemFit
#' @export
vpc.saemFit <- function(sim, ...) {
  vpc_saemFit(sim, ...)
}
