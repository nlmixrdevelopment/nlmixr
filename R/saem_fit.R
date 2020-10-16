## saem_fit.R: population PK/PD modeling library
##
## Copyright (C) 2014 - 2016  Wenping Wang
##
## This file is part of nlmixr.
##
## nlmixr is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## nlmixr is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with nlmixr.  If not, see <http://www.gnu.org/licenses/>.

# genSaemUserFunction(f$rxode.pred, f$saem.pars, f$pred, f$error)
genSaemUserFunction <- function(model, PKpars = attr(model, "default.pars"), pred = NULL, err=NULL,
                                control=saemControl(), inPars=NULL) {
  .x <- deparse(body(pred))
  .len <- length(.x)
  .x <- if(.x[1]=="{") .x[2:(.len-1)] else .x
  .len <- length(.x)
  .nendpnt <- .len
  .mod <- RxODE::RxODE(RxODE::rxGenSaem(model, function() { return(nlmixr_pred)}, PKpars,
                                        sum.prod=control$sum.prod,
                                        optExpression=control$optExpression))
  .fnPred <- bquote(function(a, b, c){
    RxODE::rxLoad(.(.mod))
    RxODE::rxLock(.(.mod))
    RxODE::rxAllowUnload(FALSE)
    on.exit({RxODE::rxUnlock(.(.mod)); RxODE::rxAllowUnload(TRUE); RxODE::rxSolveFree()})
    .Call(`_nlmixr_saem_do_pred`, a, b, c);
  })
  .fn <- bquote(function(a, b, c){
    RxODE::rxLoad(.(.mod))
    RxODE::rxLock(.(.mod))
    on.exit({RxODE::rxUnlock(.(.mod)); RxODE::rxAllowUnload(TRUE); RxODE::rxSolveFree()})
    if (missing(b) && missing(c)){
      .ret <- .Call(`_nlmixr_saem_fit`, a, PACKAGE="nlmixr")
      attr(.ret, "dopred") <- .(.fnPred)
      return(.ret)
    } else {
      .curFn <- .(.fnPred)
      return(.curFn(a, b, c))
    }
  })
  .param <- RxODE::rxParam(.mod)
  .inits <- names(RxODE::rxInits(.mod))
  .nrhs <- length(.param) - length(.inits)
  if (any(.param == "CMT")){
    inPars <- unique(c(inPars, "CMT"))
  }
  .parmUpdate <- sapply(.param, function(x) {
    if (any(x == inPars)) {
      return(0L)
    } else if (any(x == .inits)) {
      return(0L)
    } else {
      return(1L)
    }
  })
  .fn <- eval(.fn)
  attr(.fn, "form") <- "ode" ## Not sure this is necessary any more
  attr(.fn, "neq") <- length(RxODE::rxState(.mod))
  attr(.fn, "nlhs") <- length(RxODE::rxLhs(.mod))
  attr(.fn, "nrhs") <- sum(.parmUpdate)
  attr(.fn, "paramUpdate") <- .parmUpdate
  attr(.fn, "rx") <- .mod
  attr(.fn, "inPars") <- inPars
  attr(.fn, "nendpnt") <- .nendpnt # not calculated; Is this a problem?
  return(.fn)
}

#' Generate an SAEM model
#'
#' Generate an SAEM model using either closed-form solutions or ODE-based model definitions
#'
#' @param model a compiled SAEM model by gen_saem_user_fn()
#' @param PKpars PKpars function
#' @param pred pred function;  This will be a focei-style pred
#' @param inPars a character vector of parameters to be read from the
#'   input dataset (including time varying covariates)
#' @details Fit a generalized nonlinear mixed-effect model using the
#'   Stochastic Approximation Expectation-Maximization (SAEM)
#'   algorithm
#'
#' @author Matthew Fidler & Wenping Wang
#' @keywords internal
#' @export
gen_saem_user_fn <- genSaemUserFunction

#' Configure an SAEM model
#'
#' Configure an SAEM model by generating an input list to the SAEM model function
#'
#' @param model a compiled saem model by gen_saem_user_fn()
#' @param data input data
#' @param inits initial values
#' @param mcmc a list of various mcmc options
#' @param ODEopt optional ODE solving options
#' @param seed seed for random number generator
#' @param distribution one of c("normal","poisson","binomial")
#' @param fixed a character vector of fixed effect only parameters (no random effects attached) to be fixed
#' @param DEBUG Integer determining if debugging is enabled.
#' @inheritParams RxODE::rxSEinner
#' @inheritParams saemControl
#' @details
#'    Fit a generalized nonlinear mixed-effect model by he Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Wenping Wang & Matthew Fidler
#' @examples
#' \dontrun{
#' library(nlmixr)
#'
#' ode <- "d/dt(depot) =-KA*depot;
#'         d/dt(centr) = KA*depot - KE*centr;"
#' m1 = RxODE(ode, modName="m1")
#' # ode <- "C2 = centr/V;
#' #      d/dt(depot) =-KA*depot;
#' #      d/dt(centr) = KA*depot - KE*centr;"
#' # m2 = RxODE(ode, modName="m2")
#'
#' # Note: only use the '=' assignment, not the '<-' at this point
#'
#' PKpars <- function() {
#'   CL <- exp(lCL)
#'   V <- exp(lV)
#'   KA <- exp(lKA)
#'   KE <- CL / V
#'   # initCondition = c(0, KE - CL/V)
#' }
#' PRED <- function() centr / V
#' PRED2 <- function() C2
#'
#'  saem_fit <- gen_saem_user_fn(model=m1, PKpars, pred=PRED)
#' # saem_fit <- gen_saem_user_fn(model=m2, PKpars, pred=PRED2)
#'
#'
#' #--- saem cfg
#' nmdat <- theo_sd
#' model <- list(saem_mod = saem_fit, covars = "WT")
#' inits <- list(theta = c(.05, .5, 2))
#' cfg <- configsaem(model, nmdat, inits)
#' cfg$print <- 50
#'
#' # cfg$Rfn = nlmixr:::Ruser_function_cmt
#' # dyn.load("m1.d/m1.so");cfg$Rfn = nlmixr:::Ruser_function_ode
#' fit <- saem_fit(cfg)
#' df <- simple.gof(fit)
#' xyplot(DV ~ TIME | ID, df, type = c("p", "l"), lwd = c(NA, 1), pch = c(1, NA), groups = grp)
#' fit
#' }
#' @export
configsaem <- function(model, data, inits,
                       mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                       ODEopt = list(atol = 1e-6, rtol = 1e-4, method = "lsoda", transitAbs = FALSE, maxeval = 100000),
                       distribution = c("normal", "poisson", "binomial"),
                       addProp=c("combined2","combined1"),
                       seed = 99, fixed = NULL, DEBUG = 0,
                       normal=c("rnorm", "vandercorput"),
                       tol=1e-4, itmax=100L, type=c("nelder–mead", "newuoa"),
                       lambdaRange=3, powRange=10) {
  type.idx <- c("nelder–mead" = 1L, "newuoa" = 2L)
  type <- match.arg(type)
  type <- type.idx[type]
  normal <- match.arg(normal)
  names(ODEopt) <- gsub("transit_abs", "transitAbs", names(ODEopt))
  ODEopt <- do.call(RxODE::rxControl, ODEopt)
  # mcmc=list(niter=c(200,300), nmc=3, nu=c(2,2,2));ODEopt = list(atol=1e-6, rtol=1e-4, stiff=1, transit_abs=0);distribution=c("normal","poisson","binomial");seed=99;data=dat;distribution=1;fixed=NULL
  set.seed(seed)
  distribution.idx <- c("normal" = 1, "poisson" = 2, "binomial" = 3, "lnorm" = 4)
  distribution <- match.arg(distribution)
  distribution <- distribution.idx[distribution]
  .data <- data
  ## RxODE::rxTrans(data, model)
  data <- list(nmdat = data)

  neq    <- attr(model$saem_mod, "neq")
  nlhs   <- attr(model$saem_mod, "nlhs")
  inPars <- attr(model$saem_mod, "inPars")
  ninputpars <- length(inPars)
  opt <- optM <- c(list(neq = neq, nlhs = nlhs, inits = numeric(neq)),
    ninputpars = ninputpars, inPars = inPars
  )

  model$N.eta <- attr(model$saem_mod, "nrhs")
  model$nendpnt <- attr(model$saem_mod, "nendpnt")
  if (is.null(model$nendpnt)) model$nendpnt <- 1

  if (is.null(model$log.eta)) model$log.eta <- rep(TRUE, model$N.eta)
  if (is.null(model$omega)) model$omega <- diag(model$N.eta)
  if (is.null(model$res.mod)) model$res.mod <- rep(1, model$nendpnt)
  if (is.null(inits$omega)) inits$omega <- rep(1, model$N.eta) * 4
  if (is.null(inits$ares)) inits$ares <- 10
  if (is.null(inits$bres)) inits$bres <- 1
  if (is.null(inits$cres)) inits$cres <- 1
  if (is.null(inits$lres)) inits$lres <- 1
  if (is.null(mcmc$print)) mcmc$print <- 1
  if (is.null(names(inits$theta))) names(inits$theta) <- rep("", length(inits$theta))
  inits.save <- inits
  inits$theta.fix <- matrix(names(inits$theta),
                            byrow = T,
                            ncol = model$N.eta)
  inits$theta <- matrix(inits$theta, byrow = T, ncol = model$N.eta)
  model$cov.mod <- 1 - is.na(inits$theta)
  data$N.covar <- nrow(inits$theta) - 1
  inits$theta[is.na(inits$theta)] <- 0


  ###  FIXME
  mcmc$stepsize <- 0:1
  mcmc$burn.in <- 300

  ###  FIXME: chk covars as char vec
  wh <- setdiff(c(model$covars, inPars), names(data$nmdat))
  if (length(wh)) {
    msg <- paste0("covariate(s) not found: ", paste(wh, collapse = ", "))
    stop(msg)
  }
  s <- subset(data$nmdat, EVID == 0)
  data$data <- as.matrix(s[, c("ID", "TIME", "DV", c(model$covars, inPars))])

  ###  chk for no obs records
  wh <- setdiff(unique(data$nmdat$ID), unique(data$data[, "ID"]))
  if (length(wh)) {
    msg <- paste0("No data with ID: ", paste(wh, collapse = ", "))
    stop(msg)
  }

  nphi <- model$N.eta
  mcov <- model$cov.mod
  covstruct <- model$omega

  check <- sum((covstruct - t(covstruct)) != 0)
  if (check) stop("illegal covstruct")
  check <- nphi - dim(covstruct)[1]
  if (check) stop("nphi and covstruct dim mismatch")

  check <- prod(mcov[1, ])
  if (check == 0) stop("structural par(s) absent")
  check <- nphi - dim(mcov)[2]
  if (check) stop("nphi and ncol(mcov) mismatch")
  check <- sum(dim(inits$theta) - dim(mcov) != 0)
  if (check) stop("initial theta's and mcov dim mismatch")
  check <- data$N.covar + 1 - dim(mcov)[1]
  if (check) stop("dim mcov and N.covar mismatch")

  check <- length(model$log.eta) - nphi
  if (check) stop("jlog length and nphi mismatch")

  check <- length(inits$omega) - nphi
  if (check) stop("length of omega inits and nphi mismatch")

  # check = mcmc$burn.in>sum(mcmc$niter)
  # if (check) stop("#burn-in exceeds niter")

  check <- prod(is.element(covstruct, c(0, 1)))
  if (check == 0) warning("non-zero value(s) in covstruct set to 1")
  covstruct[covstruct != 0] <- 1

  check <- prod(is.element(mcov, c(0, 1)))
  if (check == 0) warning("non-zero value(s) in mcov set to 1")
  mcov[mcov != 0] <- 1

  check <- sum(inits$theta[1, model$log.eta] <= 0)
  if (check) stop("illegal initial theta's")
  check <- sum(inits$omega <= 0)
  if (check) stop("illegal initial omega")
  # check = inits$sigma2<=0
  # if (check) stop("illegal initial sigma2")
  check <- sum(diag(covstruct) == 1)
  if (!check) stop("0 ETA's")
  y <- data$data[, "DV"]
  id <- data$data[, "ID"]
  check <- any(diff(unique(id)) != 1)
  if (check) stop("saem classic UI needs sequential ID. check your data")
  ntotal <- length(id)
  N <- length(unique(id))
  covariables <- if (is.null(model$covars)) NULL else unlist(stats::aggregate(.as.data.frame(data$data[, model$covars, drop = FALSE]), list(id), unique)[, -1, drop = FALSE])
  if (!is.null(covariables)) dim(covariables) <- c(N, data$N.covar)
  nb_measures <- table(id)
  ncov <- data$N.covar + 1
  nmc <- mcmc$nmc
  nM <- mcmc$nmc * N
  yM <- rep(y, nmc)
  mlen <- max(nb_measures)
  io <- t(sapply(nb_measures, function(x) rep(1:0, c(x, mlen - x))))
  ix <- rep(1:dim(io)[1], nmc)
  ioM <- io[ix, ]
  indioM <- grep(1, t(ioM)) - 1
  ## mPars <- if (ninputpars == 0) NULL else unlist(stats::aggregate(.as.data.frame(data$data[, inPars]), list(id), unique)[, -1])
  ## if (!is.null(mPars)) {
  ##   dim(mPars) <- c(N, ninputpars)
  ##   opt$mPars <- mPars
  ##   ix <- rep(1:dim(mPars)[1], nmc)
  ##   optM$mPars <- mPars[ix, ]
  ##   dim(optM$mPars) <- c(nmc * N, ninputpars)
  ## }

  if (is.null(data$nmdat$CMT)) data$nmdat$CMT <- 1 ## CHECKME
  if (any(is.na(data$nmdat$CMT))) {
    stop("'CMT' has NA(s)")
  }
  ## CHECKME
  form <- attr(model$saem_mod, "form")
  .nobs <- 0
  dat <- RxODE::etTrans(data$nmdat, attr(model$saem_mod, "rx"), TRUE, TRUE)
  .nobs <- attr(class(dat), ".RxODE.lst")$nobs
  ## if(length(dat) !=7) stop("SAEM doesn't support time varying covariates yet.");
  .rx <- attr(model$saem_mod, "rx")
  .pars <- .rx$params
  .pars <- setNames(rep(1.1, length(.pars)), .pars)
  opt$.rx <- .rx
  opt$.pars <- .pars
  ## opt$.dat <- dat;
  dat <- .as.data.frame(dat[, -6])
  names(dat) <- toupper(names(dat))
  dat$ID <- as.integer(dat$ID)

  evt <- dat
  evt$ID <- evt$ID - 1
  ## r
  evtM <- evt[rep(1:dim(evt)[1], nmc), ]
  evtM$ID <- cumsum(c(FALSE, diff(evtM$ID) != 0))

  i1 <- grep(1, diag(covstruct))
  i0 <- grep(0, diag(covstruct))
  nphi1 <- sum(diag(covstruct))
  nphi0 <- nphi - nphi1
  na <- length(mcmc$stepsize)
  nlambda1 <- sum(mcov[, i1])
  nlambda0 <- sum(mcov[, i0])
  nlambda <- nlambda1 + nlambda0
  nd1 <- nphi1 + nlambda1 + 1
  nd2 <- nphi1 + nlambda1 + nlambda0
  nb_param <- nd2 + 1
  Mcovariables <- cbind(rep(1, N), covariables)[, 1:nrow(mcov)]
  dim(Mcovariables) <- c(length(Mcovariables) / nrow(mcov), nrow(mcov)) # FIXME

  # get fixed ix
  fixed <- inits$theta.fix != ""
  wh <- fixed[, i1][mcov[, i1] == 1]
  len <- length(wh)
  fixed.i1 <- (1:len)[wh] - 1
  wh <- fixed[, i0][mcov[, i0] == 1]
  len <- length(wh)
  fixed.i0 <- (1:len)[wh] - 1

  jlog1 <- grep(T, model$log.eta)
  jcov <- grep(T, apply(mcov, 1, sum) > 0)
  covstruct1 <- covstruct[i1, i1]
  dim(covstruct1) <- c(nphi1, nphi1)
  ind_cov <- grep(1, mcov[mcov > 0])

  mcov1 <- matrix(mcov[, i1], ncol = length(i1))
  mcov0 <- matrix(mcov[, i0], nrow = nrow(mcov), ncol = length(i0))
  ind_cov1 <- grep(1, mcov1[mcov1 > 0]) - 1
  ind_cov0 <- grep(1, mcov0[mcov0 > 0]) - 1

  pc <- apply(mcov, 2, sum)
  ipc <- cumsum(c(0, pc[1:(nphi - 1)])) + 1
  ipcl1 <- ipc[jlog1]
  for (x in jlog1) inits$theta[1, x] <- log(inits$theta[1, x])

  idx <- as.vector(mcov1 > 0)
  COV1 <- Mcovariables[, row(mcov1)[idx]]
  dim(COV1) <- c(N, sum(idx))
  COV21 <- crossprod(COV1)

  x <- mcov1 * col(mcov1)
  x <- sapply(x[idx], function(x) {
    ret <- rep(0, nphi1)
    ret[x] <- 1
    ret
  })
  LCOV1 <- t(x)
  dim(LCOV1) <- c(nlambda1, nphi1)
  pc1 <- apply(LCOV1, 2, sum)

  x1 <- diag(sum(idx))
  diag(x1) <- inits$theta[, i1][idx]
  MCOV1 <- x1 %*% LCOV1
  jcov1 <- grep(1, LCOV1) - 1

  idx <- as.vector(mcov0 > 0)
  COV0 <- Mcovariables[, row(mcov0)[idx]]
  dim(COV0) <- c(N, sum(idx))
  COV20 <- crossprod(COV0)

  x <- mcov0 * col(mcov0)
  x <- sapply(x[idx], function(x) {
    ret <- rep(0, nphi0)
    ret[x] <- 1
    ret
  })
  LCOV0 <- t(x)
  dim(LCOV0) <- c(nlambda0, nphi0)

  x1 <- diag(sum(idx))
  diag(x1) <- inits$theta[, i0][idx]
  if (dim(x1)[1] > 0) {
    MCOV0 <- x1 %*% LCOV0
  } else {
    MCOV0 <- matrix(x1, nrow = 0, ncol = dim(LCOV0)[2])
  }
  jcov0 <- grep(1, LCOV0) - 1

  mprior_phi1 <- Mcovariables %*% inits$theta[, i1]
  mprior_phi0 <- Mcovariables %*% inits$theta[, i0]

  Gamma2_phi1 <- diag(nphi1)
  diag(Gamma2_phi1) <- inits$omega[i1]
  Gamma2_phi0 <- diag(nphi0)
  diag(Gamma2_phi0) <- inits$omega[i0]

  phiM <- matrix(0, N, nphi)
  phiM[, i1] <- mprior_phi1
  phiM[, i0] <- mprior_phi0
  phiM <- phiM[rep(1:N, nmc), , drop = FALSE]
  .tmp <- diag(sqrt(inits$omega))
  if (model$N.eta == 1) .tmp <- matrix(sqrt(inits$omega))
  if (normal == "vandercorput") {
    .mat2 <- matrix(RxODE::rxnormV(n=length(phiM)), dim(phiM))
  } else {
    .mat2 <- matrix(rnorm(phiM), dim(phiM))
  }
  phiM <- phiM + .mat2 %*% .tmp


  mc.idx <- rep(1:N, nmc)
  statphi <- sapply(1:nphi, function(x) tapply(phiM[, x], mc.idx, mean))
  statphi11 <- statphi[, i1]
  dim(statphi11) <- c(N, length(i1))
  statphi01 <- statphi[, i0]
  dim(statphi01) <- c(N, length(i0))
  statphi12 <- crossprod(phiM[, i1])
  statphi02 <- crossprod(phiM[, i0])

  # x = mcov *cumsum(mcov)
  # x1 = cbind(x[,i1], x[,i0])
  # indiphi = order(x1[x1>0])	#FINDME

  niter <- sum(mcmc$niter)
  niter_phi0 <- niter * .5
  nb_sa <- mcmc$niter[1] * .75
  nb_correl <- mcmc$niter[1] * .75
  va <- mcmc$stepsize
  vna <- mcmc$niter
  na <- length(va)
  pas <- 1 / (1:vna[1])^va[1]
  for (ia in 2:na) {
    end <- length(pas)
    k1 <- pas[end]^(-1 / va[ia])
    pas <- c(pas, 1 / ((k1 + 1):(k1 + vna[ia]))^va[ia])
  }
  pash <- c(rep(1, mcmc$burn.in), 1 / (1:niter))
  minv <- rep(1e-20, nphi)

  # preserve par order when printing iter history
  mcov[mcov == 1] <- 1:nlambda
  ilambda1 <- mcov[, i1]
  ilambda1 <- ilambda1[ilambda1 > 0] - 1
  ilambda0 <- mcov[, i0]
  ilambda0 <- ilambda0[ilambda0 > 0] - 1

  i1 <- i1 - 1
  i0 <- i0 - 1
  opt$distribution <- distribution
  opt$paramUpdate <- attr(model$saem_mod, "paramUpdate")
  optM$paramUpdate <- attr(model$saem_mod, "paramUpdate")
  opt$ODEopt <- ODEopt
  optM$ODEopt <- ODEopt
  cfg <- list(
    ODEopt=ODEopt,
    inits = inits.save,
    nu = mcmc$nu,
    niter = niter,
    nb_sa = nb_sa,
    nb_correl = nb_correl,
    niter_phi0 = niter_phi0,
    nmc = nmc,
    coef_phi0 = .9638, # FIXME
    rmcmc = .5,
    coef_sa = .95,
    pas = pas,
    pash = pash,
    minv = minv,

    N = N,
    ntotal = ntotal,
    y = y,
    yM = yM,
    phiM = phiM,
    evt = as.matrix(evt),
    evtM = as.matrix(evtM),
    mlen = mlen,
    indioM = indioM,

    pc1 = pc1,
    covstruct1 = covstruct1,
    Mcovariables = Mcovariables,

    i1 = i1,
    i0 = i0,
    nphi1 = nphi1,
    nphi0 = nphi0,
    nlambda1 = nlambda1,
    nlambda0 = nlambda0,
    COV0 = COV0,
    COV1 = COV1,
    COV20 = COV20,
    COV21 = COV21,
    LCOV0 = LCOV0,
    LCOV1 = LCOV1,
    MCOV0 = MCOV0,
    MCOV1 = MCOV1,
    Gamma2_phi0 = Gamma2_phi0,
    Gamma2_phi1 = Gamma2_phi1,
    mprior_phi0 = mprior_phi0,
    mprior_phi1 = mprior_phi1,
    jcov0 = jcov0,
    jcov1 = jcov1,
    ind_cov0 = ind_cov0,
    ind_cov1 = ind_cov1,
    statphi11 = statphi11,
    statphi01 = statphi01,
    statphi02 = statphi02,
    statphi12 = statphi12,
    res.mod = model$res.mod,
    ares = inits$ares,
    bres = inits$bres,
    cres = inits$cres,
    lres = inits$lres,
    opt = opt,
    optM = optM,
    print = mcmc$print,
    distribution = distribution,
    par.hist = matrix(0, sum(niter), nlambda1 + nlambda0 + nphi1 + 1 + (model$res.mod == 2) + 2 * (model$res.mod == 4)),
    seed = seed,
    fixed.i1 = fixed.i1,
    fixed.i0 = fixed.i0,
    ilambda1 = as.integer(ilambda1),
    ilambda0 = as.integer(ilambda0),
    nobs = .nobs
  )


  ## CHECKME
  s <- cfg$evt[cfg$evt[, "EVID"] == 0, "CMT"]
  cfg$opt$cmt_endpnt <- cfg$optM$cmt_endpnt <- sort(unique(s))
  cfg$nendpnt <- length(unique(s))
  if (model$nendpnt != cfg$nendpnt) {
    msg <- sprintf("mis-match in nbr endpoints in model & in data")
    stop(msg)
  }
  t <- unlist(split(1L:length(s), s))
  cfg$ysM <- rep(cfg$y[t], cfg$nmc)
  cfg$ix_sorting <- t - 1 # c-index for sorting by endpnt
  cfg$y_offset <- c(0, cumsum(table(s)))
  s <- cfg$evtM[cfg$evtM[, "EVID"] == 0, "CMT"]
  cfg$ix_endpnt <- as.integer(as.factor(s)) - 1 # to derive vecares & vecbres
  s <- cfg$evtM[cfg$evtM[, "EVID"] == 0, "ID"]
  t <- cumsum(c(0, table(s)))
  cfg$ix_idM <- cbind(t[-length(t)], t[-1] - 1) # c-index of obs records of each subject

  cfg$ares <- rep(10, cfg$nendpnt)
  cfg$bres <- rep(1, cfg$nendpnt)
  cfg$cres <- rep(1, cfg$nendpnt)
  cfg$lres <- rep(1, cfg$nendpnt)
  cfg$yj <- rep(2L, cfg$nendpnt)
  cfg$propT <- rep(0L, cfg$nendpnt)
  cfg$lambda <- rep(1.0, cfg$nendpnt)
  cfg$low <- rep(0.0, cfg$nendpnt)
  cfg$hi <- rep(1.0, cfg$nendpnt)
  cfg$ares[cfg$res.mod == 2] <- 0
  cfg$bres[cfg$res.mod == 1] <- 0
  nres <- (1:4)[(cfg$res.mod == 10L) * 3 + (cfg$res.mod %in% c(4L, 8L, 9L)) * 2 + (cfg$res.mod %in% c(3L, 5L, 6L, 7L)) + 1]
  cfg$res_offset <- cumsum(c(0, nres))
  cfg$par.hist <- matrix(0, cfg$niter, nlambda1 + nlambda0 + nphi1 + sum(nres))
  cfg$addProp <- c("combined1" = 1L, "combined2" = 2L)[match.arg(addProp)]

  cfg$DEBUG <- cfg$opt$DEBUG <- cfg$optM$DEBUG <- DEBUG
  cfg$phiMFile <- tempfile("phi-", RxODE::rxTempDir(), ".phi")
  cfg$tol <- tol
  cfg$itmax <- itmax
  cfg$type <- type
  cfg$lambdaRange <- lambdaRange
  cfg$powRange <- powRange
  cfg
}

#' Print an SAEM model fit summary
#'
#' Print an SAEM model fit summary
#'
#' @param object a saemFit object
#' @param ... others
#' @return a list
#' @export
summary.saemFit <- function(object, ...) {
  fit <- object ## Rcheck hack

  th <- fit$Plambda
  nth <- length(th)
  H <- solve(fit$Ha[1:nth, 1:nth])
  se <- sqrt(diag(H))

  m <- cbind(exp(th), th, se) # FIXME
  ## lhsVars = scan("LHS_VARS.txt", what="", quiet=T)
  ## if (length(lhsVars)==nth) dimnames(m)[[1]] = lhsVars
  dimnames(m)[[2]] <- c("th", "log(th)", "se(log_th)")
  cat("THETA:\n")
  print(m)
  cat("\nOMEGA:\n")
  print(fit$Gamma2_phi1)
  if (any(fit$sig2 == 0)) {
    cat("\nSIGMA:\n")
    print(max(fit$sig2^2))
  } else {
    cat("\nARES & BRES:\n")
    print(fit$sig2)
  }

  invisible(list(theta = th, se = se, H = H, omega = fit$Gamma2_phi1, eta = fit$mpost_phi))
}

#' Print an SAEM model fit summary
#'
#' Print an SAEM model fit summary
#'
#' @param x a saemFit object
#' @param ... others
#' @return a list
#' @export
print.saemFit <- function(x, ...) {
  fit <- x ## Rcheck hack

  th <- fit$Plambda
  nth <- length(th)
  H <- solve(fit$Ha[1:nth, 1:nth])
  se <- sqrt(diag(H))

  m <- cbind(exp(th), th, se) # FIXME
  ## lhsVars = scan("LHS_VARS.txt", what="", quiet=T)
  ## if (length(lhsVars)==nth) dimnames(m)[[1]] = lhsVars
  dimnames(m)[[2]] <- c("th", "log(th)", "se(log_th)")
  cat("THETA:\n")
  print(m)
  cat("\nOMEGA:\n")
  print(fit$Gamma2_phi1)
  if (any(fit$sig2 == 0)) {
    cat("\nSIGMA:\n")
    print(max(fit$sig2^2))
  } else {
    cat("\nARES & BRES:\n")
    print(fit$sig2)
  }

  invisible(list(theta = th, se = se, H = H, omega = fit$Gamma2_phi1, eta = fit$mpost_phi))
}


#' Fit an SAEM model
#'
#' Fit an SAEM model using either closed-form solutions or ODE-based model definitions
#'
#' @param model an RxODE model or lincmt()
#' @param data input data
#' @param inits initial values
#' @param PKpars PKpars function
#' @param pred  pred function
#' @param covars Covariates in data
#' @param mcmc a list of various mcmc options
#' @param ODEopt optional ODE solving options
#' @inheritParams configsaem
#' @param seed seed for random number generator
#' @details
#'    Fit a generalized nonlinear mixed-effect model using the Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Matthew Fidler & Wenping Wang
#' @export
saem.fit <- function(model, data, inits,
                     PKpars = NULL, pred = NULL,
                     covars = NULL,
                     mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                     ODEopt = list(atol = 1e-06, rtol = 1e-04, method = "lsoda", transitAbs = FALSE),
                     distribution = c("normal", "poisson", "binomial", "lnorm"),
                     seed = 99) {
  UseMethod("saem.fit")
}
##' @rdname saem.fit
##' @export
saem <- saem.fit

##' @rdname saem.fit
##' @export
saem.fit.nlmixr.ui.nlme <- function(model, data, inits,
                                    PKpars = NULL, pred = NULL,
                                    covars = NULL,
                                    mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                                    ODEopt = list(atol = 1e-06, rtol = 1e-04, method = "lsoda", transitAbs = FALSE),
                                    distribution = c("normal", "poisson", "binomial", "lnorm"),
                                    seed = 99) {
  call <- as.list(match.call(expand.dots = TRUE))[-1]
  names(call)[1] <- "object"
  call$est <- "saem"
  return(do.call(getFromNamespace("nlmixr", "nlmixr"), call, envir = parent.frame(1)))
}

##' @rdname saem.fit
##' @export
saem.fit.function <- saem.fit.nlmixr.ui.nlme

##' @rdname saem.fit
##' @export
saem.fit.nlmixrUI <- saem.fit.nlmixr.ui.nlme

##' @rdname saem.fit
##' @export
saem.fit.RxODE <- function(model, data, inits,
                           PKpars = NULL, pred = NULL,
                           covars = NULL,
                           mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                           ODEopt = list(atol = 1e-06, rtol = 1e-04, method = "lsoda", transitAbs = FALSE),
                           distribution = c("normal", "poisson", "binomial", "lnorm"),
                           seed = 99) {
  saem_fit <- gen_saem_user_fn(model, PKpars, pred)
  model <- list(saem_mod = saem_fit, covars = covars)
  cfg <- configsaem(model, data, inits, mcmc, ODEopt, distribution, seed)
  fit <- saem_fit(cfg)
  ## dyn.unload("saem_main.dll")
  fit
}

##' @rdname saem.fit
##' @export
saem.fit.default <- function(model, data, inits,
                             PKpars = NULL, pred = NULL,
                             covars = NULL,
                             mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                             ODEopt = list(atol = 1e-06, rtol = 1e-04, method = "lsoda", transitAbs = FALSE),
                             distribution = c("normal", "poisson", "binomial", "lnorm"),
                             seed = 99) {
  saem_fit <- gen_saem_user_fn(model)
  model <- list(saem_mod = saem_fit, covars = covars)
  cfg <- configsaem(model, data, inits, mcmc, ODEopt, distribution, seed)
  fit <- saem_fit(cfg)
  # dyn.unload("saem_main.dll")
  fit
}

##' @export
ranef.saemFit <- function(object, ...) {
  return(object$eta)
}

##' @export
fixef.saemFit <- function(object, ...) {
  return(object$Plambda)
}

focei.theta.saemFit <- function(object, uif, ...) {
  ## Get the thetas needed for FOCEi fit.
  this.env <- environment()
  if (class(uif) == "function") {
    uif <- nlmixr(uif)
  }
  n <- uif$focei.names
  thetas <- structure(rep(NA, length(n)), .Names = n)
  sf <- structure(as.vector(fixed.effects(object)), .Names = uif$saem.theta.name)
  for (n in names(sf)) {
    thetas[n] <- sf[n]
  }
  .predDf <- uif$predDf
  .ini <- .as.data.frame(uif$ini)
  .resMat <- object$resMat
  sapply(seq_along(.predDf$cond), function(i) {
    x <- paste(.predDf$cond[i])
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(sapply(.tmp$err, function(x) any(x == c("prop", "propT", "pow", "powT"))))
    if (length(.w) == 1) {
      thetas[paste(.tmp$name[.w])] <- .resMat[i, 2]
    }
    .w <- which(sapply(.tmp$err, function(x) any(x == c("pow2", "pow2T"))))
    if (length(.w) == 1) {
      thetas[paste(.tmp$name[.w])] <- .resMat[i, 3]
    }
    .w <- which(sapply(seq_along(.tmp$err), function(.x) {
      x <- .tmp$err[.x]
      if (any(x == c(
        "add", "norm", "dnorm", "lnorm", "dlnorm",
        "dlogn", "logn"))) {
        if (!is.na(.tmp$est[.x])) {
          return(TRUE)
        }
      }
      return(FALSE)
    }))
    if (length(.w) == 1) {
      thetas[paste(.tmp$name[.w])] <- .resMat[i, 1]
    }
    .w <- which(sapply(.tmp$err, function(x) {
      any(x == c("boxCox", "yeoJohnson"))
    }))
    if (length(.w) == 1) {
      thetas[paste(.tmp$name[.w])] <- .resMat[i, 4]
    }
    assign("thetas", thetas, this.env)
  })
  return(thetas)
}

focei.eta.saemFit <- function(object, uif, ...) {
  if (class(uif) == "function") {
    uif <- nlmixr(uif)
  }
  ## Reorder based on translation
  eta.trans <- uif$saem.eta.trans
  for (i in seq(1, max(eta.trans))) {
    while (!(any(i == eta.trans)) && max(eta.trans) > i) {
      eta.trans[eta.trans >= i] <- eta.trans[eta.trans >= i] - 1
    }
  }
  ## orig eta ->  new eta
  df <- .as.data.frame(uif$ini)
  eta <- df[!is.na(df$neta1), ]
  len <- length(eta$name)
  cur.lhs <- character()
  cur.rhs <- numeric()
  ome <- character()
  cur.ome <- object$Gamma2_phi1
  for (i in seq_along(eta$name)) {
    last.block <- FALSE
    if (i == len) {
      last.block <- TRUE
    } else if (eta$neta1[i + 1] == eta$neta2[i + 1]) {
      last.block <- TRUE
    }
    if (eta$neta1[i] == eta$neta2[i]) {
      cur.lhs <- c(cur.lhs, sprintf("ETA[%d]", eta$neta1[i]))
      cur.rhs <- c(cur.rhs, cur.ome[eta.trans[eta$neta1[i]], eta.trans[eta$neta2[i]]])
      if (last.block) {
        ome[length(ome) + 1] <- sprintf(
          "%s ~ %s", paste(cur.lhs, collapse = " + "),
          paste(deparse(cur.rhs), collapse = " ")
        )
        cur.lhs <- character()
        cur.rhs <- numeric()
      }
    } else {
      cur.rhs <- c(cur.rhs, cur.ome[eta.trans[eta$neta1[i]], eta.trans[eta$neta2[i]]])
    }
  }
  ome <- eval(parse(text = sprintf("list(%s)", paste(ome, collapse = ","))))
  return(ome)
}

as.focei.saemFit <- function(object, uif, pt = proc.time(), ..., data, calcResid = TRUE, obf = NULL,
                             nnodes.gq = 1, nsd.gq = 3, adjObf = TRUE,
                             calcCov = TRUE, covMethod = NULL,
                             calcCovTime = NULL, calcTables = TRUE) {
  .saemCfg <- attr(object, "saem.cfg")
  .saemTime <- proc.time() - pt
  if (class(uif) == "function") {
    uif <- nlmixr(uif)
  }
  .dist <- ""
  if (any(uif$saem.distribution == c("poisson", "binomial"))) {
    calcResid <- NA
    .dist <- uif$saem.distribution
  }
  uif.new <- uif
  fit <- object
  mat <- random.effects(fit)
  ## Reorder based on translation
  eta.trans <- uif$saem.eta.trans
  for (i in seq(1, max(eta.trans))) {
    while (!(any(i == eta.trans)) && max(eta.trans) > i) {
      eta.trans[eta.trans >= i] <- eta.trans[eta.trans >= i] - 1
    }
  }
  mat2 <- mat[, eta.trans, drop = FALSE]
  th <- focei.theta(fit, uif)
  for (n in names(th)) {
    uif.new$est[uif.new$name == n] <- th[n]
  }
  ome <- focei.eta(fit, uif)
  init <- list(
    THTA = as.vector(th),
    OMGA = ome
  )
  saem.time <- proc.time() - pt
  if (missing(data)) {
    stop("Requires Data...")
  } else {
    dat <- data
  }
  .tn <- uif$saem.theta.name
  .ini <- .as.data.frame(uif$ini)
  .ini <- .ini[uif$ini$name %in% .tn, ]
  if (any(.ini$fix)) {
    .fixed <- paste(.ini$name[.ini$fix])
    .tn <- .tn[!(.tn %in% .fixed)]
  }
  .calcCov <- calcCov
  if (uif$env$covMethod == "") {
    .cov <- NULL
    .addCov <- FALSE
  } else if (inherits(calcCov, "matrix")) {
    .cov <- calcCov
    .addCov <- TRUE
  } else {
    .nth <- length(.tn)

    .ini <- .as.data.frame(uif$ini)
    .ini <- .ini[is.na(.ini$err), ]
    .ini <- .ini[!is.na(.ini$ntheta), ]
    .ini <- .ini[!.ini$fix, ]
    .ini <- paste(.ini$name)
    .calcCovTime <- proc.time()
    if (calcCov) {
      .covm <- object$Ha[1:.nth, 1:.nth]
      .covm <- try(calc.COV(object))
      .doIt <- !inherits(.covm, "try-error")
      if (.doIt && dim(.covm)[1] != .nth) .doIt <- FALSE
      if (.doIt) {
        .tmp <- try(chol(.covm), silent = TRUE)
        .addCov <- TRUE
        .sqrtm <- FALSE
        if (inherits(.tmp, "try-error")) {
          .tmp <- .covm
          .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
          if (inherits(.tmp, "try-error")) {
            .calcCov <- FALSE
            .covm <- object$Ha[1:.nth, 1:.nth]
            .tmp <- try(chol(.covm), silent = TRUE)
            .addCov <- TRUE
            .sqrtm <- FALSE
            if (inherits(.tmp, "try-error")) {
              .tmp <- object$Ha[1:.nth, 1:.nth]
              .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
              if (inherits(.tmp, "try-error")) {
                .addCov <- FALSE
              } else {
                .sqrtm <- TRUE
              }
            } else {
              .tmp <- object$Ha[1:.nth, 1:.nth]
            }
          } else {
            .sqrtm <- TRUE
          }
        } else {
          .tmp <- .covm
        }
      } else {
        .tmp <- object$Ha[1:.nth, 1:.nth]
        .tmp <- try(chol(.covm), silent = TRUE)
        .calcCov <- FALSE
        .addCov <- TRUE
        .sqrtm <- FALSE
        if (inherits(.tmp, "try-error")) {
          .tmp <- object$Ha[1:.nth, 1:.nth]
          .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
          if (inherits(.tmp, "try-error")) {
            .addCov <- FALSE
          } else {
            .sqrtm <- TRUE
          }
        } else {
          .tmp <- object$Ha[1:.nth, 1:.nth]
          .calcCov <- FALSE
        }
      }
    } else {
      .tmp <- try(chol(.covm), silent = TRUE)
      .addCov <- TRUE
      .sqrtm <- FALSE
      if (inherits(.tmp, "try-error")) {
        .tmp <- object$Ha[1:.nth, 1:.nth]
        .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
        if (inherits(.tmp, "try-error")) {
          .addCov <- FALSE
        } else {
          .sqrtm <- TRUE
        }
      } else {
        .tmp <- object$Ha[1:.nth, 1:.nth]
        .calcCov <- FALSE
      }
    }
    if (.addCov) {
      if (!.calcCov) {
        .cov <- RxODE::rxInv(.tmp)
      } else {
        .cov <- .tmp
      }
      attr(.cov, "dimnames") <- list(.tn, .tn)
      .cov <- .cov[.ini, .ini, drop = FALSE]
    }
    .calcCovTime <- proc.time() - .calcCovTime
    .calcCovTime <- .calcCovTime["elapsed"]
  }
  .ini <- .as.data.frame(uif$ini)
  .ini <- .ini[!is.na(.ini$ntheta), ]
  .skipCov <- !is.na(.ini$err)
  .fixed <- uif$focei.fixed
  .skipCov <- .skipCov | .fixed[seq_along(.skipCov)]
  .covMethod <- uif$env$covMethod
  if (!any(.covMethod == c("", "r", "s", "r,s"))) {
    .covMethod <- ""
  }
  if (is.na(calcResid)) .covMethod <- ""
  .allThetaNames <- c(uif$saem.theta.name, uif$saem.omega.name, uif$saem.res.name)
  .m <- object$par_hist
  if (ncol(.m) > length(.allThetaNames)) {
    .m <- .m[, seq_along(.allThetaNames)]
  }
  if (.dist == "binomial") {
    .dist <- "bernoulli"
  }
  dimnames(.m) <- list(NULL, .allThetaNames)
  .fixedNames <- paste(uif$ini$name[which(uif$ini$fix)])
  .rn <- ""
  .likTime <- 0
  if (is.na(obf)) {
    .saemObf <- NA
  } else if (is.null(obf)) {
    .likTime <- proc.time()
    .saemObf <- calc.2LL(object, nnodes.gq = nnodes.gq, nsd.gq = nsd.gq)
    .likTime <- proc.time() - .likTime
    .likTime <- .likTime["elapsed"]
    if (nnodes.gq == 1) {
      .rn <- paste0("laplace", nsd.gq)
    } else {
      .rn <- paste0("gauss", nnodes.gq, "_", nsd.gq)
    }
  } else if (is(obf, "logical")) {
    if (is.na(obf)) {
      .saemObf <- NA
    } else if (obf) {
      .likTime <- proc.time()
      .saemObf <- calc.2LL(object, nnodes.gq = nnodes.gq, nsd.gq = nsd.gq)
      .likTime <- proc.time() - .likTime
      .likTime <- .likTime["elapsed"]
      if (nnodes.gq == 1) {
        .rn <- paste0("laplace", nsd.gq)
      } else {
        .rn <- paste0("gauss", nnodes.gq, "_", nsd.gq)
      }
    } else {
      .saemObf <- NA
    }
  } else if (is(object, "numeric")) {
    .saemObf <- obf
  }
  .notCalced <- TRUE
  .cwresTime <- proc.time()
  while (.notCalced) {
    .env <- new.env(parent = emptyenv())
    .env$nobs2 <- .saemCfg$nobs
    .env$nnodes.gq <- nnodes.gq
    .env$nsd.gq <- nsd.gq
    .env$adjObf <- adjObf
    .env$method <- "SAEM"
    .env$uif <- uif
    .env$saem <- object
    if (.addCov) {
      .env$cov <- .cov
    }
    .env$parHistStacked <- .data.frame(
      val = as.vector(.m),
      par = rep(.allThetaNames, each = nrow(.m)),
      iter = rep(1:nrow(.m), ncol(.m))
    )
    .env$parHist <- .data.frame(iter = rep(1:nrow(.m)), .as.data.frame(.m))
    if (length(.fixedNames) > 0) {
      .env$parHistStacked <- .env$parHistStacked[!(.env$parHistStacked$par %in% .fixedNames), , drop = FALSE]
      .env$parHist <- .env$parHist[, !(names(.env$parHist) %in% .fixedNames), drop = FALSE]
    }
    if (is.na(calcResid)) {
      .setSaemExtra(.env, .rn)
      .env$theta <- .data.frame(
        lower = -Inf, theta = init$THTA, upper = Inf, fixed = .fixed[seq_along(init$THTA)],
        row.names = uif$focei.names
      )
      .env$fullTheta <- setNames(init$THTA, uif$focei.names)
      .om0 <- .genOM(.parseOM(init$OMGA))
      attr(.om0, "dimnames") <- list(uif$eta.names, uif$eta.names)
      .env$omega <- .om0
      .env$etaObf <- .data.frame(ID = seq_along(mat2[, 1]), setNames(.as.data.frame(mat2), uif$eta.names), OBJI = NA)
      .env$noLik <- FALSE
      .env$objective <- .saemObf
    } else if (calcResid) {
      .setSaemExtra(.env, .rn)
    } else {
      .setSaemExtra(.env, .rn)
      .env$theta <- .data.frame(
        lower = -Inf, theta = init$THTA, upper = Inf, fixed = .fixed[seq_along(init$THTA)],
        row.names = uif$focei.names
      )
      .env$fullTheta <- setNames(init$THTA, uif$focei.names)
      .om0 <- .genOM(.parseOM(init$OMGA))
      attr(.om0, "dimnames") <- list(uif$eta.names, uif$eta.names)
      .env$omega <- .om0
      .env$etaObf <- .data.frame(ID = seq_along(mat2[, 1]), setNames(.as.data.frame(mat2), uif$eta.names), OBJI = NA)
      .env$noLik <- TRUE
      .env$objective <- .saemObf
    }
    .ctl <- uif$env$ODEopt
    names(.ctl) <- sub("maxsteps", "maxstepsOde", names(.ctl))
    .ctl <- .ctl[names(.ctl) != "scale"]
    .ctl$maxOuterIterations <- 0
    .ctl$maxInnerIterations <- 0
    .ctl$covMethod <- .covMethod
    .ctl$sumProd <- uif$env$sum.prod
    .ctl$optExpression <- uif$env$optExpression
    .ctl$scaleTo <- 0
    .ctl$calcTables <- calcTables
    if (.saemCfg$addProp == 1L) {
      .ctl$addProp <- "combined1"
    } else {
      .ctl$addProp <- "combined2"
    }
    .ctl <- do.call(foceiControl, .ctl)
    if (.ctl$singleOde) {
      .pars <- NULL
      .mod <- uif$focei.rx1
    } else {
      .pars <- uif$theta.pars
      .mod <- uif$rxode.pred
    }
    fit.f <- try(foceiFit.data.frame(
      data = dat, inits = init, PKpars = .pars,
      model = .mod, pred = function() {
        return(nlmixr_pred)
      }, err = uif$error,
      lower = uif$focei.lower,
      upper = uif$focei.upper,
      thetaNames = uif$focei.names,
      etaNames = uif$eta.names,
      etaMat = mat2,
      env = .env,
      fixed = .fixed,
      skipCov = .skipCov,
      control = .ctl
    ), silent = FALSE)
    if (inherits(fit.f, "try-error")) {
      if (is.na(calcResid)) {
        warning("Error calculating nlmixr object, return classic object")
        .notCalced <- FALSE
        return(object)
      } else if (calcResid) {
        calcResid <- FALSE
      } else {
        calcResid <- NA
      }
    } else {
      .notCalced <- FALSE
    }
  }
  .cwresTime <- proc.time() - .cwresTime
  if (is.na(calcResid)) {
    .cwresTime <- 0
  } else if (!calcResid) {
    .cwresTime <- 0
  }
  .env <- fit.f$env
  if (uif$env$covMethod == "") {
  } else if (inherits(.calcCov, "matrix")) {
    if (!is.null(covMethod)) {
      .env$covMethod <- covMethod
    }
    .calcCovTime <- calcCovTime
  } else if (.calcCov) {
    .env$covMethod <- "linFim"
    if (.addCov & .sqrtm) {
      .env$covMethod <- "|linFim|"
      warning("Covariance matrix non-positive definite, corrected by sqrtm(linFim %*% linFim)")
    }
  } else {
    if (calcCov) {
      warning("Linearization of FIM could not be used to calculate covariance.")
    }
    if (.addCov & .sqrtm) {
      .env$covMethod <- "|fim|"
      warning("Covariance matrix non-positive definite, corrected by sqrtm(fim %*% fim)")
    } else if (!.addCov) {
      warning("FIM non-positive definite and cannot be used to calculate the covariance")
    }
  }

  if (is.null(.env$time)) {
    .env$time <- .data.frame(saem = .saemTime["elapsed"], check.names = FALSE, row.names = c(""))
  } else {
    .time <- .env$time
    .time <- .time[, !(names(.time) %in% c("optimize", "covariance"))]
    .saemTime <- .saemTime["elapsed"]
    if (calcResid) {
      .saemTime <- .saemTime - .cwresTime["elapsed"]
      .time <- .data.frame(.time, cwres = .cwresTime["elapsed"], check.names = FALSE)
    }
    if (.likTime > 0) {
      .time <- .data.frame(.time, logLik = .likTime, check.names = FALSE)
      .saemTime <- .saemTime - .likTime
    }
    if (uif$env$covMethod != "") {
      .saemTime <- .saemTime - .calcCovTime
      .time <- .data.frame(.time, covariance = .calcCovTime, check.names = FALSE)
    }
    .env$time <- .data.frame(saem = .saemTime, .time, check.names = FALSE, row.names = c(""))
  }
  .env$message <- ""
  if (is.na(calcResid)) {
    row.names(.env$objDf) <- .rn
  } else if (calcResid) {
    if (!is.na(.saemObf)) {
      .llik <- -.saemObf / 2
      .nobs <- .env$nobs
      attr(.llik, "df") <- attr(get("logLik", .env), "df")
      .objf <- ifelse(.env$adjObf, .saemObf - .nobs * log(2 * pi), .saemObf)
      ## for (.t in c("OBJF","objective", "objf")){
      ##   assign(.t,.objf,.env);
      ## }
      .tmp <- .data.frame(
        OBJF = .objf,
        AIC = .saemObf + 2 * attr(get("logLik", .env), "df"),
        BIC = .saemObf + log(.env$nobs) * attr(get("logLik", .env), "df"),
        "Log-likelihood" = as.numeric(.llik), check.names = FALSE
      )
      if (any(names(.env$objDf) == "Condition Number")) {
        .tmp <- .data.frame(.tmp,
          "Condition Number" = .env$objDf[, "Condition Number"],
          check.names = FALSE
        )
      }
      .env$objDf <- rbind(.env$objDf, .tmp)
      row.names(.env$objDf) <- c("FOCEi", .rn)
    }
    .setSaemExtra(.env, "FOCEi")
  } else {
    row.names(.env$objDf) <- .rn
  }
  if (inherits(fit.f, "nlmixrFitData")) {
    .cls <- class(fit.f)
    .env <- attr(.cls, ".foceiEnv")
    .cls[1] <- "nlmixrSaem"
    class(fit.f) <- .cls
  }
  assign("uif", .syncUif(fit.f$uif, fit.f$popDf, fit.f$omega), fit.f$env)
  return(fit.f)
}

.saemResidF <- function(x) {
  .Call(`_saemResidF`, x)
}

.newuoa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- control
  if (is.null(.ctl$npt)) .ctl$npt <- length(par) * 2 + 1
  if (is.null(.ctl$rhobeg)) .ctl$rhobeg <- 0.2
  if (is.null(.ctl$rhoend)) .ctl$rhoend <- 1e-4
  .ctl$iprint <- 0L
  .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")]
  .ret <- minqa::newuoa(par, fn,
    control = .ctl
  )
  .ret$x <- .ret$par
  .ret$message <- .ret$msg
  .ret$convergence <- .ret$ierr
  .ret$value <- .ret$fval
  return(.ret)
}


## FIXME: coef_phi0, rmcmc, coef_sa
## FIXME: Klog, rho, sa, nmc
## FIXME: N.design
## FIXME: g = gc = 1
## FIXME: ODE inits
## FIXME: Tinf for ODE
## FIXME: chk infusion poor fit

## Local Variables:
## ess-indent-level: 2
## indent-tabs-mode: nil
## End:
