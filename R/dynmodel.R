thresh <- function(x, cut = .Machine$double.xmin) {
  x <- abs(x)
  ifelse(x > cut, x, cut)
}

err.msg <- function(x, pre = "", post = "") {
  msg <- paste0(x, collapse = ", ")
  paste0(pre, msg, post)
}
##' Convert fit to classic dynmodel object
##'
##' @param x nlmixr object to convert to dynmodel object
##' @return dynmodel
##' @author Matthew Fidler
##' @export
as.dynmodel <- function(x) {
  .ret <- try(x$dynmodelObject, silent = TRUE)
  if (inherits(.ret, "try-error") || is.null(.ret)) {
    stop("Cannot convert to dynmodel object")
  }
  return(.ret)
}

# plot.dyn.ID() -----------------------------------------------------------
#' Plot of a non-population dynamic model fit
#'
#' @param x a dynamodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
gof <- function(x, ...) {
  gof_str <- "
  {
  theta = th[1:npar]
  names(theta) = names(inits)[1:npar]
  theta = c(theta, fixPars)
  system$solve(theta, ev, atol = 1e-08, rtol = 1e-08)
  }
  "
  f <- x$obj
  dat <- x$data
  body(f) <- parse(text = gof_str)
  x <- f(x$par)

  rows <- environment(f)$rows
  m <- environment(f)$model
  for (s in m) {
    time <- x[, "time"]
    yo <- dat[, s["dv"]] # FIXME
    yp <- x[, s["pred"]]

    # dv vs pred
    plot(time, yp, type = "n", xlab = "time", ylab = s["dv"])
    points(dat[, "time"], yo, ...)
    lines(time, yp)

    # pred vs res
    yp <- yp[rows]
    res <- yo - yp
    plot(yp, yp, xlab = "pred", ylab = s["dv"], ...)
    abline(0, 1, col = "red", lty = 2)

    # pred vs res
    plot(yp, res, xlab = "pred", ylab = "res", ...)
    abline(h = 0, col = "red", lty = 2)
  }
}

#' @export
#' @rdname gof
plot.dyn.ID <- gof
# #########################################################################

# print.dyn.ID() ----------------------------------------------------------
#' Print a non-population dynamic model fit object
#'
#' @param x a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
print.dyn.ID <- function(x, ...) {
  print(x$res[, c(1, 3)])
  cat("\n")
  aic <- 2 * x$value + 2 * x$npar
  bic <- 2 * x$value + log(x$nobs) * x$npar
  print(c("-loglik" = x$value, "AIC" = aic, "BIC" = bic))
  cat("\n")
}

# #########################################################################

# summary.dyn.ID() --------------------------------------------------------
#' Summary of a non-population dynamic model fit
#'
#' @param object a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
summary.dyn.ID <- function(object, ...) {
  print(object$res)
  cat("\n")
  aic <- 2 * object$value + 2 * object$npar
  bic <- 2 * object$value + log(object$nobs) * object$npar
  print(c("-loglik" = object$value, "AIC" = aic, "BIC" = bic))
  cat("\niter:", object$iter, "\n")
  cat(object$message, "\n")
}

# nmsimplex() and mymin() -------------------------------------------------
#' Nelder-Mead simplex search
#'
#' @param start initials
#' @param fr objective function
#' @param rho evaluation environment
#' @param control additional optimization options
#' @return a list of ...
#' @export
nmsimplex <- function(start, fr, rho = NULL, control = list()) {
  if (is.null(rho)) rho <- environment(fr)
  step <- -.2 * start

  con <- list(maxeval = 999, reltol = 1e-6, rcoeff = 1., ecoeff = 2., ccoeff = .5, trace = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) {
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }

  .Call(neldermead_wrap, fr, rho, length(start), start, step,
    as.integer(con$maxeval), con$reltol, con$rcoeff, con$ecoeff, con$ccoeff,
    as.integer(con$trace),
    PACKAGE = "nlmixr"
  )
}

mymin <- function(start, fr, rho = NULL, control = list()) {
  if (is.null(rho)) rho <- environment(fr)
  step <- -.2 * start

  con <- list(maxeval = 999, ftol_rel = 1e-6, rcoeff = 1., ecoeff = 2., ccoeff = .5, trace = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  # if (length(noNms <- namc[!namc %in% nmsC]))
  #    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  .Call(neldermead_wrap, fr, rho, length(start), start, step,
    as.integer(con$maxeval), con$reltol, con$rcoeff, con$ecoeff, con$ccoeff,
    as.integer(con$trace),
    PACKAGE = "nlmixr"
  )
}


# #########################################################################

# as.focei.dynmodel() -----------------------------------------------------------
#' Output nlmixr format for dynmodel
#'
#' Convert dynmodel output to nlmixr focei style output
#'
#' @param .dynmodelObject dynmodel object
#' @param .nlmixrObject nlmixr object
#' @param .data RxODE data set
#' @param .time proc.time object used in dynmodel timing
#' @param .theta estimated terms excluding error terms
#' @param .fit optimized parameters
#' @param .message optimization messages
#' @param .inits.err initial estimates for error terms
#' @param .cov covariance matrix from solved Hessian
#' @param .sgy variance of Y given the model
#' @param .dynmodelControl control options for dynmodel
#' @param .nobs2 -
#' @param .pt proc.time object used in dynmodel timing
#' @param .rxControl control options for RxODE
#' @author Mason McComb and Matt Fidler
#' @export
as.focei.dynmodel <- function(.dynmodelObject, .nlmixrObject, .data,
                              .time, .theta, .fit, .message, .inits.err,
                              .cov, .sgy, .dynmodelControl, .nobs2 = 0,
                              .pt = proc.time(),
                              .rxControl = RxODE::rxControl()) {
  .Estimates <- .dynmodelObject$res[, 1]
  .model <- RxODE::rxSymPySetupPred(.nlmixrObject$rxode.pred,
    function() {
      return(nlmixr_pred)
    },
    .nlmixrObject$theta.pars,
    .nlmixrObject$error,
    grad = FALSE,
    pred.minus.dv = TRUE, sum.prod = FALSE, # control$sumProd,
    theta.derivs = FALSE, optExpression = TRUE, # control$optExpression,
    run.internal = TRUE, only.numeric = TRUE
  )
  .temp.model <- .model
  .nlmixrObject.df <- as.data.frame(.nlmixrObject$ini)
  # assign fixed terms
  .fix.index <-
    if (length(which(.nlmixrObject.df$fix == TRUE)) == 0) {
      NULL
    } else {
      which(.nlmixrObject.df$fix == TRUE)
    } # obtain row location for fixed terms
  .ref.fix <- substring(.nlmixrObject.df$name[.fix.index], 2)
  .temp.log.fixPars.index <- intersect(.fix.index, .temp.model$log.etas)
  .temp.log.fixPars <- .nlmixrObject.df$est[.temp.log.fixPars.index]
  names(.temp.log.fixPars) <- substring(.nlmixrObject.df$name[.temp.log.fixPars.index], 2)

  .temp.nonlog.fixPars.index <- setdiff(.fix.index, .temp.model$log.etas)
  .temp.nonlog.fixPars <- .nlmixrObject.df$est[.temp.nonlog.fixPars.index]
  names(.temp.nonlog.fixPars) <- substring(.nlmixrObject.df$name[.temp.nonlog.fixPars.index], 2)

  .fixPars <- c(.temp.log.fixPars, .temp.nonlog.fixPars)
  .fixPars <- .fixPars[order(factor(names(.fixPars), levels = .ref.fix))]

  # assign theta terms(estimated terms excluding error terms)
  .theta.index <-
    if (is.null(.fix.index)) {
      which(!is.na(.nlmixrObject.df$ntheta) & is.na(.nlmixrObject.df$err), TRUE)
    } else {
      which(!is.na(.nlmixrObject.df$ntheta) & is.na(.nlmixrObject.df$err), TRUE)[-which(.nlmixrObject.df$fix == TRUE)]
    }

  ## $nlmixrObject ----
  .temp.theta.index <-
    c(1:sum(
      as.data.frame(
        .nlmixrObject$ini
      )$fix == FALSE &
        is.na(as.data.frame(.nlmixrObject$ini)$err)
    ))
  .temp.replacements <- .dynmodelObject$res[.theta.index, 1]
  ini <- as.data.frame(.nlmixrObject$ini)
  ini$est[.temp.theta.index] <- c(.temp.replacements)
  class(ini) <- c("nlmixrBounds", "data.frame")
  .nlmixrObject$ini <- ini
  uif <- .nlmixrObject

  ## now use focei to get the same estimates
  .env <- new.env(parent = emptyenv())
  ## put in covariance information
  labels <- names(.Estimates)
  dimnames(.cov) <- list(labels, labels)
  .env$cov <- .cov
  .env$nobs <- .dynmodelObject$nobs
  .env$covMethod <- .dynmodelControl$covMethod
  .ctl <- .dynmodelControl
  class(.ctl) <- NULL
  .ctl$maxOuterIterations <- 0
  .ctl$maxInnerIterations <- 0
  .ctl$covMethod <- ""
  .ctl$sumProd <- uif$env$sum.prod
  .ctl$optExpression <- uif$env$optExpression
  .ctl$scaleTo <- 0
  .ctl$calcTables
  .ctl$method <-  "liblsoda"
  .ctl <- do.call(foceiControl, .ctl)
  if (.ctl$singleOde) {
    .mod <- uif$focei.rx1
    .pars <- NULL
  } else {
    .mod <- uif$rxode.pred
    .pars <- uif$theta.pars
  }
  fit.f <- try(foceiFit.data.frame(
      data = .data, inits = uif$focei.inits, PKpars = .pars,
      model = .mod, pred = function() {
        return(nlmixr_pred)
      }, err = uif$error,
      lower = uif$focei.lower,
      upper = uif$focei.upper,
      thetaNames = uif$focei.names,
      etaNames = uif$eta.names,
      fixed = uif$focei.fixed,
      env = .env,
      control = .ctl
  ), silent = FALSE)
  assign("method", "dynmodel", .env)
  assign("extra", sprintf("(method: %s)",  .dynmodelControl$method), .env)
  assign("dynmodelObject",  .dynmodelObject, .env)
  .cls <- c("nlmixrDynmodel", class(fit.f))
  attr(.cls, ".foceiEnv") <- .env
  class(fit.f) <- .cls
  return(fit.f)
}

# nlmixrDynmodelConvert() #################################################
#' Converting nlmixr objects to dynmodel objects
#'
#' Convert nlmixr Objects to dynmodel objects for use in fitting non-population dynamic models
#'
#' @param .nmf nlmixr object
#' @return list containing inputs for the dynmodel()
#' \itemize{
#' \item \code{$fixPars} - fixed parameters defined as \code{fixed()} in the nlmixr object
#' \item \code{$sigma} - error model parameters
#' \item \code{$inits} - initial estimates for parameters in the model
#' \item \code{$lower} - lower boundaries for estimated parameters
#' \item \code{$upper} - upper boundaries for estimated parameters
#' \item \code{$system} - RxODE object that defines the structural model
#' \item \code{$model} - error model
#' }
#'
#' @author Mason McComb and Matt Fidler
#' @export
nlmixrDynmodelConvert <- function(.nmf) {
  # initialize list for output
  .return <- list()

  # convert nlmixr model to data.frame
  .nmf.original <- .nmf
  .nmf <- as.data.frame(.nmf$ini)
  .temp.model <-
    RxODE::rxSymPySetupPred(
      .nmf.original$rxode.pred,
      function() {
        return(nlmixr_pred)
      },
      .nmf.original$theta.pars,
      .nmf.original$error,
      grad = FALSE,
      pred.minus.dv = TRUE,
      sum.prod = FALSE, # control$sumProd,
      theta.derivs = FALSE,
      optExpression = TRUE, # control$optExpression,
      run.internal = TRUE,
      only.numeric = TRUE
    )

  # assign fixed terms
  .fix.index <- # obtain row location for fixed terms
    if (length(which(.nmf$fix == TRUE)) == 0) {
      NULL
    } else {
      which(.nmf$fix == TRUE)
    }
  .ref.fix <- substring(.nmf$name[.fix.index], 2)

  .temp.log.fixPars.index <- intersect(.fix.index, .temp.model$log.etas)
  .temp.log.fixPars <- .nmf$est[.temp.log.fixPars.index]
  names(.temp.log.fixPars) <- substring(.nmf$name[.temp.log.fixPars.index], 2)

  .temp.nonlog.fixPars.index <- setdiff(.fix.index, .temp.model$log.etas)
  .temp.nonlog.fixPars <- .nmf$est[.temp.nonlog.fixPars.index]
  names(.temp.nonlog.fixPars) <- substring(.nmf$name[.temp.nonlog.fixPars.index], 2)

  .fixPars <- c(.temp.log.fixPars, .temp.nonlog.fixPars)
  .fixPars <- .fixPars[order(factor(names(.fixPars), levels = .ref.fix))]

  .return <- c(.return, fixPars = list(.fixPars))

  # assign theta terms(estimated terms excluding error terms)
  .theta.index <-
    if (is.null(.fix.index)) {
      which(!is.na(.nmf$ntheta) & is.na(.nmf$err), TRUE)
    } else {
      which(!is.na(.nmf$ntheta) & is.na(.nmf$err), TRUE)[-which(.nmf$fix == TRUE)]
    } # row location for theta values

  .ref.theta <- .nmf$name[.theta.index]

  .temp.log.theta.index <- intersect(.theta.index, .temp.model$log.etas)
  .temp.log.theta <- .nmf$est[.temp.log.theta.index]
  names(.temp.log.theta) <- .nmf$name[.temp.log.theta.index]

  .temp.nonlog.theta.index <- setdiff(.theta.index, .temp.model$log.etas)
  .temp.nonlog.theta <- .nmf$est[.temp.nonlog.theta.index]
  names(.temp.nonlog.theta) <- .nmf$name[.temp.nonlog.theta.index]

  .theta <- c(.temp.nonlog.theta, .temp.log.theta)
  .theta <- .theta[order(factor(names(.theta), levels = .ref.theta))]

  # .theta.name.index <-  .theta.index + 2 + length(.theta.index) + length(.fix.index)

  # names(.theta) <- .temp.model$pred.only$lhs[.theta.name.index]

  # assign sigma terms (estimated)
  .sigma.index <-
    if (is.null(.fix.index)) {
      which(!is.na(.nmf$ntheta) & !is.na(.nmf$err), TRUE)
    } else {
      which(!is.na(.nmf$ntheta) & !is.na(.nmf$err), TRUE)[-which(.nmf$fix == TRUE)]
    } # row location for theta values

  .sigma <- .nmf$est[.sigma.index] # initial estimate for theta values, back-transformed

  names(.sigma) <- .nmf$err[.sigma.index]
  .return <- c(.return, sigma = list(.sigma))

  # assign "inits", vector of theta and sigma terms # (will be used when the likelihood function is changed)
  .inits <- .theta
  .return$inits <- .inits

  # assign boundaries
  mask_theta_term <- !is.na(.nmf$ntheta) & !.nmf$fix & is.na(.nmf$err)
  mask_error_term <- !is.na(.nmf$ntheta) & !.nmf$fix & !is.na(.nmf$err)
  .lower <-
    c(
      .nmf$lower[mask_theta_term], # theta terms
      .nmf$lower[mask_error_term] # error terms
    )
  .upper <-
    c(
      .nmf$upper[mask_theta_term], # theta terms
      .nmf$upper[mask_error_term] # error terms
    )

  .return <- c(.return, lower = list(.lower), upper = list(.upper))

  # obtain system
  .system <- RxODE(.nmf.original$rxode.pred) # (use nlmixr_pred)
  .system$stateExtra <- NULL # remove extraState, the error model term should not be inclcuded
  .system$lhs <- .system$lhs[-length(.system$lhs)] # remove the error model term
  .return <- c(.return, system = .system)

  # create error model
  .DV <- .nmf$condition[!is.na(.nmf$condition) & .nmf$condition != "ID"]
  .PRED <- "nlmixr_pred" # need to obtain from data? id dont know

  .formula <- list()
  for (i in 1:length(.sigma.index)) {
    if (i == 1) {
      .temp <- paste(.nmf$err[.sigma.index[i]], "(", .nmf$est[.sigma.index[i]], ")", sep = "")
    } else {
      .temp <- paste("+ ", .nmf$err[.sigma.index[i]], "(", .nmf$est[.sigma.index[i]], ")", sep = "")
    }
    .formula <- paste(.formula, .temp)
  }

  .model <- as.formula(paste(.DV, "~", .PRED, "+", .formula))
  .return <- c(.return, model = .model)

  # Output
  return(.return)
}

# dynmodelControl() #######################################################
#' Control Options for dynmodel
#'
#' @inheritParams RxODE::rxSolve
#' @inheritParams foceiControl
#'
#' @param nlmixrOutput Option to change output style to nlmixr output. By default
#'   this is FALSE.
#' @param digs Option for the number of significant digits of the output. By
#'   default this is 3.
#' @param lower Lower bounds on the parameters used in optimization. By default
#'   this is -Inf.
#' @param upper Upper bounds on the parameters used in optimization. By default
#'   this is Inf.
#' @param maxeval Maximum number of iterations for Nelder-Mead of simplex
#'   search. By default this is 999.
#' @param maxit Maximum number of iterations for lbfgsb3c. See
#'   \code{\link[lbfgsb3c]{lbfgsb3c}} for more details. By default this is
#'   100000L.
#' @param maxfun The maximum allowed number of function evaluations. If this is
#'   exceeded, the method will terminate. See \code{\link[minqa]{bobyqa}} for
#'   more details. By default this value is NULL.
#' @param iprint Print option for optimization. See \code{\link[minqa]{bobyqa}},
#'   \code{\link[lbfgsb3c]{lbfgsb3c}}, and \code{\link[lbfgs]{lbfgs}} for more
#'   details. By default this is 0.
#' @param trace Tracing information on the progress of the optimization is
#'   produced. See \code{\link[minqa]{bobyqa}},
#'   \code{\link[lbfgsb3c]{lbfgsb3c}}, and \code{\link[lbfgs]{lbfgs}} for more
#'   details. By default this is 0.
#' @param factr Controls the convergence of the "L-BFGS-B" method.  Convergence
#'   occurs when the reduction in the objective is within this factor of the
#'   machine tolerance. Default is 1e10, which gives a tolerance of about
#'   \code{2e-6}, approximately 4 sigdigs.  You can check your exact tolerance
#'   by multiplying this value by \code{.Machine$double.eps}
#' @param pgtol is a double precision variable.
#'
#'     On entry pgtol >= 0 is specified by the user.  The iteration
#'     will stop when:
#'
#'        \code{max(\| proj g_i \| i = 1, ..., n) <= lbfgsPgtol}
#'
#'     where pg_i is the ith component of the projected gradient.
#'
#'     On exit pgtol is unchanged.  This defaults to zero, when the
#'     check is suppressed.
#' @param lmm An integer giving the number of BFGS updates retained in the
#'   "L-BFGS-B" method, It defaults to 7.
#' @param abs.tol Used in Nelder-Mead optimization and PORT optimization.
#'   Absolute tolerance. Defaults to 0 so the absolute convergence test is not
#'   used. If the objective function is known to be non-negative, the previous
#'   default of 1e-20 would be more appropriate.
#' @param xf.tol Used in Nelder-Mead optimization and PORT optimization. false
#'   convergence tolerance. Defaults to 2.2e-14. See \code{\link[stats]{nlminb}}
#'   for more details.
#' @param step.min Used in Nelder-Mead optimization and PORT optimization.
#'   Minimum step size. By default this is 1. See \code{\link[stats]{nlminb}}
#'   for more details.
#' @param step.max Used in Nelder-Mead optimization and PORT optimization.
#'   Maximum step size. By default this is 1. See \code{\link[stats]{nlminb}}
#'   for more details.
#' @param sing.tol Used in Nelder-Mead optimization and PORT optimization.
#'   Singular convergence tolerance; defaults to rel.tol. See
#'   \code{\link[stats]{nlminb}} for more details.
#' @param scale.init Used in Nelder-Mead optimization and PORT optimization.
#'   See \code{\link[stats]{nlminb}} for more details.
#' @param diff.g Used in Nelder-Mead optimization and PORT optimization. An
#'   estimated bound on the relative error in the objective function value. See
#'   \code{\link[stats]{nlminb}} for more details.
#' @param covMethod Method for calculating covariance.  In this discussion, R is
#'   the Hessian matrix of the objective function. The S matrix is the sum of
#'   individual gradient cross-product (evaluated at the individual empirical
#'   Bayes estimates).
#' @param rxControl This uses RxODE family of objects, file, or model
#'   specification to solve a ODE system. See \code{\link[RxODE]{rxControl}} for
#'   more details. By default this is NULL.
#'
#' @inheritParams RxODE::rxSolve
#'
#' @author Mason McComb and Matthew L. Fidler
#' @export
dynmodelControl <- function(...,
                            ci = 0.95,
                            nlmixrOutput = FALSE,
                            digs = 3,
                            lower = -Inf,
                            upper = Inf,
                            ## mma doesn't work
                            ## lbfgsbLG
                            ## slsqp
                            method = c(
                              "bobyqa", "Nelder-Mead", "lbfgsb3c", "L-BFGS-B", "PORT",
                              "mma", "lbfgsbLG", "slsqp", "Rvmmin"
                            ),
                            maxeval = 999,
                            scaleTo = 1.0,
                            scaleObjective = 0,
                            normType = c("rescale2", "constant", "mean", "rescale", "std", "len"),
                            scaleType = c("nlmixr", "norm", "mult", "multAdd"),
                            scaleCmax = 1e5,
                            scaleCmin = 1e-5,
                            scaleC = NULL,
                            scaleC0 = 1e5,
                            # RxODE
                            atol = NULL,
                            rtol = NULL,
                            ssAtol = NULL,
                            ssRtol = NULL,
                            # bobyqaControl
                            npt = NULL,
                            rhobeg = 0.2,
                            rhoend = NULL,
                            iprint = 0,
                            print = 1,
                            maxfun = NULL,
                            # lbfgsb3c
                            trace = 0,
                            factr = NULL,
                            pgtol = NULL,
                            abstol = NULL,
                            reltol = NULL,
                            lmm = NULL,
                            maxit = 100000L,
                            # ,iprint=NULL repreated above
                            # nlminb (PORT)
                            eval.max = NULL,
                            iter.max = NULL,
                            # trace=NULL,
                            abs.tol = NULL,
                            rel.tol = NULL,
                            x.tol = NULL,
                            xf.tol = NULL,
                            step.min = NULL,
                            step.max = NULL,
                            sing.tol = NULL,
                            scale.init = NULL,
                            diff.g = NULL,
                            ## Sigdig
                            boundTol = NULL,
                            epsilon = NULL,
                            derivSwitchTol = NULL,
                            sigdig = 4,
                            covMethod = c("nlmixrHess", "optimHess"),
                            # rxControl
                            gillK = 10L,
                            gillStep = 4,
                            gillFtol = 0,
                            gillRtol = sqrt(.Machine$double.eps),
                            gillKcov = 10L,
                            gillStepCov = 2,
                            gillFtolCov = 0,
                            rxControl = NULL) {
  if (is.null(boundTol)) {
    boundTol <- 5 * 10^(-sigdig + 1)
  }
  if (is.null(epsilon)) {
    epsilon <- 10^(-sigdig - 1)
  }
  if (is.null(abstol)) {
    abstol <- 10^(-sigdig - 1)
  }
  if (is.null(reltol)) {
    reltol <- 10^(-sigdig - 1)
  }

  if (is.null(rhoend)) {
    rhoend <- 10^(-sigdig - 1)
  }
  if (is.null(factr)) {
    factr <- 10^(-sigdig - 1) / .Machine$double.eps
  }
  if (is.null(atol)) {
    atol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(rtol)) {
    rtol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(ssAtol)) {
    ssAtol <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(ssRtol)) {
    ssRtol <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(rel.tol)) {
    rel.tol <- 10^(-sigdig - 1)
  }
  if (is.null(x.tol)) {
    x.tol <- 10^(-sigdig - 1)
  }
  if (is.null(derivSwitchTol)) {
    derivSwitchTol <- 2 * 10^(-sigdig - 1)
  }
  if (inherits(rxControl, "list")) {
    if (!any(names(rxControl) == "atol")) {
      rxControl$atol <- 0.5 * 10^(-sigdig - 2)
    }
    if (!any(names(rxControl) == "rtol")) {
      rxControl$rtol <- 0.5 * 10^(-sigdig - 2)
    }
    if (!any(names(rxControl) == "ssAtol")) {
      rxControl$ssAtol <- 0.5 * 10^(-sigdig - 1.5)
    }
    if (!any(names(rxControl) == "ssRtol")) {
      rxControl$ssRtol <- 0.5 * 10^(-sigdig - 1.5)
    }
    rxControl <- do.call(RxODE::rxControl, rxControl)
  } else {
    atol <- 0.5 * 10^(-sigdig - 2)
    rtol <- 0.5 * 10^(-sigdig - 2)
    ssAtol <- 0.5 * 10^(-sigdig - 1.5)
    ssRtol <- 0.5 * 10^(-sigdig - 1.5)
    rxControl <- RxODE::rxControl(atol = atol, rtol = rtol, ssAtol = ssAtol, ssRtol = ssRtol)
  }

  if (missing(method)) {
    method <- "bobyqa"
  }
  if (missing(normType)) {
    normType <- "rescale2"
  } # normType= 'constant'
  if (missing(scaleType)) {
    scaleType <- "nlmixr"
  } # scaleType= 'norm'

  .ret <- list(
    ci = ci,
    nlmixrOutput = nlmixrOutput,
    digs = digs,
    lower = lower,
    upper = upper,
    method = method,
    maxeval = maxeval,
    scaleTo = scaleTo,
    scaleObjective = scaleObjective,
    normType = normType,
    scaleType = scaleType,
    scaleCmax = scaleCmax,
    scaleCmin = scaleCmin,
    scaleC = scaleC, # NULL,
    scaleC0 = scaleC0,
    # RxODE
    atol = atol,
    rtol = rtol,
    # bobyqa
    npt = npt,
    rhobeg = rhobeg,
    rhoend = rhoend,
    iprint = iprint, # also used in lbfgsb3c
    print = print,
    maxfun = maxfun,
    # lbfgsb3c
    trace = trace,
    factr = factr,
    pgtol = pgtol,
    abstol = abstol,
    reltol = reltol,
    lmm = lmm,
    maxit = maxit,
    # ,iprint = iprint,
    # nlminb (PORT)
    eval.max = eval.max,
    iter.max = iter.max,
    # trace=NULL,
    abs.tol = abs.tol,
    rel.tol = rel.tol,
    x.tol = x.tol,
    xf.tol = xf.tol,
    step.min = step.min,
    step.max = step.max,
    sing.tol = sing.tol,
    scale.init = scale.init,
    diff.g = diff.g,
    covMethod = match.arg(covMethod),
    rxControl = rxControl,
    gillK = as.integer(gillK),
    gillKcov = as.integer(gillKcov),
    gillRtol = as.double(gillRtol),
    gillStep = as.double(gillStep),
    gillStepCov = as.double(gillStepCov),
    gillFtol = as.double(gillFtol),
    gillFtolCov = as.double(gillFtolCov)
  )
  .w <- which(sapply(.ret, is.null))
  .ret <- .ret[-.w]
  class(.ret) <- "dynmodelControl"
  return(.ret)
}


# dynmodel()  #############################################################
#' Fit a non-population dynamic model
#'
#' @param system RxODE object. See \code{\link[RxODE]{RxODE}} for more details.
#' @param model Error model.
#' @param inits Initial values of system parameters.
#' @param data Dataset to estimate. Needs to be RxODE compatible in EVIDs.
#' @param fixPars Fixed system parameters. Default is NULL.
#' @param nlmixrObject nlmixr object. See \code{\link[nlmixr]{nlmixr}} for more
#'   details. Default is NULL.
#' @param control Control options for dynmodel
#'   \code{\link[nlmixr]{dynmodelControl}} .
#' @param ... Other parameters (ignored)
#' @author Wenping Wang, Mason McComb and Matt Fidler
#' @examples
#'
#' \donttest{
#'
#' # dynmodel example --------------------------------------------------------
#' ode <- "
#'       kel = CL/V;
#'       d/dt(X) = -kel*X;
#'       C=X/V;
#'       PRED = C
#'       "
#' ode_system <- RxODE(model = ode)
#' model_error_structure <- cp ~ C + add(0.01) + prop(0.01)
#' inits <- c(CL = 1, V = 10)
#' control <- dynmodelControl(method = "Nelder-Mead")
#' fit <-
#'   dynmodel(
#'     system = ode_system,
#'     model = model_error_structure,
#'     data = Bolus_1CPT,
#'     inits = inits,
#'     control = control
#'   )
#'
#' # nlmixr model example ----------------------------------------------------------
#' model_onecmt_bolus <- function() {
#'   ini({
#'     CL <- c(0, 5, 10) # Clearance (L/hr)
#'     V <- c(0, 50, 100) # Volume of Distribution
#'     prop.err <- c(0, 0.01, 1)
#'   })
#'   model({
#'     kel <- CL / V
#'     d / dt(X) <- -kel * X
#'     cp <- X / V
#'     cp ~ prop(prop.err)
#'   })
#' }
#'
#' fit <- nlmixr(object = model_onecmt_bolus, data = Bolus_1CPT, est = "dynmodel")
#'
#' as.dynmodel(fit)
#'
#' # method = "focei" is slightly more flexible and well tested
#'
#' fit <- nlmixr(object = model_onecmt_bolus, data = Bolus_1CPT, est = "focei")
#'
#' }
#' @export
dynmodel <- function(system, model, inits, data, fixPars = NULL, nlmixrObject = NULL, control = list(), ...) {
  # Timing and environment --------------------------------------------------
  .pt <- proc.time()
  .time <- c()
  .norm <- .dnorm <- .logn <- .dlnorm <- NULL

  .dynmodel.env <- new.env(parent = emptyenv())

  # dynmodelControl Handling ------------------------------------------------
  if (!RxODE::rxIs(control, "dynmodelControl")) {
    control <- do.call(dynmodelControl, control)
  }

  # reassign control names
  for (i in 1:length(control)) {
    assign(names(control[i]), control[[i]])
  }

  # Error model  -------------------------------------------------------------
  modelList <-
    if (class(model) == "formula") {
      list(model)
    } else {
      model
    }
  inits.err <- NULL
  .env <- environment()

  modelParsed <-
    lapply(modelList, function(.model) {
      .model <- unlist(lapply(attr(terms(.model), "variables"), as.list))
      .model <- sapply(.model, deparse)

      # assign error terms
      .sigma.add <-
        if ("add" %in% .model) {
          as.numeric(.model[which(.model == "add") + 1])
        } else {
          NULL
        }
      .sigma.prop <-
        if ("prop" %in% .model) {
          as.numeric(.model[which(.model == "prop") + 1])
        } else {
          NULL
        }
      .sigma.pow <-
        if ("pow" %in% .model) {
          as.numeric(.model[which(.model == "pow") + 1])
        } else {
          NULL
        }
      .sigma.pow2 <-
        if ("pow2" %in% .model) {
          as.numeric(.model[which(.model == "pow2") + 1])
        } else {
          NULL
        }
      .sigma.yeoJohnson <-
        if ("yeoJohnson" %in% .model) {
          as.numeric(.model[which(.model == "yeoJohnson") + 1])
        } else {
          NULL
        }
      .sigma.boxCox <-
        if ("boxCox" %in% .model) {
          as.numeric(.model[which(.model == "boxCox") + 1])
        } else {
          NULL
        }
      .sigma.norm <-
        if ("norm" %in% .model) {
          assign(".norm", TRUE, .env)
          as.numeric(.model[which(.model == "norm") + 1])
        } else {
          assign(".norm", NULL, .env)
          NULL
        }
      .sigma.dnorm <-
        if ("dnorm" %in% .model) {
          assign(".dnorm", TRUE, .env)
          as.numeric(.model[which(.model == "dnorm") + 1])
        } else {
          assign(".dnorm", NULL, .env)
          NULL
        }
      .sigma.logn <-
        if ("logn" %in% .model) {
          assign(".logn", TRUE, .env)
          as.numeric(.model[which(.model == "logn") + 1])
        } else {
          assign(".logn", NULL, .env)
          NULL
        }
      .sigma.dlnorm <-
        if ("dlnorm" %in% .model) {
          assign(".dlnorm", TRUE, .env)
          as.numeric(.model[which(.model == "dlnorm") + 1])
        } else {
          assign(".dlnorm", NULL, .env)
          NULL
        }
      .sigma.tbs <-
        if ("tbs" %in% .model) {
          as.numeric(.model[which(.model == "tbs") + 1])
        } else {
          NULL
        }
      .sigma.tbsYj <-
        if ("tbsYj" %in% .model) {
          as.numeric(.model[which(.model == "tbsYj") + 1])
        } else {
          NULL
        }
      # keep error model terms
      inits.err <-
        c(
          add = .sigma.add,
          prop = .sigma.prop,
          pow = .sigma.pow,
          pow2 = .sigma.pow2,
          yeoJohnson = .sigma.yeoJohnson,
          boxCox = .sigma.boxCox,
          norm = .sigma.norm,
          dnorm = .sigma.dnorm,
          tbs = .sigma.tbs,
          tbsYj = .sigma.tbsYj
        )
      inits.err <- inits.err[which(names(inits.err) %in% intersect(names(inits.err), .model))]
      assign("inits.err", inits.err, .env)
      .model <- c("dv" = .model[2], "pred" = .model[3], inits.err)
    })

  inits <- c(inits, inits.err)
  if ("pow2" %in% names(inits) & !("pow" %in% names(inits))) {
    stop("Error Model: pow must be defined when using pow2")
  }

  # Check dynmodel() inputs, Define vars, modelVars, pars,  ------------

  # NOTES: Check to make sure all there is consistency between error model, data. inits, and ODE model

  # data handling
  ## data <- RxODE::etTrans(data, system, addCmt = TRUE, dropUnits = TRUE, allTimeVar = TRUE)
  .original.data <- data

  # warn DV must be in data
  if (!("DV" %in% names(data))) {
    stop("data must contain column named DV")
  }
  # remove row that has time zero with zero observation, and no event
  if (data$TIME[1] == 0 & data$EVID[1] == 0 & data$DV[1] == 0) {
    rows <- -1
  } else {
    rows <- TRUE
  }
  # add model variable to data
  .dv.name <- unlist(modelParsed)["dv"][[1]]
  if (any(names(data) %in% .dv.name) == FALSE) {
    data[[.dv.name]] <- data$DV
  }
  if (any(data$EVID %in% 2)) {
    ## warning("Removing EVID==2 rows from the data; they are not supported in dynmodel.")
    data <- data[data$EVID != 2, ]
    if (nrow(data) == 0) {
      stop("Zero rows of data remain after removing EVID == 2.")
    }
  }
  data0 <- data[data$EVID == 0, ]
  nobst <- length(data0$DV)

  # Error "model" contains "data" variables?
  # check to see if there is a discrepency between error model names and data
  nodef <- setdiff(sapply(modelParsed, function(x) x["dv"]), names(data))
  # print error message
  if (length(nodef) & is.null(nlmixrObject)) {
    msg <- err.msg(nodef, pre = "var(s) not found in data: ")
    stop(msg)
  }

  # "system" variables contain error "model" variables?
  # obtain all variables from the system
  modelVars <- system$cmpMgr$get.modelVars()
  # reassign vars to combine state and lhs variables
  vars <- c(modelVars$state, modelVars$lhs)
  # Check to see if the prediction term is in the error model
  nodef <- setdiff(sapply(modelParsed, function(x) x["pred"]), vars)
  # print error message
  if (length(nodef)) {
    msg <- err.msg(nodef, pre = "modelVar(s) not found in model: ")
    stop(msg)
  }

  #  "system" variables contain estimated "init" variables and fixed "fixPars" variables?
  # obtain fixed and estimated parameters
  pars <-
    if (is.null(nlmixrObject)) {
      modelVars$params
    } else {
      names(
        nlmixrDynmodelConvert(nlmixrObject)$inits
      )
    }

  # Check to see if there are values in pars, that are not in the initial conditions and fixed parameters
  # nodef = setdiff(pars, c(names(inits), names(fixPars)))
  # # print error message
  # if (length(nodef)) {
  #  msg = err.msg(nodef, pre="par(s) not found: ")
  #  stop(msg)
  # }

  # Additional assignment ---------------------------------------------------
  # number of estimated parameters, excluding the error terms
  npar <- length(pars) - length(fixPars)

  # Objective Function ------------------------------------------------------
  .time$setupTime <- (proc.time() - .pt)["elapsed"]

  .funs <- list()
  sgy <- c()
  nobs <- c()

  rxControl <- control$rxControl

  reducedTol <- FALSE

  obj <- function(th) {
    # unscale
    unscaled.th <- numeric(length(th))
    th <- sapply(seq_along(th), function(x) {
      unscalePar(th, x)
    })
    .funs$unscaled(th)

    # define parameters used for simulation, all parameters except the error terms
    .ixpar <- npar
    theta <- th[1:npar]
    names(theta) <- names(inits)[1:npar]
    theta <- c(theta, fixPars)

    if (!is.null(nlmixrObject)) {
      theta <- nlmixrObject$dynmodel.fun(theta)
    }
    .rxControl <- rxControl
    .rxControl$returnType <- "data.frame"
    # call rxODE for simulation
    s <-
      do.call(
        RxODE::rxSolve,
        c(
          list(object = system, params = theta, events = data),
          .rxControl
        )
      )
    i <- 1
    i.max <- 10
    while (any(is.na(s$nlmixr_pred)) & i < i.max) {
      .rxControl$atol <- .rxControl$atol * 100
      .rxControl$rtol <- .rxControl$rtol * 100

      s <-
        do.call(
          RxODE::rxSolve,
          c(
            list(object = system, params = theta, events = data),
            .rxControl
          )
        )
      assign("reducedTol", TRUE, .env)
      i <- i + 1
    }

    # sum of log-likelihood function:
    l <- lapply(modelParsed, function(x) {
      # name the inputs
      names(th) <- names(inits)

      # initialize sigmas for objective function
      add <- 0
      prop <- 0
      pow <- 1
      pow2 <- 1
      norm <- NULL
      dnorm <- NULL
      boxCox <- NULL
      tbs <- NULL
      yeoJohnson <- NULL
      tbsYj <- NULL

      # assign names sigma names
      if (any(names(th) %in% names(modelParsed[[1]]))) {
        for (i in 1:sum((names(th) %in% names(modelParsed[[1]])))) {
          assign(
            names(th[names(th) %in% names(modelParsed[[1]])])[i],
            as.numeric(th[names(th) %in% names(modelParsed[[1]])])[i]
          )
        }
      }

      if (!is.null(.norm)) add <- norm
      if (!is.null(.dnorm)) add <- dnorm
      if (!is.null(.logn)) {
        lambda <- 0
        th <- th[names(th) != "logn"]
      }
      if (!is.null(dlnorm)) {
        lambda <- 0
        th <- th[names(th) != "dlnorm"]
      }
      if (!is.null(boxCox)) lambda <- boxCox
      if (!is.null(tbs)) lambda <- tbs
      if (!is.null(yeoJohnson)) lambda <- yeoJohnson
      if (!is.null(tbsYj)) lambda <- tbsYj

      # predicted and observed values from RxODE
      yo <- data0[, "DV"]
      yp <- s[, x["pred"]]
      assign("nobs", length(yo), .env)
      if (length(yo) != length(yp)) {
        stop("length mismatch ", length(yo), "!=", length(yp))
      }

      # log normal transformation ----
      if (!is.null(.logn) | !is.null(.dlnorm)) {
        .h.x <- boxCox(yo, lambda) # log(yo) # obs
        .h.y <- boxCox(yp, lambda) # log(yp)  # pred
        if ("pow" %in% names(modelParsed[[1]])) {
          .h.y.var <- yp^(2 * pow2) * thresh(pow)^2 + thresh(add)^2 # variance of pred
        } else {
          .h.y.var <- yp^(2 * pow2) * thresh(prop)^2 + thresh(add)^2 # variance of pred
        }
        # boxCox transformed -2 log-likelihood
        .logn.n2ll <- log(.h.y.var) + ((.h.x - .h.y)^2) / .h.y.var
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll <- .logn.n2ll - 2 * (0 - 1) * log(yo) - 2 * log(2 * pi) # lambda is zero here
        assign("sgy", .h.y.var, .env)
        # negative log-likelihood function for output
        ll <- .5 * (.n2ll)
      }
      # boxCox Transform ----
      else if ("boxCox" %in% names(modelParsed[[1]]) | "tbs" %in% names(modelParsed[[1]])) {
        .h.x <- boxCox(yo, lambda) # obs
        .h.y <- boxCox(yp, lambda) # pred
        if ("pow" %in% names(modelParsed[[1]])) {
          .h.y.var <- (yp^(2 * pow2)) * thresh(pow)^2 + thresh(add)^2 # variance of pred
        } else {
          .h.y.var <- yp^(2 * pow2) * thresh(prop)^2 + thresh(add)^2 # variance of pred
        }
        # boxCox transformed -2 log-likelihood
        .boxCox.n2ll <- log(.h.y.var) + ((.h.x - .h.y)^2) / .h.y.var
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll <- .boxCox.n2ll - 2 * (lambda - 1) * log(yo) - 2 * log(2 * pi)
        assign("sgy", .h.y.var, .env)
        # negative log-likelihood function for output
        ll <- .5 * (.n2ll)
      }
      # yeoJohnson Transform ----
      else if ("yeoJohnson" %in% names(modelParsed[[1]]) | "tbsYj" %in% names(modelParsed[[1]])) {
        .h.x <- yeoJohnson(yo, lambda) # obs
        .h.y <- yeoJohnson(yp, lambda) # pred
        if ("pow" %in% names(modelParsed[[1]])) {
          .h.y.var <- (yp^(2 * pow2)) * thresh(pow)^2 + thresh(add)^2 # variance of pred
        } else {
          .h.y.var <- yp^(2 * pow2) * thresh(prop)^2 + thresh(add)^2 # variance of pred
        }
        # yeoJohnson transformed -2 log-likelihood
        .yeoJohnson.n2ll <- log(.h.y.var) + ((.h.x - .h.y)^2) / .h.y.var
        # back-transformed  -2 log-likelihood function, with penalty added
        .n2ll <-
          ifelse(
            yo >= 0,
            .yeoJohnson.n2ll - 2 * (lambda - 1) * log(yo + 1) - 2 * log(2 * pi),
            .yeoJohnson.n2ll - 2 * (1 - lambda) * log(-yo + 1) - 2 * log(2 * pi)
          )
        assign("sgy", .h.y.var, .env)
        # negative log-likelihood function for output
        ll <- .5 * (.n2ll)
      }
      # power model ----
      else if ("pow2" %in% names(modelParsed[[1]])) {
        sgy <- thresh(add) + thresh(pow) * yp^(pow2)
        assign("sgy", sgy, envir = .dynmodel.env)
        ll <- .5 * ((yo - yp)^2 / (sgy^2) + log(sgy^2) + log(2 * pi))
      }
      # all other error models ----
      else {
        #  if (identical(c("dv","pred"),names(model[[1]]))){
        if (length(names(modelParsed[[1]])) == 2) {
          sgy <- 1
          assign("sgy", sgy, envir = .dynmodel.env)
        } else {
          sgy <- thresh(add) + thresh(prop) * yp
          assign("sgy", sgy, envir = .dynmodel.env)
        }
        ll <- .5 * ((yo - yp)^2 / (sgy^2) + log(sgy^2) + log(2 * pi))
      }
      sum(ll, na.rm = TRUE)
    })
    assign("sgy", sgy, envir = .dynmodel.env)
    do.call("sum", l) # same as return(as.numeric(l)), l is a list for each value in the model?
  }

  # FIXME: Put options from control here gillK etc
  .funs <-
    nlmixrGradFun(
      obj,
      print = control$print,
      gillRtol = control$gillRtol,
      gillK = control$gillK,
      gillStep = control$gillStep,
      gillFtol = control$gillFtol,
      thetaNames = names(inits)
    )

  # Scaling functions -----------------------------------------------------------------------
  # normType assignment for scaling (normalization type)
  normType <- control$normType
  if (is.null(nlmixrObject)) {
    normType <- "constant"
    scaleType <- "norm"
  }
  if (normType == "constant") {
    C1 <- 0
    C2 <- 1
  } else if (normType == "rescale2") {
    C1 <- (max(inits) + min(inits)) / 2
    C2 <- (max(inits) - min(inits)) / 2
  } else if (normType == "mean") {
    C1 <- mean(inits)
    C2 <- max(inits) - min(inits)
  } else if (normType == "rescale") {
    C1 <- min(inits)
    C2 <- max(inits) - min(inits)
  } else if (normType == "std") {
    C1 <- mean(inits)
    C2 <- sd(inits)
  } else if (normType == "len") {
    C1 <- 0
    C2 <- sqrt(sum(inits * inits))
  }

  # produces vector of scaleC values if missing, handles incorrect length of scaleC values.
  scaleC <- control$scaleC
  if (is.null(scaleC)) {
    scaleC <- rep(1, length(inits))
  } else {
    scaleC <- as.double(scaleC)
    if (length(inits) > length(scaleC)) {
      scaleC <- c(scaleC, rep(1, length(inits) - length(scaleC)))
    } else if (length(inits) < length(scaleC)) {
      scaleC <- scaleC[seq(1, length(inits))]
      warning("scaleC control option has more options than estimated parameters, please check.")
    }
  }

  # assign value to scaleC. If log value, scaleC <- 1, else scalec <- 1/abs(inits[i])
  # need to known when dynmodel has log values, and which ones, or else there
  # will be a problem scaling.
  theta <- inits[seq_len(length(inits) - length(inits.err))]
  names(scaleC) <- names(inits)

  if (!is.null(nlmixrObject)) {
    .model <-
      RxODE::rxSymPySetupPred(
        nlmixrObject$rxode.pred,
        function() {
          return(nlmixr_pred)
        },
        nlmixrObject$theta.pars,
        nlmixrObject$error,
        grad = FALSE,
        pred.minus.dv = TRUE,
        sum.prod = FALSE, # control$sumProd,
        theta.derivs = FALSE,
        optExpression = TRUE, # control$optExpression,
        run.internal = TRUE,
        only.numeric = TRUE
      )

    scaleC[.model$log.thetas] <- 1
    scaleC[setdiff(1:length(theta), .model$log.thetas)] <-
      1 / abs(theta[setdiff(1:length(theta), .model$log.thetas)])
  }

  # assign value to scaleC based on the error model used
  n.inits.err <- names(inits.err)
  for (i in seq_along(n.inits.err)) {
    if (is.na(scaleC[n.inits.err[i]])) {
      if (n.inits.err[i] %in% c("boxCox", "yeoJohnson", "pow2", "tbs", "tbsYj")) {
        scaleC[n.inits.err[i]] <- 1
      } else if (n.inits.err[i] %in% c("prop", "add", "norm", "dnorm", "logn", "dlogn", "lnorm", "dlnorm")) {
        scaleC[n.inits.err[i]] <- 0.5 * abs(inits.err[i])
      }
    }
  }

  # assign value to scaleC based on additional functions
  if (!is.null(nlmixrObject)) {
    for (i in .model$extraProps$powTheta) {
      if (is.na(scaleC[i])) scaleC[i] <- 1 ## Powers are log-scaled
    }

    .ini <- as.data.frame(nlmixrObject$ini)

    for (i in .model$extraProps$factorial) {
      if (is.na(scaleC[i])) scaleC[i] <- abs(1 / digamma(.ini$est[i] + 1))
    }

    for (i in .model$extraProps$gamma) {
      if (is.na(scaleC[i])) scaleC[i] <- abs(1 / digamma(.ini$est[i]))
    }

    for (i in .model$extraProps$log) {
      if (is.na(scaleC[i])) scaleC[i] <- log(abs(.ini$est[i])) * abs(.ini$est[i])
    }
  }
  ## FIXME: needs to be based on actual initial values in sin because typically change to correct scale
  ## Ctime is also usually used for circadian rhythm models
  ## for (.i in .ini$model$extraProps$sin){
  ##     if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- fabs(tan(.ini$est[.i]));
  ## }
  ## for (.i in .ini$model$extraProps$cos){
  ##     if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- fabs(1 / tan(.ini$est[.i]));
  ## }
  ## for (.i in .ini$model$extraProps$tan){
  ##     if (is.na(.ret$scaleC[.i])) .ret$scaleC[.i] <- fabs(2 * sin(2 * .ini$est[.i]));
  ## }

  # Function for scaling parameters based on scaleType
  scalePar <- function(x, i) {
    if (scaleType == "norm") { # simple scaling
      return((x[i] - C1) / C2)
    } else if (scaleType == "nlmixr") { # nlmixr
      scaleTo <- (inits[i] - C1) / C2
      return((x[i] - inits[i]) / scaleC[i] + scaleTo)
    } else if (scaleType == "mult") { # simple multiplicatice scaling
      if (scaleTo > 0) {
        return(x[i] / inits[i] * scaleTo)
      } else {
        return(x[i])
      }
    } else if (scaleType == "multAdd") { # log non-log multiplicative scaling
      if (scaleTo > 0) {
        return((x[i] - inits[i]) + scaleTo)
      } else {
        return(x[i] / inits[i] * scaleTo)
      }
    } else { # When should this be used? "norm" scaling is essentially no scaling when normType specified.
      if (scaleTo > 0) {
        return((x[i] - inits[i]) + scaleTo)
      } else {
        return(x[i])
      }
    }
  }

  # Function for unscaling parameters based on scaleType
  unscalePar <- function(x, i) {
    if (scaleType == "norm") { # simple scaling
      return(x[i] * C2 + C1)
    } else if (scaleType == "nlmixr") { # nlmixr
      scaleTo <- (inits[i] - C1) / C2
      return((x[i] - scaleTo) * scaleC[i] + inits[i])
    } else if (scaleType == "mult") { # simple multiplicatice scaling
      if (scaleTo > 0) {
        return(x[i] * inits[i] / scaleTo)
      } else {
        return(x[i])
      }
    } else if (scaleType == "multAdd") { # log non-log multiplicative scaling
      if (scaleTo > 0) {
        return((x[i] - scaleTo) + inits[i])
      } else {
        return(x[i] * inits[i] / scaleTo)
      }
    } else { # When should this be used? "norm" scaling is essentially no scaling when normType specified.
      if (scaleTo > 0) {
        return((x[i] - scaleTo) * 1 + inits[i])
      } else {
        return(x[i])
      }
    }
  }

  # Scale --------------
  scaleType <- control$scaleType
  lower <- control$lower
  upper <- control$upper
  if (normType != "constant" & scaleType != "norm") {
    .st <- proc.time()
    # scaled the initial conditions
    inits.temp <- numeric(length(inits))
    for (i in 1:length(inits)) {
      inits.temp[i] <- scalePar(inits, i)
    }
    .inits <- inits.temp
    # scaled lower boundary
    if (!is.null(lower)) {
      if (any(lower == -Inf)) {
        lower[which(lower == -Inf)] <- 0
        warning("Lower boundary of -Inf set to 0.")
      }
      lower.temp <- numeric(length(lower))
      for (i in 1:length(lower)) {
        lower.temp[i] <- scalePar(lower, i)
      }
      .lower <- lower.temp
    } else {
      .lower <- NULL
    }

    # scaled upper boundary
    if (!is.null(upper)) {
      upper.temp <- numeric(length(upper))
      for (i in 1:length(upper)) {
        upper.temp[i] <- scalePar(upper, i)
      }
      .upper <- upper.temp
    } else {
      .upper <- NULL
    }
    .time$scalingTime <- (proc.time() - .st)["elapsed"]
  } else {
    .st <- proc.time()
    .inits <- inits
    .lower <- lower
    .upper <- upper
    .time$scalingTime <- (proc.time() - .st)["elapsed"]
  }
  method <- control$method
  # Optimization -----------------------------------------------------------------------
  if (method == "bobyqa") {
    .optFun <- .bobyqa
  } else if (method %in% c("nlminb", "PORT")) {
    .optFun <- .nlminb
  } else if (method == "mma") {
    .optFun <- .nloptr
  } else if (method == "slsqp") {
    .optFun <- .slsqp
  } else if (method == "lbfgsbLG") {
    .optFun <- .lbfgsbLG
  } else if (method == "Rvmmin") {
    .optFun <- .Rvmmin
  } else if (method == "Nelder-Mead") {
    .optFun <- .mymin
  } else if (method == "lbfgsb3c") {
    .optFun <- .lbfgsb3c
  } else if (method == "L-BFGS-B") {
    .optFun <- .lbfgsbO
  } else {
    stop("Optimization method unknown.")
  }
  .ot <- proc.time()

  fit <- .optFun(as.vector(.inits), fn = .funs$eval, gr = .funs$gr, lower = .lower, upper = .upper, control = control)
  .message <- fit$message
  ## fit$value
  assign("fit", fit, envir = .dynmodel.env)

  .time$optimizationTime <- (proc.time() - .ot)["elapsed"]

  # Hessian -----------------------------------------------------------------------
  .ht <- proc.time()

  fit$unscaledPar <- fit$par
  fit$par <- sapply(seq_along(fit$par), function(x) {
    unscalePar(fit$par, x)
  })

  fit$normType <- normType
  fit$scaleType <- scaleType

  normType <- "constant"
  C1 <- 0
  C2 <- 1
  scaleType <- "norm"

  if (control$covMethod == "optimHess") {
    fit$hessian <- try(optimHess(fit$par, obj, control = control), silent = TRUE)
  } else {
    # FIXME: Put options from control here gillK etc
    fit$hessian <- try(nlmixrHess(fit$par, obj,
      gillRtol = control$gillRtol,
      gillK = control$gillKcov,
      gillStep = control$gillStepCov,
      gillFtol = control$gillFtolCov
    ), silent = TRUE)
  }


  if (inherits(fit$hessian, "try-error")) {
    se <- rep(NA, length(fit$par))
    warning("standard error of the Hessian has failed")
  } else {
    cov.matrix <- solve(fit$hessian)
    se <- sqrt(diag(solve(fit$hessian)))
  }

  # reassign the negative values to positive for add, prop/pow since they are standard deviations
  if (!is.na(match("add", names(inits)))) {
    fit$par[match("add", names(inits))] <- abs(fit$par[match("add", names(inits))])
  }
  if (!is.na(match("prop", names(inits)))) {
    fit$par[match("prop", names(inits))] <- abs(fit$par[match("prop", names(inits))])
  }
  if (!is.na(match("pow", names(inits)))) {
    fit$par[match("pow", names(inits))] <- abs(fit$par[match("pow", names(inits))])
  }
  if (!is.na(match("norm", names(inits)))) {
    fit$par[match("norm", names(inits))] <- abs(fit$par[match("norm", names(inits))])
  }
  if (!is.na(match("dnorm", names(inits)))) {
    fit$par[match("dnorm", names(inits))] <- abs(fit$par[match("dnorm", names(inits))])
  }

  .time$hessianTime <- (proc.time() - .ht)["elapsed"]
  # dynmodel Output -------------------------------------------------------
  # unscale optimized parameters here if scaling was used:

  # create table for output
  res <- cbind(fit$par, abs(se), abs(se / fit$par * 100))

  dimnames(res) <- list(names(inits), c("est", "se", "%cv"))

  if (!is.null(fit$objective)) fit$value <- fit$objective

  # Output
  res <- c(list(res = res, obj = obj, npar = length(fit$par), nobs = nobs, data = data), fit)

  class(res) <- "dyn.ID"
  # Final Output ----------------------------------------------------------
  .time$totalTime <- (proc.time() - .pt)["elapsed"]

  if (reducedTol) warning("atol and rtol reduced to complete optimization")

  .time <- as.data.frame(.time)
  names(.time) <- c("setup", "scaling", "optimization", "Hessian", "total")

  if (!is.null(nlmixrObject) & control$nlmixrOutput) {
    nlmixr.ouptut <- as.focei.dynmodel(
      .dynmodelObject = res, .nlmixrObject = nlmixrObject, .data = .original.data, .time = .time, .fit = fit, .message = .message, .inits.err = inits.err, .cov = cov.matrix, .sgy = sgy,
      .dynmodelControl = control, .nobs2 = 0, .pt = proc.time(), .rxControl = rxControl
    )
    ## .hist <- .funs$hist()
    ## assign("parHistData", .hist, nlmixr.ouptut$env)
    ## .tmp <- .hist
    ## .tmp <- .tmp[.tmp$type == "Unscaled", names(.tmp) != "type"]
    ## .iter <- .tmp$iter
    ## .tmp <- .tmp[, names(.tmp) != "iter"]
    ## .ret <- nlmixr.ouptut$env
    ## .ret$parHistStacked <- data.frame(stack(.tmp), iter = .iter)
    ## names(.ret$parHistStacked) <- c("val", "par", "iter")
    ## .ret$parHist <- data.frame(iter = .iter, .tmp)
    return(nlmixr.ouptut)
  }
  else {
    return(res)
  }
}

# ####################################################################### #
#
## MCMC Section
#
# ####################################################################### #
uni_slice <- function(x0, fr, rho = NULL, w = 1, m = 1000, lower = -1.0e20, upper = 1.0e20) {
  if (is.null(rho)) rho <- environment(fr)
  .Call(slice_wrap, fr, rho, x0, w, as.integer(m), lower, upper, PACKAGE = "nlmixr")$x1
}

# Error model  -------------------------------------------------------------
genobj <- function(system, model, evTable, inits, data, fixPars = NULL,
                   squared = TRUE) {

  # Error model  -------------------------------------------------------------
  if (class(model) == "formula") {
    model <- list(model)
  }
  inits.err <- NULL
  .env <- environment()
  model <- lapply(model, function(f) {
    s <- unlist(lapply(attr(terms(f), "variables"), as.list))
    s <- sapply(s, deparse)

    ix.add <- match("add", s, nomatch = 0)
    ix.pro <- match("prop", s, nomatch = 0)
    err.type <- c("add", "prop", "combo")[(ix.add > 0) + 2 * (ix.pro > 0)]

    sig.add <- if (ix.add > 0) as.numeric(s[ix.add + 1]) else NULL
    sig.pro <- if (ix.pro > 0) as.numeric(s[ix.pro + 1]) else NULL

    assign("inits.err", c(inits.err, sig.add, sig.pro), envir=.env)

    if (any(is.na(inits.err) | inits.err <= 0)) {
      stop("error model misspecification")
    }

    s <- c(s[2:3], err.type)
    names(s) <- c("dv", "pred", "err")
    s
  })
  names(inits.err) <- paste0("err", 1:length(inits.err))
  inits <- c(inits, inits.err)

  # Check dynmodel() inputs, Define vars, modelVars, pars,  ------------
  vars <- names(data)
  nodef <- setdiff(sapply(model, function(x) x["dv"]), vars)
  if (length(nodef)) {
    msg <- err.msg(nodef, pre = "var(s) not found in data: ")
    stop(msg)
  }

  modelVars <- system$cmpMgr$get.modelVars()
  vars <- c(modelVars$state, modelVars$lhs)
  nodef <- setdiff(sapply(model, function(x) x["pred"]), vars)
  if (length(nodef)) {
    msg <- err.msg(nodef, pre = "var(s) not found in model: ")
    stop(msg)
  }

  pars <- modelVars$params
  nodef <- setdiff(pars, c(names(inits), names(fixPars)))
  if (length(nodef)) {
    msg <- err.msg(nodef, pre = "par(s) not found: ")
    stop(msg)
  }

  npar <- length(pars) - length(fixPars)

  # Additional assignment ---------------------------------------------------
  ## is this necessary ##
  have_zero <- min(data$time) <= 0
  rows <- if (have_zero) TRUE else -1 # used in line 304 in obj()
  ## ---------------- ##
  s.save <- NULL
  .env <- environment()

  # Objective Function ------------------------------------------------------
  obj <- function(th, do.ode.solving = TRUE, negation = FALSE) {
    .ixpar <- npar
    theta <- th[1:npar]
    names(theta) <- names(inits)[1:npar]
    theta <- c(theta, fixPars)
    if (do.ode.solving) {
      s <- system$solve(theta, evTable, atol = 1e-06, rtol = 1e-06)
      assign("s.save", s, .env)
    } else {
      s <- s.save
    }

    l <- lapply(model, function(x) {
      err.combo <- (x["err"] == "combo") + 0
      .ixpar <<- .ixpar + 1
      sig <- th[.ixpar:(.ixpar + err.combo)]
      sig <- if (x["err"] == "add") {
        c(sig, 0)
      } else if (x["err"] == "prop") {
        c(0, sig)
      } else {
        .ixpar <<- .ixpar + 1
        sig
      }

      yp <- s[rows, x["pred"]]
      sgy <- thresh(sig[1] + yp * sig[2])
      yo <- data[, x["dv"]]
      ll <- .5 * ((yo - yp)^2 / sgy^2 + 2 * log(sgy) + log(2 * pi))
      sum(ll)
    })

    res <- do.call("sum", l)
    if (negation) -res else res
  }
  list(obj = obj, inits = inits)
}

#-- mcmc
error.terms <- paste0("err", 1:40)

# slice sampling  ---------------------------------------------------------
do.slice <- function(pars, fr0) {
  rho <- environment()
  lapply(
    names(pars),
    function(wh, fr0) {
      do.ode.solving <- match(wh, error.terms, nomatch = 0) == 0
      pars.cp <- get("pars", rho)
      x0 <- pars.cp[wh]
      fr <- function(x) {
        pars.cp[wh] <- x
        fr0(pars.cp, do.ode.solving = do.ode.solving, negation = TRUE)
      }

      # pars.cp as a data frame, run in do all the ode solving at the end in
      # parallel
      pars.cp[wh] <- uni_slice(x0, fr, lower = 0)
      assign("pars", pars.cp, rho)
      NULL
    },
    fr0 = fr0
  )

  pars
}

# dynmodel.mcmc() ---------------------------------------------------------
#' Fit a non-population dynamic model using mcmc
#'
#' @param system an RxODE object
#' @param model a list of statistical measurement models
#' @param evTable an Event Table object
#' @param inits initial values of system parameters
#' @param data input data
#' @param fixPars fixed system parameters
#' @param nsim number of mcmc iterations
#' @param squared if parameters be squared during estimation
#' @param seed random number seed
#' @author Wenping Wang
#' @return NULL
#' @examples
#' \donttest{
#'
#' ode <- "
#'    dose=200;
#'    pi = 3.1415926535897931;
#'
#'    if (t<=0) {
#'       fI = 0;
#'    } else {
#'       fI = F*dose*sqrt(MIT/(2.0*pi*CVI2*t^3))*exp(-(t-MIT)^2/(2.0*CVI2*MIT*t));
#'    }
#'
#'    C2 = centr/V2;
#'    C3 = peri/V3;
#'    d/dt(centr) = fI - CL*C2 - Q*C2 + Q*C3;
#'    d/dt(peri)  =              Q*C2 - Q*C3;
#' "
#' sys1 <- RxODE(model = ode)
#'
#'
#' ## ------------------------------------------------------------------------
#' dat <- invgaussian
#' mod <- cp ~ C2 + prop(.1)
#' inits <- c(MIT = 190, CVI2 = .65, F = .92)
#' fixPars <- c(CL = .0793, V2 = .64, Q = .292, V3 = 9.63)
#' ev <- eventTable()
#' ev$add.sampling(c(0, dat$time))
#' (fit <- dynmodel.mcmc(sys1, mod, ev, inits, dat, fixPars))
#' }
#' @export
dynmodel.mcmc <- function(system, model, evTable, inits, data,
                          fixPars = NULL, nsim = 500, squared = TRUE,
                          seed = NULL) {
  calls <- match.call()

  # Objective Function ------------------------------------------------------
  l <- genobj(system, model, evTable, inits, data, fixPars, squared)
  rho <- environment()
  pars <- l$inits
  fr0 <- l$obj

  if (is.null(seed)) seed <- 99
  set.seed(seed)

  # progress
  on.exit(RxODE::rxProgressAbort("Aborted MCMC Calculation"))
  RxODE::rxProgress(nsim)

  # slice sampling
  s <-
    t(sapply(
      1:nsim,
      function(k, rho) {
        pars <- do.slice(get("pars", rho), fr0)
        RxODE::rxTick()
        assign("pars", pars, rho)
      },
      rho = rho
    ))

  if (squared) s <- s * s
  attr(s, "calls") <- calls
  attr(s, "obj") <- fr0
  attr(s, "class") <- "dyn.mcmc"
  RxODE::rxProgressStop()
  s
}

#' @rdname print.dyn.mcmc
#' @export
summary.dyn.mcmc <- function(object, ...) {
  print.dyn.mcmc(x = object, ...)
}

#' Print summary of a non-population dynamic model fit using mcmc
#'
#' @param x,object a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
print.dyn.mcmc <- function(x, ...) {
  s <- t(apply(x, 2, function(x) c(mean(x), sd(x), sd(x) / mean(x) * 100)))
  dimnames(s)[[2]] <- c("mean", "sd", "cv%")
  print(s)
  cat("\n# samples:", dim(x)[1], "\n")
}

#' Plot of a non-population dynamic model fit using mcmc
#'
#' Plot of a non-population dynamic model fit using mcmc
#'
#' @param x a dynmodel fit object
#' @param ... additional arguments
#' @return NULL
#' @export
plot.dyn.mcmc <- function(x, ...) {
  fit <- list(obj = attr(x, "obj"), par = apply(x, 2, mean))
  gof(fit)
}
