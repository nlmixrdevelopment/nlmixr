## gnlmm.R: population PK/PD modeling library
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
## along with nlmixr.  If not, see <http:##www.gnu.org/licenses/>.

# to suppress Rcheck warning

parseOM <- function(OMGA) {
  re <- "\\bETA\\[(\\d+)\\]\\b"
  .offset <- as.integer(0)
  lapply(1:length(OMGA), function(k) {
    s <- OMGA[[k]]
    f <- eval(parse(text = (sprintf("y~%s", deparse(s[[2]])))))
    r <- unlist(lapply(attr(terms(f), "variables"), deparse))[-(1:2)]
    nr <- length(r)

    ix <- grep(re, r)
    if (nr - length(ix)) stop("invalid OMGA specs")

    ix <- as.integer(sub(re, "\\1", r))
    if (any(ix - (.offset + 1:nr))) stop("invalid OMGA specs")
    .offset <<- .offset + nr
    eval(s[[3]])
  })
}

genOM <- function(s) {
  getNR <- function(a) round(sqrt(2 * length(a) + 0.25) - 0.1)
  nr <- sum(sapply(s, getNR))
  .mat <- matrix(0, nr, nr)
  .offset <- as.integer(0)
  j <- lapply(1:length(s), function(k) {
    a <- s[[k]]
    p <- getNR(a)
    starts <- row(.mat) > .offset & col(.mat) > .offset
    .mat[col(.mat) >= row(.mat) & col(.mat) <= .offset + p & starts] <<- a
    .offset <<- .offset + p
  })
  a <- .mat[col(.mat) >= row(.mat)]
  .mat <- t(.mat)
  .mat[col(.mat) >= row(.mat)] <- a
  .mat
}

genOMinv.5 <- function(s) {
  getNR <- function(a) round(sqrt(2 * length(a) + 0.25) - 0.1)
  nr <- sum(sapply(s, getNR))
  .mat <- matrix(0, nr, nr)
  .offset <- as.integer(0)
  j <- lapply(1:length(s), function(k) {
    a <- s[[k]]
    p <- getNR(a)
    starts <- row(.mat) > .offset & col(.mat) > .offset
    .mat[col(.mat) >= row(.mat) & col(.mat) <= .offset + p & starts] <<- a
    .offset <<- .offset + p
  })
  .mat
}


#' Print a gnlmm fit
#'
#' Print a generalized non-linear mixed effect model fit
#'
#' @param x a gnlmm fit object
#' @param ... additional arguments
#' @return NULL
#' @export
print.gnlmm.fit <- function(x, ...) {
  x$ETA <- NULL
  x$con <- NULL
  x$calls <- NULL
  x$nsplt <- NULL
  x$osplt <- NULL
  x$obj <- NULL
  x$diag.xform <- NULL
  attr(x, "class") <- NULL
  print.default(x)
}



multi2 <- function(mu, vmat, n) {
  eta <- matrix(rnorm(length(mu) * n), ncol = n, nrow = length(mu))
  Q <- chol(vmat, pivot = TRUE)
  pivot <- attr(Q, "pivot")
  oo <- order(pivot)
  para <- t(Q[, oo]) %*% eta
  sweep(para, 1, mu, "+")
}

#' Prediction after a gnlmm fit
#'
#' Generate predictions after a generalized non-linear mixed effect model fit
#'
#' @param fit a gnlmm fit object
#' @param pred prediction function
#' @param data new data
#' @param mc.cores number of cores (for Linux only)
#' @return observed and predicted
#' @examples
#' \dontrun{
#'
#' ode <- "
#' d/dt(depot) =-KA*depot;
#' d/dt(centr) = KA*depot - KE*centr;
#' "
#' sys1 <- RxODE(ode)
#'
#' pars <- function() {
#'   CL <- exp(THETA[1] + ETA[1]) # ; if (CL>100) CL=100
#'   KA <- exp(THETA[2] + ETA[2]) # ; if (KA>20) KA=20
#'   KE <- exp(THETA[3])
#'   V <- CL / KE
#'   sig2 <- exp(THETA[4])
#' }
#' llik <- function() {
#'   pred <- centr / V
#'   dnorm(DV, pred, sd = sqrt(sig2), log = TRUE)
#' }
#' inits <- list(THTA = c(-3.22, 0.47, -2.45, 0))
#' inits$OMGA <- list(ETA[1] ~ .027, ETA[2] ~ .37)
#' theo <- read.table("theo_md.txt", head = TRUE)
#'
#' fit <- gnlmm(llik, theo, inits, pars, sys1,
#'   control = list(trace = TRUE, nAQD = 5)
#' )
#'
#' pred <- function() {
#'   pred <- centr / V
#' }
#'
#' s <- prediction(fit, pred)
#' plot(s$p, s$dv)
#' abline(0, 1, col = "red")
#' }
#' @export
# new prediction function with new data
prediction <- function(fit, pred, data = NULL, mc.cores = 1) {
  if (!is.null(data)) { # new data
    fit$ETA <- NULL
  } else {
    data <- fit$calls$data
  }
  syspar <- fit$calls$syspar
  system <- fit$calls$system
  th <- fit$par
  con <- fit$con
  diag.xform <- fit$diag.xform
  nsplt <- fit$nsplt
  osplt <- fit$osplt

  square <- function(x) x * x
  diag.xform.inv <- c("sqrt" = "square", "log" = "exp", "identity" = "identity")[diag.xform]

  data.obs <- subset(data, data$EVID == 0)
  names(data) <- tolower(names(data)) # needed in ev

  # options
  ID.all <- unique(data[, "id"])
  ID.ord <- order(ID.all)
  names(ID.ord) <- ID.all
  nsub <- length(ID.all)


  lth <- tapply(th, nsplt, identity)
  THETA <- lth[[1]]
  if (is.null(fit$ETA)) { # with new data
    lD <- tapply(lth[[2]], osplt, identity)
    Dinv.5 <- genOMinv.5(lD)
    diag(Dinv.5) <- eval(call(diag.xform.inv, diag(Dinv.5)))
    D.5 <- solve(Dinv.5)
    D <- crossprod(D.5)
    fit$ETA <- t(multi2(rep(0, dim(D)[1]), D, nsub))
    # print(fit$ETA)
  }

  bpar <- body(syspar)
  blik <- body(pred)

  ep <- environment()
  ..lik.sub <- mclapply(ID.all, function(ix) {
    #-- data
    if (!is.null(system)) {
      ev <- RxODE::eventTable()
      ev$import.EventTable(data[data$id == ix, ])
    }

    sel <- data.obs$ID == ix
    list2env(data.obs[sel, ], envir = ep)

    #-- inner optim: empirical bayesian
    ..g.fn <- function(ETA, dolog = T) {
      if (!is.null(syspar)) {
        eval(bpar)
        pars <- as.list(environment())
        inits <- pars$initCondition
        pars$initCondition <- NULL
        pars <- unlist(pars)
        # print(pars)
      }
      else {
        inits <- NULL
      }

      if (!is.null(system)) {
        m <- system$run(pars, ev, inits, transit_abs = con$transit_abs, atol = con$atol.ode, rtol = con$rtol.ode)
        df <- lapply(2:ncol(m), function(ix) m[, ix])
        names(df) <- dimnames(m)[[2]][-1]
        list2env(df, envir = ep)
      }

      eval(blik)
    }

    .wh <- ID.ord[as.character(ix)]
    ..g.fn(fit$ETA[.wh, ])
  }, mc.cores = mc.cores)

  list(p = unlist(..lik.sub), dv = data.obs$DV)
}


#' Fit a generalized nonlinear mixed-effect model
#'
#' Fit a generalized nonlinear mixed-effect model by adaptive Gaussian quadrature (AQD)
#'
#' @param llik log-likelihood function
#' @param data data to be fitted
#' @param inits initial values
#' @param system an optional (compiled) RxODE object
#' @param syspar function: calculation of PK parameters
#' @param diag.xform transformation to diagonal elements of OMEGA during fitting
#' @param ... additional options
#' @param control additional optimization options
#' @return NULL
#' @details
#'    Fit a generalized nonlinear mixed-effect model by adaptive Gaussian quadrature (AGQ)
#'
#' @author Wenping Wang
#' @examples
#' llik <- function() {
#'   lp <- THETA[1] * x1 + THETA[2] * x2 + (x1 + x2 * THETA[3]) * ETA[1]
#'   p <- pnorm(lp)
#'   dbinom(x, m, p, log = TRUE)
#' }
#' inits <- list(THTA = c(1, 1, 1), OMGA = list(ETA[1] ~ 1))
#'
#' gnlmm(llik, rats, inits, control = list(nAQD = 3))
#' \dontrun{
#' llik <- function() {
#'   if (group == 1) {
#'     lp <- THETA[1] + THETA[2] * logtstd + ETA[1]
#'   } else {
#'     lp <- THETA[3] + THETA[4] * logtstd + ETA[1]
#'   }
#'   lam <- exp(lp)
#'   dpois(y, lam, log = TRUE)
#' }
#' inits <- list(THTA = c(1, 1, 1, 1), OMGA = list(ETA[1] ~ 1))
#'
#' fit <- gnlmm(llik, pump, inits,
#'   control = list(
#'     reltol.outer = 1e-4,
#'     optim.outer = "nmsimplex",
#'     nAQD = 5
#'   )
#' )
#'
#'
#'
#' ode <- "
#' d/dt(depot) =-KA*depot;
#' d/dt(centr) = KA*depot - KE*centr;
#' "
#' sys1 <- RxODE(ode)
#'
#' pars <- function() {
#'   CL <- exp(THETA[1] + ETA[1]) # ; if (CL>100) CL=100
#'   KA <- exp(THETA[2] + ETA[2]) # ; if (KA>20) KA=20
#'   KE <- exp(THETA[3])
#'   V <- CL / KE
#'   sig2 <- exp(THETA[4])
#' }
#' llik <- function() {
#'   pred <- centr / V
#'   dnorm(DV, pred, sd = sqrt(sig2), log = TRUE)
#' }
#' inits <- list(THTA = c(-3.22, 0.47, -2.45, 0))
#' inits$OMGA <- list(ETA[1] ~ .027, ETA[2] ~ .37)
#' # inits$OMGA=list(ETA[1]+ETA[2]~c(.027, .01, .37))
#' theo <- read.table("theo_md.txt", head = TRUE)
#'
#' fit <- gnlmm(llik, theo, inits, pars, sys1,
#'   control = list(trace = TRUE, nAQD = 5)
#' )
#'
#' cv <- calcCov(fit)
#' cbind(fit$par[fit$nsplt == 1], sqrt(diag(cv)))
#' }
#' @export
# this version uses lbfgs instead of optim(); ~ 1/3+ improvement in speed
# require(lbfgs)
gnlmm <- function(llik, data, inits, syspar = NULL,
                  system = NULL, diag.xform = c("sqrt", "log", "identity"),
                  ..., control = list())
                  # TODO:
                  #-- chk ode pars
                  #-- post process
                  #-- gof
                  #-- chk ode_inits
                  #-- t-dist
                  #-- start.zero.inner
                  #-- data
#-- optim.outer
{
  # data
  if (is.null(data$ID)) stop('"ID" not found in data')
  if (is.null(data$EVID)) data$EVID <- 0
  data.obs <- subset(data, data$EVID == 0)
  data.sav <- data
  names(data) <- tolower(names(data)) # needed in ev

  ## model
  RxODE::rxReq("lbfgs")
  if (is.null(system)) {}
  else if (class(system) == "RxODE") {}
  else if (class(system) == "character") {
    obj <- basename(tempfile())
    system <- RxODE(model = system, modName = obj)
  }
  else {
    stop("invalid system input")
  }

  # options
  con <- list(
    trace = 0,
    maxit = 100L,
    atol.ode = 1e-08,
    rtol.ode = 1e-08,
    reltol.inner = 1.0e-6,
    reltol.outer = 1.0e-3,
    optim.inner = "lbfgs",
    optim.outer = "Nelder-Mead",
    start.zero.inner = FALSE,
    mc.cores = 1,
    nAQD = 1,
    transit_abs = FALSE,
    cov = FALSE,
    eps = c(1e-8, 1e-3) # finite difference step
  )
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) {
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }

  square <- function(x) x * x
  diag.xform <- match.arg(diag.xform)
  diag.xform.inv <- c("sqrt" = "square", "log" = "exp", "identity" = "identity")[diag.xform]

  # process inits
  lh <- parseOM(inits$OMGA)
  nlh <- sapply(lh, length)
  osplt <- rep(1:length(lh), nlh)

  lini <- list(inits$THTA, unlist(lh))
  nlini <- sapply(lini, length)
  nsplt <- rep(1:length(lini), nlini)

  om0 <- genOM(lh)
  th0.om <- lapply(1:length(lh), function(k) {
    m <- genOM(lh[k])
    nr <- nrow(m)
    mi <- tryCatch(
      backsolve(chol(m), diag(nr)),
      error = function(e) {
        stop("OMEGA not positive-definite")
      }
    )
    diag(mi) <- eval(call(diag.xform, diag(mi)))
    mi[col(mi) >= row(mi)]
  })
  inits.vec <- c(inits$THTA, unlist(th0.om), inits$SGMA)
  names(inits.vec) <- NULL

  nTHTA <- nlini[1]
  nETA <- nrow(om0)
  ID.all <- unique(data[, "id"])
  ID.ord <- order(ID.all)
  names(ID.ord) <- ID.all
  nSUB <- length(ID.all)

  # gaussian quadrature nodes & wts
  nAQD <- con$nAQD
  nw <- gauss.quad(nAQD)
  mij <- as.matrix(
    do.call("expand.grid", lapply(1:nETA, function(x) 1:nAQD))
  )
  nij <- nrow(mij)


  # obj fn by AQD
  if (!is.null(syspar)) {
    bpar <- body(syspar)
  }
  blik <- body(llik)
  starts <- matrix(0.1, nSUB, nETA)
  omga_save <- NULL
  update_starts <- TRUE

  obj <- function(th, noomga = FALSE) {
    if (noomga) th <- c(th, omga_save)

    lth <- tapply(th, nsplt, identity)
    THETA <- lth[[1]]
    lD <- tapply(lth[[2]], osplt, identity)
    Dinv.5 <- genOMinv.5(lD)
    diag(Dinv.5) <- eval(call(diag.xform.inv, diag(Dinv.5)))
    detDinv.5 <- prod(diag(Dinv.5))

    ep <- environment()
    ..lik.sub <- mclapply(ID.all, function(ix) {
      #-- data
      if (!is.null(system)) {
        ev <- RxODE::eventTable()
        ev$import.EventTable(data[data$id == ix, ])
      }

      sel <- data.obs$ID == ix
      list2env(data.obs[sel, ], envir = ep)

      #-- inner optim: empirical bayesian
      ..g.fn <- function(ETA, dolog = T) {
        if (!is.null(syspar)) {
          eval(bpar)
          pars <- as.list(environment())
          inits <- pars$initCondition
          pars$initCondition <- NULL
          pars <- unlist(pars)

          if (any(is.na(pars))) {
            cat("Error: Na/NaN in pars.\n")
            print(pars)
            cat("Consider to put bounds on pars and/or change optim method.\n")
            stop()
          }
        }
        else {
          inits <- NULL
        }

        if (!is.null(system)) {
          m <- system$run(pars, ev, inits, transit_abs = con$transit_abs, atol = con$atol.ode, rtol = con$rtol.ode)
          if (any(is.na(m))) {
            cat("Error: missing values resulted from ODE solving. ")
            cat("check ODEs and/or parameters.\n")
            print(pars[system$get.modelVars()$params])
            print(m)
            stop()
          }
          df <- lapply(2:ncol(m), function(ix) m[, ix])
          names(df) <- dimnames(m)[[2]][-1]
          list2env(df, envir = ep)
        }

        llik.dat <- sum(eval(blik))
        llik.eta <- -crossprod(Dinv.5 %*% ETA) / 2 - nETA / 2 * log(2 * pi) + log(detDinv.5)
        llik <- llik.dat + llik.eta

        if (is.infinite(llik) || is.na(llik)) {
          msg <- "Infinite/NaN likelihood. Consider to adjust ODE solver tolerance and/or code the likelihood more carefully, e.g., put bounds on some quantities.\n"
          stop(msg)
        }

        if (dolog) {
          -llik
        } else {
          if (llik > 400) llik <- 400
          if (llik < (-700)) llik <- -700
          exp(llik)
        }
      }

      .wh <- ID.ord[as.character(ix)]
      if (con$optim.inner == "lbfgs" || nETA == 1) {
        fg <- function(par) {
          if (identical(par, pvd[[1]])) {
            return(pvd)
          }
          ym <- ..g.fn(par)
          pvd <<- list(par, ym)
        }
        f <- function(x) fg(x)[[2]]
        g <- function(x) {
          eps <- con$eps
          fx <- f(x)
          sapply(1:length(x), function(k) {
            xc <- x
            dx <- abs(xc[k]) * eps[1] + eps[2]
            xc[k] <- xc[k] + dx
            (f(xc) - fx) / dx
          })
        }

        pvd <- NULL
        nfcall <- 0
        ..fit.inner <- lbfgs::lbfgs(f, g, starts[.wh, ], invisible = T, epsilon = 10000 * con$reltol.inner)
        ..fit.inner$hessian <- optimHess(..fit.inner$par, f, g)
      } else {
        ..fit.inner <- optim(
          par = starts[.wh, ], ..g.fn, method = con$optim.inner,
          control = list(reltol = con$reltol.inner), hessian = T
        )
      }

      Ginv.5 <- tryCatch(
        {
          .m <- chol(..fit.inner$hessian)
          backsolve(.m, diag(nETA))
        },
        error = function(e) {
          cat("Warning: Hessian not positive definite\n")
          print(..fit.inner$hessian)
          .m <- ..fit.inner$hessian
          # .m <- chol(.m+diag(nETA)*100)
          # .m[col(.m)!=row(.m)] = .001*.m[col(.m)!=row(.m)]
          .md <- matrix(0, nETA, nETA)
          diag(.md) <- abs(diag(.m)) * 1.1
          .m <- chol(.md)
          backsolve(.m, diag(nETA))
        }
      )
      det.Ginv.5 <- prod(diag(Ginv.5))

      #-- AQD
      ..lik.ij <- mclapply(1:nij, function(ix) {
        ij <- mij[ix, ]
        w <- nw$weights[ij]
        z <- nw$nodes[ij]
        a <- ..fit.inner$par + sqrt(2) * Ginv.5 %*% z
        f1 <- ..g.fn(a, dolog = F)
        f2 <- prod(w * exp(z^2))
        f1 * f2
      }, mc.cores = 1)

      ..lik <- 2^(nETA / 2) * det.Ginv.5 * do.call("sum", ..lik.ij)
      c(-2 * log(..lik), .wh, ..fit.inner$par)
    }, mc.cores = con$mc.cores)

    m <- matrix(unlist(..lik.sub), ncol = 2 + nETA, byrow = T)

    if (update_starts) starts[m[, 2], ] <<- m[, 3:(2 + nETA)]
    r <- sum(m[, 1])
    attr(r, "subj") <- m[, 1]
    r
  }

  args <- list(inits.vec, obj, control = list(trace = con$trace, reltol = con$reltol.outer))

  fit <- if (con$optim.outer == "nmsimplex") {
    do.call("nmsimplex", args)
  } else {
    args <- c(args, method = con$optim.outer)
    do.call("optim", args)
  }

  fit <- c(fit, obj = obj, list(ETA = starts, con = con, diag.xform = diag.xform, nsplt = nsplt, osplt = osplt, calls = list(data = data.sav, system = system, syspar = syspar)))
  attr(fit, "class") <- "gnlmm.fit"
  fit
}


linesearch_secant <- function(f, d, x, maxIter = 5, trace = F) {
  # Line search using secant method
  # Note: I'm not checking for alpha > 0.

  epsilon <- 10^(-3) # line search tolerance
  # maxIter = 5; #maximum number of iterations
  alpha_curr <- 0
  alpha <- 10^(-5)
  s <- f(x)
  y <- s$y
  grad <- s$grad
  dphi_zero <- crossprod(grad, d)
  dphi_curr <- dphi_zero

  i <- 0
  while (all(abs(dphi_curr) > epsilon * abs(dphi_zero))) {
    alpha_old <- alpha_curr
    alpha_curr <- alpha
    dphi_old <- dphi_curr
    s <- f(x + alpha_curr * d)
    y <- s$y
    grad <- s$grad
    dphi_curr <- crossprod(grad, d)
    alpha <- (dphi_curr * alpha_old - dphi_old * alpha_curr) / (dphi_curr - dphi_old)
    dim(alpha) <- NULL
    if (!is.finite(alpha)) break
    i <- i + 1
    if ((i >= maxIter) && (abs(dphi_curr) > epsilon * abs(dphi_zero))) {
      s <- "Line search terminating with number of iterations: %d"
      warning(sprintf(s, i))
      break
    }
  } # while
  if (trace) print(i)

  alpha
}

#' Calculate gnlmm variance-covariance matrix of fixed effects
#'
#' Calculate variance-covariance matrix of fixed effects after a gnlmm() fit
#'
#' @param fit a gnlmm fit object
#' @param method method for calculating variance-covariance matrix
#' @param trace logical whether to trace the iterations
#' @return variance-covariance matrix of model parameters
#' @export
calcCov <- function(fit, method = 1, trace = FALSE) {
  lth <- tapply(fit$par, fit$nsplt, identity)
  e1 <- environment(fit$obj)
  e1$omga_save <- lth[[2]]
  e1$update_starts <- FALSE
  # e1$starts = .1+0*e1$starts

  pvd <- NULL
  nfcall <- 0 ## CHECKME
  fg <- function(par) {
    if (identical(par, pvd[[1]])) {
      return(pvd)
    }
    ym <- fit$obj(par, noomga = T)
    pvd <<- list(par, ym)
  }
  f <- function(x) fg(x)[[2]]
  g <- function(x) {
    # eps = con$eps
    eps <- c(1.6e-6, 1.6e-4)
    fx <- f(x)
    fxi <- attr(fx, "subj")
    lx <- length(x)
    sm <- sapply(1:lx, function(k) {
      xc <- x
      dx <- abs(xc[k]) * eps[1] + eps[2]
      xc[k] <- xc[k] + dx
      fxc <- f(xc)
      fxci <- attr(fxc, "subj")
      da <- (fxc - fx) / dx
      di <- (fxci - fxi) / dx
      c(da, di)
    })
    r <- sm[1, ]
    attr(r, "subj") <- t(sm[-1, ])
    r
  }

  FuncGrad <- function(x) {
    pvd <<- NULL
    nfcall <- 0
    list(y = f(x), grad = g(x))
  }

  x <- lth[[1]]
  lx <- length(x)
  lx <- length(x)
  pvd <- NULL
  nfcall <- 0 ## CHECKME
  x <- g(x)
  x <- attr(x, "subj")
  sm <- matrix(apply(apply(x, 2, function(y) y %*% t(y)), 1, sum), lx, lx)
  sinv <- solve(sm)
  H <- 2 * diag(diag(sinv))
  # H = pmax(2*diag(diag(sinv)), diag(lx));
  H <- 2 * sinv
  x <- lth[[1]]

  for (k in 1:10) {
    if (trace) print(k)
    l <- FuncGrad(x)
    value <- l$y
    grad <- l$grad
    p <- -H %*% grad
    alpha <- linesearch_secant(FuncGrad, p, x, trace = trace)
    if (trace) print(alpha)
    x <- x + alpha * p
    s <- alpha * p
    dist <- sqrt(sum(s^2))
    l <- FuncGrad(x)
    newvalue <- l$y
    newgrad <- l$grad
    y <- newgrad - grad
    rho <- 1 / (t(y) %*% s)
    dim(rho) <- NULL
    H <- (diag(lx) - rho * s %*% t(y)) %*% H %*% (diag(lx) - rho * y %*% t(s)) + rho * s %*% t(s)
    # print(c(alpha, rho, x)); #
    if (trace) print(diag(H))
    if (dist < .01) break
  }

  cv <- if (method == 1) {
    2 * H
  } else if (method == 2) {
    4 * sinv
  } else {
    (H %*% sinv %*% H)
  }
  attr(cv, "RinvS") <- list(Rinv = H, S = sm)
  cv
}


#' Calculate gnlmm variance-covariance matrix of random effects
#'
#' Calculate variance-covariance matrix of random effects after a gnlmm() fit
#'
#' @param fit a gnlmm fit object
#' @return variance-covariance matrix of random effects
#' @export
getOMEGA <- function(fit) {
  th <- fit$par
  diag.xform <- fit$diag.xform
  nsplt <- fit$nsplt
  osplt <- fit$osplt

  square <- function(x) x * x
  diag.xform.inv <- c("sqrt" = "square", "log" = "exp", "identity" = "identity")[diag.xform]

  lth <- tapply(th, nsplt, identity)
  THETA <- lth[[1]]
  lD <- tapply(lth[[2]], osplt, identity)
  Dinv.5 <- genOMinv.5(lD)
  diag(Dinv.5) <- eval(call(diag.xform.inv, diag(Dinv.5)))
  D.5 <- solve(Dinv.5)
  D <- crossprod(D.5)

  shrinkageETA <- 100 * (1 - apply(fit$ETA, 2, var) / diag(D))
  shrinkageETA <- ifelse(shrinkageETA > 0, shrinkageETA, 0)

  list(OMEGA = D, shrinkageETA = shrinkageETA)
}
