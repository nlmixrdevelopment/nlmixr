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

## require(lbfgs)
## require(lbfgsb3)
## require(madness)
## require(Rcpp)
# require(nlmixr)
## require(parallel)
## require(minqa)
## require(Deriv)

llik_binomial <- function(y, n, params) {
  r <- llik_binomial_c(y, n, params)
  r$J <- diag(r$J)
  return(r)
}


#-- for new gnlmm
getModelVars <- function(blik, bpar, m1) {
  argsList <- list(
    dpois = list(ix = 2, dvdx = ".arg1"),
    dbinom = list(ix = 2:3, dvdx = ".arg2"),
    dnorm = list(ix = 2:3, dvdx = c(".arg1", ".arg2")),
    dbeta = list(ix = 2:3, dvdx = c(".arg1", ".arg2")),
    dneg_binomial = list(ix = 2:3, dvdx = c(".arg1", ".arg2")),
    dbetabinomial = list(ix = 2:4, dvdx = c(".arg2", ".arg3")),
    dt = list(ix = 2:4, dvdx = c(".arg1", ".arg2", ".arg3"))
  )

  blik.txt <- deparse(blik)
  len <- length(blik.txt)

  s <- gsub("\\s+", "", blik.txt[len - 1], perl = T) # FIXME
  lp <- regexpr("\\(", s)
  dist <- substr(s, 1, lp - 1)
  s <- strsplit(substr(s, lp + 1, 200), ",")[[1]]

  args <- argsList[[dist]]
  args.ix <- args$ix # position of dist pars
  args.dvdx <- args$dvdx # args need dvdx
  narg <- length(args.ix)
  blik.new.text <- paste(c(
    blik.txt[2:(len - 2)],
    paste0(".arg", 1:narg, "=", s[args.ix])
  ), collapse = "\n")
  blik.new <- parse(text = blik.new.text)

  dist.df <- NULL
  if (dist == "dbinom") dist.df <- s[2] # binomial size
  if (dist == "dt") dist.df <- s[2] # t df

  #----------------------------
  f <- deparse(blik)
  len <- length(f)
  f <- parse(text = f[c(-1, -len)])
  out <- utils::getParseData(f)
  s <- out$text[out$token %in% c("SYMBOL", "LEFT_ASSIGN", "EQ_ASSIGN")]
  ix <- s %in% c("=", "<-")
  lhs <- (1:length(ix))[ix] - 1
  ix[lhs] <- TRUE
  rhs <- setdiff(s, s[ix])

  ixLlik <- len - 1
  ss <- sub("^\\s*\\w+\\(", "", deparse(blik)[len - 1], perl = T) # FIXME
  prob <- strsplit(ss, ",")[[1]][3] # FIXME

  states <- m1$get.modelVars()$state
  state.llik <- intersect(states, s[!ix]) # state used in llik

  f <- deparse(bpar)
  len <- length(f)
  f <- parse(text = f[c(-1, -len)])
  out <- utils::getParseData(f)
  s <- out$text[out$token %in% c("SYMBOL", "LEFT_ASSIGN", "EQ_ASSIGN")]
  ix <- s %in% c("=", "<-")
  lhs <- (1:length(ix))[ix] - 1
  pars.llik <- intersect(s[lhs], rhs) # pars used in llik
  vars.par <- s[lhs] # vars def'ed in pars

  list(
    state.llik = state.llik, pars.llik = pars.llik,
    vars.par = vars.par, # prob=prob, ixLlik=ixLlik,
    dist = dist, dist.df = dist.df,
    blik.new = blik.new, blik.new.text = blik.new.text,
    args.dvdx = args.dvdx
  )
}

##' @rdname gnlmm
##' @export
gnlmm2 <- function(llik, data, inits, syspar = NULL,
                   system = NULL, diag.xform = c("sqrt", "log", "identity"),
                   ..., control = list()) {
  ## data
  if (is.null(data$ID)) stop('"ID" not found in data')
  if (is.null(data$EVID)) data$EVID <- 0
  data.obs <- data[data$EVID == 0, ]
  data.sav <- data
  names(data) <- tolower(names(data)) # needed in ev

  # model
  if (is.null(system)) {}
  else if (class(system) == "RxODE") {
    system <- RxODE(system, calcSens = TRUE)
  }
  else if (class(system) == "character") {
    system <- RxODE(model = system, calcSens = TRUE)
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
    reltol.inner = 1.0e-4,
    reltol.outer = 1.0e-3,
    optim.inner = "lbfgs",
    optim.outer = "newuoa",
    start.zero.inner = FALSE,
    mc.cores = 1,
    nAQD = 1,
    transit_abs = FALSE,
    cov = FALSE,
    eps = c(1e-8, 1e-3), # finite difference step
    NOTRUN = F,
    DEBUG.INNER = F,
    rhobeg = .2,
    rhoend = 1e-3,
    iprint = 2,
    npt = NULL,
    do.optimHess = TRUE
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
  modVars <- getModelVars(blik, bpar, system)


  # === start of dvdx code
  proc.deriv <- function() {
    getDeriv <- function(pars) {
      npar <- length(pars)
      s <- deparse(bpar)
      for (i in 1:nETA) {
        s <- gsub(sprintf("\\bETA\\[%d\\]", i), sprintf("ETA%d", i), s, perl = T)
      }
      s <- gsub("initCondition", "#initCondition", s)
      len <- length(s)
      s <- s[-len]

      ix <- matrix(unlist(expand.grid(1:npar, 1:nETA)), ncol = 2)
      m <- sapply(1:nrow(ix), function(k) {
        i <- ix[k, 1]
        j <- ix[k, 2]
        s1 <- c("Deriv::Deriv(~", s, pars[i], sprintf("}, \"ETA%d\")", j))
        a <- paste(s1, collapse = "\n")
        e <- eval(parse(text = a))
        s <- if (class(e) == "call") {
          s <- deparse(e)
          gsub(sprintf("\\bETA%d\\b", j), sprintf("ETA[%d]", j), s, perl = T)
        } else {
          "0"
        }
        # parse(text=s)
        s
      })
      m
    }

    env <- environment()
    dati <- data.sav[data.sav$ID == ID.all[1], ]
    list2env(dati, env)
    THETA <- tapply(inits.vec, nsplt, identity)[[1]]
    ETA <- madness::madness(array(0, c(nETA, 1))) ## Is madness still needed...?
    eval(bpar)

    px <- as.list(env)
    madIx <- sapply(px, function(x) {
      if (class(x) == "madness") TRUE else FALSE
    })
    madVars <- names(px)[madIx]

    pars <- system$get.modelVars()$params
    px <- as.list(env)[pars]
    madIx <- sapply(px, function(x) {
      if (class(x) == "madness") TRUE else FALSE
    })

    pars <- pars[madIx]
    matode <- getDeriv(pars)
    pars <- setdiff(madVars, c(pars, "ETA"))
    matllk <- getDeriv(pars)

    list(matode = matode, madIx = madIx, matllk = matllk, madllk = pars) # deriv expr for ode pars
    # idx of ode pars that need deriv
    # llik pars that need deriv
    # deriv expr llik pars
  }
  s <- proc.deriv()

  # FIXME: chk par order of below 2 para against those in d(llik)/d(ETA) by chain rule
  # d(pars)/d(ETA) for ode
  madIx <- s$madIx
  npar <- sum(madIx)
  m <- s$matode
  idx.dpde.ode <- m != "0"
  expr.dpde.ode <- lapply(m[idx.dpde.ode], function(k) parse(text = k))
  dpde.ode <- matrix(0, npar, nETA)

  # d(pars)/d(ETA) for llik
  madVars.llk <- s$madllk
  npar <- length(madVars.llk)
  m <- s$matllk
  idx.dpde.llk <- m != "0"
  expr.dpde.llk <- lapply(m[idx.dpde.llk], function(k) parse(text = k))
  dpde.llk <- matrix(0, npar, nETA)
  dimnames(dpde.llk) <- list(madVars.llk, NULL)

  # d(args)/d(pars) for llik
  pars <- c(modVars$state.llik, madVars.llk)
  npar <- length(pars)
  lexpr <- unlist(lapply(modVars$args.dvdx, function(arg) {
    s <- sprintf("{\n%s\n%s", modVars$blik.new.text, arg) # FIXME

    expr.dldp <- lapply(1:npar, function(k) {
      s1 <- c("Deriv::Deriv(~", s, sprintf("}, \"%s\")", pars[k]))
      a <- paste(s1, collapse = "\n")
      e <- eval(parse(text = a))
      e
    })
    names(expr.dldp) <- pars
    expr.dldp
  }))

  llik.narg <- length(modVars$args.dvdx) # llik.narg = # of args in density that need dvdx
  # may have efficiency gain to rm args that do not need dvdx
  llik.npar <- npar
  m <- matrix(lexpr, llik.npar, llik.narg)
  dadp.expr <- c(t(m)) # chg the order of args & pars in llik; see dldp in ..fg()
  # === ends of dvdx code


  starts <- matrix(0., nSUB, nETA)
  omga_save <- NULL
  update_starts <- TRUE


  # algo starts
  obj.vec <- function(th, noomga = FALSE) {
    th <- th * inits.vec
    # if (con$DEBUG.INNER) print(th)

    if (noomga) th <- c(th, omga_save)

    lth <- tapply(th, nsplt, identity)
    THETA <- lth[[1]]
    lD <- tapply(lth[[2]], osplt, identity)
    Dinv.5 <- genOMinv.5(lD)
    diag(Dinv.5) <- eval(call(diag.xform.inv, diag(Dinv.5)))
    detDinv.5 <- prod(diag(Dinv.5))

    ep <- environment()

    llik2.subj <- function(ix) {
      # ix=4; ETA=rep(0, nETA)
      dati <- data.sav[data.sav$ID == ix, ]
      evi <- dati[, c("TIME", "EVID", "AMT")]
      names(evi) <- tolower(names(evi))
      ev <- RxODE::eventTable()
      ev$import.EventTable(evi)
      dati <- dati[dati$EVID == 0, ]


      # ETA => pars & stateVar => .args => llik
      # d(llik)/d(ETA) = d(llik)/d(args) * d(args)/d(ETA)
      # d(args)/d(ETA) = d(args)/d(pars) * d(pars)/d(ETA)
      # stateVar is parallel to pars when forming args, however, stateVar = f(pars(ETA))
      # hence, we need d(State)/d(ETA). d(State)/d(ETA) = d(State)/d(pars) * d(pars)/d(ETA)
      ..g.fn <- function(ETA) {
        env <- environment()
        list2env(dati, env)
        eval(bpar)

        pars <- system$get.modelVars()$params
        po <- unlist(as.list(env)[pars])
        x <- system$run(po, ev, initCondition)
        if (any(is.na(x))) {
          print(ID[1])
          print(po)
          print(head(x, 10))
          stop("NA in ODE solution")
        }

        # d(State)/d(ETA)
        whState <- modVars$state.llik
        senState <- paste0("rx__sens_", whState, "_BY_", pars[madIx], "__")
        fxJ <- list(fx = x[, whState], J = x[, senState]) # FIXME, t()

        dvdx <- sapply(expr.dpde.ode, eval, envir = env) # FIXME
        dpde.ode[idx.dpde.ode] <- dvdx
        # d(State)/d(ETA) = d(State)/d(pars) * d(pars)/d(ETA)
        dvdxState <- fxJ$J %*% dpde.ode # FIXME

        # make state var & other symbols available for calc
        valState <- fxJ$fx
        assign(whState, valState, envir = env)
        eval(modVars$blik.new) # FIXME: here or at llik_binomial?

        # d(pars)/d(ETA)
        if (length(expr.dpde.llk) > 0) {
          dvdx <- sapply(expr.dpde.llk, eval, envir = env) # FIXME
          dpde.llk[idx.dpde.llk] <- dvdx
        } else {
          dpde.llk <- NULL
        }

        # d(args)/d(pars)
        ni <- dim(x)[1]
        dadp <- sapply(1:length(dadp.expr), function(k) { # why lapply(m, eval) doesn't work?
          s <- eval(dadp.expr[[k]])
          if (length(s) == 1) rep(s, ni) else s
        })
        dim(dadp) <- c(ni, llik.narg, llik.npar) # FIXME

        # d(args)/d(ETA) = d(args)/d(pars) * d(pars)/d(ETA)
        dade <- sapply(1:ni, function(k) { # FIXME: vectorize?
          dpde <- rbind(dvdxState[k, ], dpde.llk)
          dadp[k, , ] %*% dpde # FIXME: need t()?
        })
        # dade = t(s)										#FIXME: t() can be rm'ed?


        dim(dade) <- c(llik.narg, nETA, ni)


        # d(llik)/d(ETA) = d(llik)/d(args) * d(args)/d(ETA)
        if (modVars$dist == "dt") {
          if (length(.arg1) == 1 && ni > 1) .arg1 <- rep(.arg1, ni)
          if (length(.arg2) == 1 && ni > 1) .arg2 <- rep(.arg2, ni)
          if (length(.arg3) == 1 && ni > 1) .arg3 <- rep(.arg3, ni)
          s <- sapply(1:ni, function(k) { # FIXME: vectorize?
            unlist(llik_student_t(DV[k], c(.arg1[k], .arg2[k], .arg3[k])))
          })
          s <- list(fx = s[1, ], J = t(s[-1, ]))
        } else if (modVars$dist == "dbetabinomial") {
          if (length(.arg1) == 1 && ni > 1) .arg1 <- rep(.arg1, ni)
          if (length(.arg2) == 1 && ni > 1) .arg2 <- rep(.arg2, ni)
          if (length(.arg3) == 1 && ni > 1) .arg3 <- rep(.arg3, ni)
          s <- sapply(1:ni, function(k) { # FIXME: vectorize?
            unlist(llik_betabinomial(DV[k], .arg1[k], c(.arg2[k], .arg3[k])))
          })
          s <- list(fx = s[1, ], J = t(s[-1, ]))
        } else if (modVars$dist == "dneg_binomial") {
          if (length(.arg1) == 1 && ni > 1) .arg1 <- rep(.arg1, ni)
          if (length(.arg2) == 1 && ni > 1) .arg2 <- rep(.arg2, ni)
          s <- sapply(1:ni, function(k) { # FIXME: vectorize?
            unlist(llik_neg_binomial(DV[k], c(.arg1[k], .arg2[k])))
          })
          s <- list(fx = s[1, ], J = t(s[-1, ]))
        } else if (modVars$dist == "dbeta") {
          if (length(.arg1) == 1 && ni > 1) .arg1 <- rep(.arg1, ni)
          if (length(.arg2) == 1 && ni > 1) .arg2 <- rep(.arg2, ni)
          s <- sapply(1:ni, function(k) { # FIXME: vectorize?
            unlist(llik_beta(DV[k], c(.arg1[k], .arg2[k])))
          })
          s <- list(fx = s[1, ], J = t(s[-1, ]))
        } else if (modVars$dist == "dnorm") {
          if (length(.arg1) == 1 && ni > 1) .arg1 <- rep(.arg1, ni)
          if (length(.arg2) == 1 && ni > 1) .arg2 <- rep(.arg2, ni)
          s <- sapply(1:ni, function(k) { # FIXME: vectorize?
            if (.arg2[k] < 0.0001) {
              # print("HAY!");print(.arg1[k]);print(.arg2[k])
            }
            unlist(llik_normal(DV[k], c(.arg1[k], .arg2[k])))
          })
          s <- list(fx = s[1, ], J = t(s[-1, ]))
        } else if (modVars$dist == "dbinom") {
          if (length(.arg1) == 1 && ni > 1) .arg1 <- rep(.arg1, ni)
          s <- llik_binomial(DV, .arg1, c(.arg2))
        } else if (modVars$dist == "dpois") {
          s <- llik_poisson(DV, c(.arg1))
        } else {
          stop("dist not supported")
        }

        dim(s$J) <- c(ni, llik.narg)
        dlde <- sapply(1:ni, function(k) { # FIXME: vectorize?
          s$J[k, ] %*% dade[, , k] # FIXME: need t()?
        })
        s <- rbind(s$fx, dlde) # FIXME t()?
        s <- apply(s, 1, sum) # FIXME sum index; chg'ed w/ vec stan call

        # llik.dat = madness::madness(val=matrix(s[1], 1, 1), dvdx=matrix(s[-1], 1, nETA))
        # llik.eta = -crossprod(Dinv.5 %*% ETA)/2 -nETA/2*log(2*pi)+log(detDinv.5)
        # llik.eta = madness::madness(val=val(llik.eta), dvdx=dvdx(llik.eta))
        llik.eta.val <- -crossprod(Dinv.5 %*% ETA) / 2 - nETA / 2 * log(2 * pi) + log(detDinv.5)
        llik.eta.dvd <- -t(ETA) %*% crossprod(Dinv.5)

        r <- s[1] + c(llik.eta.val)
        attr(r, "dvdx") <- s[-1] + c(llik.eta.dvd)
        r
      }
      fg <- function(par) {
        if (identical(par, pvd[[1]])) {
          return(pvd)
        }
        ym <- ..g.fn(par)
        pvd <<- list(par, ym)
      }
      f <- function(ETA) -as.vector(fg(ETA)[[2]])
      g <- function(ETA) -as.vector(attr(fg(ETA)[[2]], "dvdx"))

      pvd <- NULL

      .wh <- ID.ord[as.character(ix)]
      ETA.val <- starts[.wh, ]
      ..fit.inner <- nlminb(ETA.val, f, g, control = list(trace = FALSE, rel.tol = 1e-4))
      # ..fit.inner = lbfgs(f, g, ETA.val, invisible=T, ftol=1e-4) # epsilon=1e-3)
      if (con$do.optimHess) {
        ..fit.inner$hessian <- optimHess(..fit.inner$par, f, g)
      }
      if (con$DEBUG.INNER) {
        # print(..fit.inner$message)
      }

      # =========================================================
      if (con$do.optimHess) {
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
      } else {
        Ginv.5 <- chol(..fit.inner$Hessian.inv)
      }
      det.Ginv.5 <- prod(diag(Ginv.5))

      ## -- AQD
      ..lik.ij <- lapply(1:nij, function(ix) {
        ij <- mij[ix, ]
        w <- nw$weights[ij]
        z <- nw$nodes[ij]
        a <- ..fit.inner$par + sqrt(2) * Ginv.5 %*% z
        f1 <- exp(as.vector(..g.fn(a))) # FIXME
        f2 <- prod(w * exp(z^2))
        f1 * f2
      })

      ..lik <- 2^(nETA / 2) * det.Ginv.5 * do.call("sum", ..lik.ij)
      c(-2 * log(..lik), .wh, ..fit.inner$par)
    }
    s <- mclapply(ID.all, llik2.subj, mc.cores = con$mc.cores) # FIXME
    m <- matrix(unlist(s), ncol = 2 + nETA, byrow = T)

    if (update_starts) starts[m[, 2], ] <<- m[, 3:(2 + nETA)]
    m[, 1]
  }

  nobjcall <- 0
  obj <- function(th, noomga = FALSE) {
    nobjcall <<- nobjcall + 1
    s <- obj.vec(th, noomga)
    r <- sum(s)
    if (con$DEBUG.INNER) {
      print(rbind(
        c(nobjcall, r, th),
        c(nobjcall, r, th * inits.vec)
      ))
    }
    attr(r, "subj") <- s
    r
  }

  np <- length(inits.vec)
  start <- rep(1, np)
  args <- list(start, obj, control = list(trace = con$trace, reltol = con$reltol.outer))

  if (!con$NOTRUN) {
    fit <- if (con$optim.outer == "nmsimplex") {
      do.call("nmsimplex", args)
    } else if (con$optim.outer == "Nelder-Mead") {
      args <- c(args, method = con$optim.outer)
      do.call("optim", args)
    }
    else if (con$optim.outer == "nlminb") {
      args <- list(start, obj, control = list(trace = con$trace, rel.tol = con$reltol.outer))
      do.call("nlminb", args)
    }
    else {
      if (!is.null(con$npt)) {
        npt <- con$npt
      } else {
        npt <- 2 * np + 1
      }
      minqa::newuoa(start, obj, control = list(rhobeg = con$rhobeg, rhoend = con$rhoend, npt = npt, iprint = con$iprint))
    }
  } else {
    fit <- NULL
  }

  fit <- c(fit, obj = obj, list(ETA = starts, con = con, diag.xform = diag.xform, nsplt = nsplt, osplt = osplt, calls = list(data = data.sav, system = system, syspar = syspar)))
  fit$par.unscaled <- fit$par * inits.vec
  attr(fit, "class") <- "gnlmm.fit"
  fit
}


#------------------
mat.indices <- function(nETA) {
  idx <- do.call(
    "rbind",
    lapply(1:nETA, function(k) cbind(k:nETA, k))
  )
  H <- matrix(1:(nETA^2), nETA, nETA)
  Hlo.idx <- row(H) >= col(H)
  lo.idx <- H[row(H) > col(H)]
  hi.idx <- t(H)[row(H) > col(H)]

  list(
    idx = idx, # (r, c) of lo-half
    Hlo.idx = Hlo.idx, # index of lo-half
    lo.idx = lo.idx, # index of strict lo-half
    hi.idx = hi.idx
  ) # index of strict hi-half
}
