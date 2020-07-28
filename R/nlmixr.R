.onLoad <- function(libname, pkgname) {
}

orig.onAttach <- function(libname, pkgname) {
  ## nocov start
  ## Setup RxODE.prefer.tbl
  nlmixrSetupMemoize()
  options(keep.source = TRUE)
  ## nocov end
}

##' Convert data to RxODE format (depreciated)
##'
##' @param data Data to "convert"
##'
##' @return Exact same data as was input
##'
##' @keywords internal
##' @export
nmDataConvert <- function(data) {
  warning("nmDataConvert is depreciated and no longer needed.")
  data
}

nlmixrSetupMemoize <- function() {
  reSlow <- rex::rex(".slow", end)
  f <- sys.function(-1)
  ns <- environment(f)
  .slow <- ls(pattern = reSlow, envir = ns)
  for (slow in .slow) {
    fast <- sub(reSlow, "", slow)
    if (!memoise::is.memoised(get(fast, envir = ns)) && is.null(get(slow, envir = ns))) {
      utils::assignInMyNamespace(slow, get(fast, envir = ns))
      utils::assignInMyNamespace(fast, memoise::memoise(get(slow, envir = ns)))
    }
  }
}

##' Clear memoise cache for nlmixr
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nlmixrForget <- function() {
  reSlow <- rex::rex(".slow", end)
  f <- sys.function(-1)
  ns <- environment(f)
  .slow <- ls(pattern = reSlow, envir = ns)
  for (slow in .slow) {
    fast <- sub(reSlow, "", slow)
    memoise::forget(get(fast, envir = ns))
  }
}


##' @importFrom stats predict logLik na.fail pchisq
##' @importFrom n1qn1 n1qn1
##' @importFrom brew brew
##' @importFrom nlme nlme fixed.effects random.effects
##' @importFrom nlme groupedData
##' @importFrom nlme getData
##' @importFrom nlme pdDiag
##' @importFrom RxODE RxODE
##' @importFrom graphics abline lines matplot plot points title
##' @importFrom stats as.formula nlminb optimHess rnorm terms predict anova optim sd var AIC BIC asOneSidedFormula coef end fitted resid setNames start simulate nobs qnorm quantile time
##' @importFrom utils assignInMyNamespace getFromNamespace head stack sessionInfo tail str
##' @importFrom parallel mclapply
##' @importFrom methods is
##' @importFrom Rcpp evalCpp
##' @importFrom ggplot2 ggplot aes geom_point facet_wrap geom_line geom_abline xlab geom_smooth aes_string
##' @importFrom RcppArmadillo armadillo_version
##' @useDynLib nlmixr, .registration=TRUE


rex::register_shortcuts("nlmixr")
## GGplot use and other issues...
utils::globalVariables(c("DV", "ID", "IPRED", "IRES", "PRED", "TIME", "grp", "initCondition", "values", "nlmixr_pred", "iter", "val", "EVID"))

nlmixr.logo <- "         _             _             \n        | | %9s (_) %s\n  _ __  | | _ __ ___   _ __  __ _ __\n | '_ \\ | || '_ ` _ \\ | |\\ \\/ /| '__|\n | | | || || | | | | || | >  < | |\n |_| |_||_||_| |_| |_||_|/_/\\_\\|_|\n"

##' Messages the nlmixr logo...
##'
##' @param str String to print
##' @param version Version information (by default use package version)
##' @author Matthew L. Fidler
nlmixrLogo <- function(str = "", version = sessionInfo()$otherPkgs$nlmixr$Version) {
  message(sprintf(nlmixr.logo, str, version))
}
##' Display nlmixr's version
##'
##' @author Matthew L. Fidler
##' @export
nlmixrVersion <- function() {
  nlmixrLogo()
}

armaVersion <- function() {
  nlmixrLogo(str = "RcppArmadiilo", RcppArmadillo::armadillo_version())
}

.nlmixrTime <- NULL

##' nlmixr fits population PK and PKPD non-linear mixed effects models.
##'
##' nlmixr is an R package for fitting population pharmacokinetic (PK)
##' and pharmacokinetic-pharmacodynamic (PKPD) models.
##'
##' The nlmixr generalized function allows common access to the nlmixr
##' estimation routines.
##'
##' @template uif
##'
##' @param object Fitted object or function specifying the model.
##' @inheritParams nlmixr_fit
##' @param ... Other parameters
##' @param save Boolean to save a nlmixr object in a rds file in the
##'     working directory.  If \code{NULL}, uses option "nlmixr.save"
##' @return Either a nlmixr model or a nlmixr fit object
##' @author Matthew L. Fidler, Rik Schoemaker
##' @export
nlmixr <- function(object, data, est = NULL, control = list(),
                   table = tableControl(), ..., save = NULL,
                   envir = parent.frame()) {
  assignInMyNamespace(".nlmixrTime", proc.time())
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  force(est)
  ## verbose?
  ## https://tidymodels.github.io/model-implementation-principles/general-conventions.html
  UseMethod("nlmixr")
}

##' @rdname nlmixr
##' @export
nlmixr.function <- function(object, data, est = NULL, control = list(), table = tableControl(), ...,
                            save = NULL, envir = parent.frame()) {
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .modName <- deparse(substitute(object))
  .uif <- nlmixrUI(object)
  class(.uif) <- "list"
  .uif$nmodel$model.name <- .modName
  if (missing(data) && missing(est)) {
    class(.uif) <- "nlmixrUI"
    return(.uif)
  } else {
    .uif$nmodel$data.name <- deparse(substitute(data))
    class(.uif) <- "nlmixrUI"
    .args$data <- data
    .args$est <- est
    .args <- c(list(uif = .uif), .args[-1])
    if (is.null(est)) {
      stop("Need to supply an estimation routine with est=.")
    }
    return(do.call(nlmixr_fit, .args, envir = envir))
  }
}

##' @rdname nlmixr
##' @export
nlmixr.nlmixrFitCore <- function(object, data, est = NULL, control = list(), table = tableControl(), ...,
                                 save = NULL, envir = parent.frame()) {
  .uif <- .getUif(object)
  if (missing(data)) {
    data <- getData(object)
  }
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .args$data <- data
  .args$est <- est
  .args <- c(list(uif = .uif), .args[-1])
  return(do.call(nlmixr_fit, .args, envir = envir))
}

##' @rdname nlmixr
##' @export
nlmixr.nlmixrUI <- function(object, data, est = NULL, control = list(), ...,
                            save = NULL, envir = parent.frame()) {
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .uif <- object
  if (missing(data) && missing(est)) {
    return(.uif)
  } else {
    .args <- c(list(uif = .uif), .args[-1])
    if (missing(data) && !is.null(.getPipedData())) {
      data <- .getPipedData()
      .args$data <- data
      .args$est <- est
    } else {
      .uif$nmodel$data.name <- deparse(substitute(data))
      .args$data <- data
      .args$est <- est
    }
    return(do.call(nlmixr_fit, .args, envir = envir))
  }
}

##' Convert/Format the data appropriately for nlmixr
##'
##' @param data is the name of the data to convert.  Can be a csv file
##'     as well.
##' @param model This is the RxODE model to use to translate against
##'     when parsing the data.
##' @return Appropriately formatted data
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nlmixrData <- function(data, model = NULL) {
  UseMethod("nlmixrData")
}
##' @export
##' @rdname nlmixrData
nlmixrData.character <- function(data, model = NULL) {
  if (!file.exists(data)) {
    stop(sprintf("%s does not exist.", data))
  }
  if (regexpr(rex::rex(".csv", end), data) != -1) {
    return(nlmixrData.default(utils::read.csv(data, na.strings = c(".", "NA", "na", ""))))
  } else {
    stop(sprintf("Do not know how to read in %s", data))
  }
}
##' @export
##' @rdname nlmixrData
nlmixrData.default <- function(data, model = NULL) {
  if (!is.null(model)) {
    dat <- RxODE::etTrans(data, model, addCmt = TRUE, dropUnits = TRUE, allTimeVar = TRUE)
  } else {
    dat <- .as.data.frame(data)
  }
  return(dat)
}
nlmixr_fit0 <- function(uif, data, est = NULL, control = list(), ...,
                        sum.prod = FALSE, table = tableControl(),
                        envir = parent.frame()) {
  if (is.null(est)) {
    stop("Estimation type must be specified by est=''")
  }
  .clearPipedData()
  .tmp <- deparse(body(uif$theta.pars))[-1]
  .tmp <- .tmp[-length(.tmp)]
  .origData <- data
  data <- RxODE::etTrans(data, paste(paste(.tmp, collapse = "\n"), "\n", uif$rxode), TRUE, TRUE, TRUE)

  .nTv <- attr(class(data), ".RxODE.lst")$nTv
  .lab <- attr(class(data), ".RxODE.lst")$idLvl
  .nid <- attr(class(data), ".RxODE.lst")$nid
  .modelId <-
    digest::digest(
      list(
        sessionInfo()$otherPkgs$nlmixr$Version,
        uif, data, est, control, sum.prod, table, ...
      )
    )
  .meta <- uif$meta
  .missingEst <- is.null(est)
  if (.missingEst & exists("est", envir = .meta)) {
    est <- .meta$est
  }
  if (.missingEst & missing(control) & exists("control", envir = .meta)) {
    control <- .meta$control
    if (is(control, "foceiControl")) {
      est <- "focei"
      if (.missingEst && est != "focei") {
        warning(sprintf("Using focei instead of %s since focei controls were specified.", est))
      }
    } else if (is(control, "saemControl")) {
      est <- "saem"
      if (.missingEst && est != "saem") {
        warning(sprintf("Using saem instead of %s since saem controls were specified.", est))
      }
    }
  }
  if (missing(table) && exists("table", envir = .meta)) {
    table <- .meta$table
  }
  start.time <- Sys.time()
  if (!is(table, "tableControl")) {
    if (is(table, "list")) {
      table <- do.call(tableControl, table, envir = envir)
    } else {
      table <- tableControl()
    }
  }
  dat <- nlmixrData(data)
  nobs2 <- sum(dat$EVID == 0)
  up.covs <- toupper(uif$all.covs)
  up.names <- toupper(names(dat))
  for (i in seq_along(up.covs)) {
    w <- which(up.covs[i] == up.names)
    if (length(w) == 1) {
      names(dat)[w] <- uif$all.covs[i]
    }
  }
  ## backSort <- attr(dat, "backSort")
  ## backSort2 <- attr(dat, "backSort2")
  ## attr(dat, "backSort") <- NULL
  ## attr(dat, "backSort2") <- NULL
  if (!is.null(uif$nmodel$lin.solved) == 1L) {
    uif$env$infusion <- max(dat$EVID) > 10000
  }
  bad.focei <- "Problem calculating residuals, returning fit without residuals."
  fix.dat <- function(x) {
    .cls <- class(x)
    class(x) <- "data.frame"
    x$ID <- as.integer(x$ID)
    attr(x$ID, "levels") <- .lab
    class(x$ID) <- "factor"
    class(x) <- .cls
    .etaO <- x$etaObf
    .etaO$ID <- as.integer(.etaO$ID)
    attr(.etaO$ID, "levels") <- .lab
    class(.etaO$ID) <- "factor"
    .eta <- x$eta
    .eta$ID <- as.integer(.eta$ID)
    attr(.eta$ID, "levels") <- .lab
    class(.eta$ID) <- "factor"
    .ranef <- x$ranef
    .ranef$ID <- as.integer(.ranef$ID)
    attr(.ranef$ID, "levels") <- .lab
    class(.ranef$ID) <- "factor"
    .uif <- x$uif
    .thetas <- x$theta
    for (.n in names(.thetas)) {
      .uif$ini$est[.uif$ini$name == .n] <- .thetas[.n]
    }
    .omega <- x$omega
    for (.i in seq_along(.uif$ini$neta1)) {
      if (!is.na(.uif$ini$neta1[.i])) {
        .uif$ini$est[.i] <- .omega[.uif$ini$neta1[.i], .uif$ini$neta2[.i]]
      }
    }
    .env <- x$env
    .env$etaObf <- .etaO
    .env$eta <- .eta
    .env$origData <- .origData
    .env$uif <- .uif
    .env$ranef <- .ranef
    .predDf <- .uif$predDf
    if (any(.predDf$cond != "") & any(names(x) == "CMT")) {
      .cls <- class(x)
      class(x) <- "data.frame"
      x$CMT <- factor(x$CMT, levels = .predDf$cmt, labels = .predDf$cond)
      class(x) <- .cls
    }
    return(x)
  }
  .addNpde <- function(x) {
    .doIt <- table$npde
    if (is.null(.doIt)) {
      if (est == "saem") {
        .doIt <- table$saemNPDE
      } else if (est == "focei") {
        .doIt <- table$foceiNPDE
      } else if (est == "foce") {
        .doIt <- table$foceNPDE
      } else if (est == "nlme") {
        .doIt <- table$nlmeNPDE
      } else {
        .doIt <- FALSE
      }
    }
    if (!is.logical(.doIt)) {
      return(x)
    }
    if (is.na(.doIt)) {
      return(x)
    }
    if (!.doIt) {
      return(x)
    }
    .ret <- try(addNpde(x, nsim = table$nsim, ties = table$ties, seed = table$seed, updateObject = FALSE), silent = TRUE)
    if (inherits(.ret, "try-error")) {
      return(x)
    }
    return(.ret)
  }
  calc.resid <- table$cwres
  if (is.null(calc.resid)) {
    if (est == "saem") {
      calc.resid <- table$saemCWRES
    } else if (est == "nlme") {
      calc.resid <- table$nlmeCWRES
    }
  }
  if (est == "saem") {
    if (.nid <= 1) stop("SAEM is for mixed effects models, try 'dynmodel' (need more than 1 individual)")
    if (.nTv != 0) stop("SAEM does not support time-varying covariates (yet)")
    pt <- proc.time()
    if (length(uif$noMuEtas) > 0) {
      stop(sprintf("Cannot run SAEM since some of the parameters are not mu-referenced (%s)", paste(uif$noMuEtas, collapse = ", ")))
    }
    args <- as.list(match.call(expand.dots = TRUE))[-1]
    default <- saemControl()
    .getOpt <- function(arg, envir = parent.frame(1)) {
      if (arg %in% names(args)) {
        assign(paste0(".", arg), args[[arg]], envir = envir)
      } else if (arg %in% names(control)) {
        assign(paste0(".", arg), control[[arg]], envir = envir)
      } else if (arg %in% names(default)) {
        assign(paste0(".", arg), default[[arg]], envir = envir)
      }
    }
    for (a in c(
      "mcmc", "ODEopt", "seed", "print",
      "DEBUG", "covMethod", "calcTables",
      "logLik", "nnodes.gq",
      "nsd.gq", "nsd.gq", "adjObf",
      "optExpression"
    )) {
      .getOpt(a)
    }
    uif$env$optExpression <- .optExpression
    .addCov <- .covMethod == "linFim"

    if (uif$saemErr != "") {
      stop(paste0("For SAEM:\n", uif$saemErr))
    }
    if (is.null(uif$nlme.fun.mu)) {
      stop("SAEM requires all ETAS to be mu-referenced")
    }
    .err <- uif$ini$err
    .low <- uif$ini$lower
    .low <- .low[!is.na(.low) & is.na(.err)]
    .up <- uif$ini$upper
    .up <- .up[!is.na(.up) & is.na(.err)]
    if (any(.low != -Inf) | any(.up != Inf)) {
      warning("Bounds are ignored in SAEM")
    }
    uif$env$mcmc <- .mcmc
    uif$env$ODEopt <- .ODEopt
    uif$env$sum.prod <- sum.prod
    uif$env$covMethod <- .covMethod
    .dist <- uif$saem.distribution
    model <- uif$saem.model
    inits <- uif$saem.init
    if (length(uif$saem.fixed) > 0) {
      nphi <- attr(model$saem_mod, "nrhs")
      m <- cumsum(!is.na(matrix(inits$theta, byrow = TRUE, ncol = nphi)))
      fixid <- match(uif$saem.fixed, t(matrix(m, ncol = nphi)))

      names(inits$theta) <- rep("", length(inits$theta))
      names(inits$theta)[fixid] <- "FIXED"
    }
    .cfg <- configsaem(
      model = model, data = dat, inits = inits,
      mcmc = .mcmc, ODEopt = .ODEopt, seed = .seed,
      distribution = .dist, DEBUG = .DEBUG
    )
    if (is(.print, "numeric")) {
      .cfg$print <- as.integer(.print)
    }
    .fit <- model$saem_mod(.cfg)
    .ret <-
      try(
        as.focei.saemFit(
          .fit, uif, pt,
          data = dat, calcResid = calc.resid, obf = .logLik,
          nnodes.gq = .nnodes.gq, nsd.gq = .nsd.gq, adjObf = .adjObf,
          calcCov = .addCov, calcTables = .calcTables
        ),
        silent = TRUE
      )
    if (inherits(.ret, "try-error")) {
      warning("Error converting to nlmixr UI object, returning saem object")
      return(.fit)
    }
    if (inherits(.ret, "nlmixrFitData")) {
      .ret <- fix.dat(.ret)
      .ret <- .addNpde(.ret)
    }
    if (inherits(.ret, "nlmixrFitCore")) {
      .env <- .ret$env
      .env$adjObj <- .adjObf
      .env$nnodes.gq <- .nnodes.gq
      .env$nsd.gq <- .nsd.gq
      assign("startTime", start.time, .env)
      assign("est", est, .env)
      assign("stopTime", Sys.time(), .env)
      assign("origControl", control, .env)
      assign("modelId", .modelId, .env)
    }
    return(.ret)
  } else if (est == "nlme" || est == "nlme.mu" || est == "nlme.mu.cov" || est == "nlme.free") {
    if (.nid <= 1) stop("nlme is for mixed effects models, try 'dynmodel' (need more than 1 individual)")
    if (.nTv != 0) stop("nlme does not support time-varying covariates (yet)")
    data <- .as.data.frame(data)
    if (length(uif$predDf$cond) > 1) stop("nlmixr nlme does not support multiple endpoints.")
    pt <- proc.time()
    est.type <- est
    if (est == "nlme.free") {
      fun <- uif$nlme.fun
      specs <- uif$nlme.specs
    } else if (est == "nlme.mu") {
      fun <- uif$nlme.fun.mu
      specs <- uif$nlme.specs.mu
    } else if (est == "nlme.mu.cov") {
      fun <- uif$nlme.fun.mu.cov
      specs <- uif$nlme.specs.mu.cov
    } else {
      if (!is.null(uif$nlme.fun.mu.cov)) {
        est.type <- "nlme.mu.cov"
        fun <- uif$nlme.fun.mu.cov
        specs <- uif$nlme.specs.mu.cov
      } else if (!is.null(uif$nlme.fun.mu)) {
        est.type <- "nlme.mu"
        fun <- uif$nlme.fun.mu
        specs <- uif$nlme.specs.mu
      } else {
        est.type <- "nlme.free"
        fun <- uif$nlme.fun
        specs <- uif$nlme.fun.specs
      }
    }
    grp.fn <- uif$grp.fn
    dat$nlmixr.grp <-
      factor(apply(dat, 1, function(x) {
        cur <- x
        names(cur) <- names(dat)
        with(as.list(cur), {
          return(grp.fn())
        })
      }))
    dat$nlmixr.num <- seq_along(dat$nlmixr.grp)
    weight <- uif$nlme.var
    if (sum.prod) {
      rxode <- RxODE::rxSumProdModel(uif$rxode.pred)
    } else {
      rxode <- uif$rxode.pred
    }
    .atol <- 1e-8
    if (!is.null(control$atol)) .atol <- control$atol
    .rtol <- 1e-8
    if (!is.null(control$rtol)) .rtol <- control$rtol
    fit <- nlme_ode(dat,
      model = rxode,
      par_model = specs,
      par_trans = fun,
      response = "nlmixr_pred",
      weight = weight,
      verbose = TRUE,
      control = control,
      atol = .atol,
      rtol = .rtol,
      ...
    )
    class(fit) <- c(est.type, class(fit))
    .ret <- try({
      as.focei.nlmixrNlme(fit, uif, pt, data = dat, calcResid = calc.resid, nobs2 = nobs2)
    })
    if (inherits(.ret, "try-error")) {
      warning("Error converting to nlmixr UI object, returning nlme object")
      return(fit)
    }
    if (inherits(.ret, "nlmixrFitData")) {
      .ret <- fix.dat(.ret)
      .ret <- .addNpde(.ret)
    }
    if (inherits(.ret, "nlmixrFitCore")) {
      .env <- .ret$env
      assign("startTime", start.time, .env)
      assign("est", est, .env)
      assign("stopTime", Sys.time(), .env)
      assign("origControl", control, .env)
      assign("modelId", .modelId, .env)
    }
    return(.ret)
  } else if (any(est == c("foce", "focei", "fo", "foi"))) {
    if (.nid <= 1) stop(sprintf("%s is for mixed effects models, try 'dynmodel' (need more than 1 individual)", est))
    if (any(est == c("foce", "fo"))) {
      control$interaction <- FALSE
    }
    env <- new.env(parent = emptyenv())
    env$uif <- uif
    if (any(est == c("fo", "foi"))) {
      control$maxInnerIterations <- 0
      control$fo <- TRUE
      control$boundTol <- 0
      env$skipTable <- TRUE
    }
    fit <- foceiFit(dat,
      inits = uif$focei.inits,
      PKpars = uif$theta.pars,
      ## par_trans=fun,
      model = uif$rxode.pred,
      pred = function() {
        return(nlmixr_pred)
      },
      err = uif$error,
      lower = uif$focei.lower,
      upper = uif$focei.upper,
      fixed = uif$focei.fixed,
      thetaNames = uif$focei.names,
      etaNames = uif$eta.names,
      control = control,
      env = env,
      ...
    )
    if (any(est == c("fo", "foi"))) {
      ## Add posthoc.
      .default <- foceiControl()
      control$maxInnerIterations <- .default$maxInnerIterations
      control$maxOuterIterations <- 0L
      control$covMethod <- 0L
      control$fo <- 0L
      .uif <- fit$uif
      .thetas <- fit$theta
      for (.n in names(.thetas)) {
        .uif$ini$est[.uif$ini$name == .n] <- .thetas[.n]
      }
      .omega <- fit$omega
      for (.i in seq_along(.uif$ini$neta1)) {
        if (!is.na(.uif$ini$neta1[.i])) {
          .uif$ini$est[.i] <- .omega[.uif$ini$neta1[.i], .uif$ini$neta2[.i]]
        }
      }
      .time <- fit$time
      .objDf <- fit$objDf
      .message <- fit$env$message
      .time <- fit$time
      env <- new.env(parent = emptyenv())
      for (.w in c("cov", "covR", "covS", "covMethod")) {
        if (exists(.w, fit$env)) {
          assign(.w, get(.w, envir = fit$env), envir = env)
        }
      }
      env$time2 <- time
      env$uif <- .uif
      env$method <- "FO"
      fit0 <-
        try(
          foceiFit(
            dat,
            inits = .uif$focei.inits,
            PKpars = .uif$theta.pars,
            ## par_trans=fun,
            model = .uif$rxode.pred,
            pred = function() {
              return(nlmixr_pred)
            },
            err = .uif$error,
            lower = .uif$focei.lower,
            upper = .uif$focei.upper,
            fixed = .uif$focei.fixed,
            thetaNames = .uif$focei.names,
            etaNames = .uif$eta.names,
            control = control,
            env = env,
            ...
          ),
          silent = TRUE
        )
      if (inherits(fit0, "try-error")) {
      } else {
        fit <- fit0
        assign("message2", fit$env$message, env)
        assign("message", .message, env)
        .tmp1 <- env$objDf
        if (any(names(.objDf) == "Condition Number")) {
          .tmp1 <- .data.frame(.tmp1, "Condition Number" = NA, check.names = FALSE)
        }
        if (any(names(.tmp1) == "Condition Number")) {
          .objDf <- .data.frame(.objDf, "Condition Number" = NA, check.names = FALSE)
        }
        env$objDf <- rbind(.tmp1, .objDf)
        row.names(env$objDf) <- c(ifelse(est == "fo", "FOCE", "FOCEi"), "FO")
        .tmp1 <- env$time
        .tmp1$optimize <- .time$optimize
        .tmp1$covariance <- .time$covariance
        assign("time", .tmp1, envir = env)
        env$objDf <- env$objDf[, c("OBJF", "AIC", "BIC", "Log-likelihood", "Condition Number")]
        setOfv(env, "fo")
      }
    }
    fit <- .addNpde(fit)
    fit <- fix.dat(fit)
    assign("start.time", start.time, env)
    assign("est", est, env)
    assign("stop.time", Sys.time(), env)
    assign("origControl", control, env)
    assign("modelId", .modelId, env)
    return(fit)
  } else if (est == "posthoc") {
    if (.nid <= 1) stop("'posthoc' estimation is for mixed effects models, try 'dynmodel' (need more than 1 individual)")
    if (class(control) != "foceiControl") control <- do.call(nlmixr::foceiControl, control)
    control$covMethod <- 0L
    control$maxOuterIterations <- 0L
    .env <- new.env(parent = emptyenv())
    .env$uif <- uif
    fit <- foceiFit(dat,
      inits = uif$focei.inits,
      PKpars = uif$theta.pars,
      ## par_trans=fun,
      model = uif$rxode.pred,
      pred = function() {
        return(nlmixr_pred)
      },
      err = uif$error,
      lower = uif$focei.lower,
      upper = uif$focei.upper,
      thetaNames = uif$focei.names,
      etaNames = uif$eta.names,
      control = control,
      env = .env,
      ...
    )
    if (inherits(fit, "nlmixrFitData")) {
      .cls <- class(fit)
      .env <- attr(.cls, ".foceiEnv")
      .cls[1] <- "nlmixrPosthoc"
      class(fit) <- .cls
    }
    assign("uif", uif, fit$env)
    ## assign("start.time", start.time, env)
    ## assign("est", est, env)
    ## assign("stop.time", Sys.time(), env)
    ## assign("start.time", start.time, env)
    ## assign("est", est, env)
    ## assign("stop.time", Sys.time(), env)
    fit <- fix.dat(fit)
    assign("origControl", control, fit$env)
    assign("modelId", .modelId, fit$env)
    return(fit)
  } else if (est == "dynmodel") {
    if (class(control) != "dynmodelControl") control <- do.call(dynmodelControl, control)
    env <- new.env(parent = emptyenv())

    env$uif <- NULL

    # update data to merge for origData and data. first add zeros or whatever is filled in for DV when there is no observations
    # to match the lengths, then merge observed data for both origData and data, and send to RxODE.

    # .dynmodelData <- data
    # nlmixr Object ---
    .nmf <- uif
    # Conversion ---
    .dynNlmixr <- nlmixrDynmodelConvert(.nmf)
    # Model ---
    .system <- .dynNlmixr$system
    # Initial Estimates ---
    .inits <- .dynNlmixr$inits
    # Error Model ---
    .model <- .dynNlmixr$model
    # Optional Control ---
    control$nlmixrOutput <- T
    control$fixPars <- if (!is.null(.dynNlmixr$fixPars)) .dynNlmixr$fixPars else NULL
    control$lower <- if (!is.null(.dynNlmixr$lower)) .dynNlmixr$lower else NULL
    control$upper <- if (!is.null(.dynNlmixr$upper)) .dynNlmixr$upper else NULL

    fit <-
      dynmodel(
        system = .system,
        model = .model,
        inits = .inits,
        data = .origData,
        nlmixrObject = .nmf,
        control = control
      )
    assign("origData", .origData, fit$env)
    assign(".fit", fit, envir = .GlobalEnv)
    fit <- fix.dat(fit)
    return(fit)
  }
  else {
    stop(sprintf("Unknown estimation method est=\"%s\"", est))
  }
}

##' Fit a nlmixr model
##'
##' @param data Dataset to estimate.  Needs to be RxODE compatible (see
##'   \url{https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-event-types.html}
##'   for detailed dataset requirements).
##' @param uif Parsed nlmixr model (by \code{nlmixr(mod.fn)}).
##' @param est Estimation method
##' @param control Estimation control options.  They could be
##'   \code{\link[nlme]{nlmeControl}}, \code{\link{saemControl}} or
##'   \code{\link{foceiControl}}
##' @param ... Parameters passed to estimation method.

##' @param sum.prod Take the RxODE model and use more precise products/sums.
##'   Increases solving accuracy and solving time.
##' @param table A list controlling the table options (i.e. CWRES, NPDE etc).
##'   See \code{\link{tableControl}}.
##' @param save This option determines if the fit will be saved to be reloaded
##'   if already run.  If NULL, get the option from
##'   \code{options("nlmixr.save")};
##' @param envir Environment that nlmixr is evaluated in.
##' @return nlmixr fit object
##' @author Matthew L. Fidler
##' @export
nlmixr_fit <- function(uif, data, est = NULL, control = list(), ...,
                       sum.prod = FALSE, table = tableControl(),
                       save = NULL, envir = parent.frame()) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  if (is.null(save)) {
    save <- getOption("nlmixr.save", FALSE)
  }
  if (save) {
    .modName <- ifelse(is.null(uif$model.name), "", paste0(uif$model.name, "-"))
    if (.modName == ".-") .modName <- ""
    .dataName <- ifelse(is.null(uif$data.name), "", paste0(uif$data.name, "-"))
    if (.dataName == ".-") .dataName <- ""
    .digest <-
      digest::digest(
        list(
          gsub("<-", "=", gsub(" +", "", uif$fun.txt)),
          .as.data.frame(uif$ini),
          data,
          est,
          control,
          sum.prod,
          table,
          ...,
          as.character(utils::packageVersion("nlmixr")),
          as.character(utils::packageVersion("RxODE"))
        )
      )
    .saveFile <- file.path(
      getOption("nlmixr.save.dir", getwd()),
      paste0("nlmixr-", .modName, .dataName, est, "-", .digest, ".rds")
    )
    if (file.exists(.saveFile)) {
      message(sprintf("Loading model already run (%s)", .saveFile))
      .ret <- readRDS(.saveFile)
      if (!is.null(.ret$warnings)) {
        sapply(.ret$warnings, warning)
      }
      return(.ret)
    }
  }
  .ret <-
    .collectWarnings(
      nlmixr_fit0(
        uif = uif, data = data, est = est, control = control, ...,
        sum.prod = sum.prod, table = table, envir = envir
      ),
      TRUE
    )
  .ws <- .ret[[2]]
  .ret <- .ret[[1]]
  if (inherits(.ret, "nlmixrFitCore")) {
    .env <- .ret$env
    .env$warnings <- .ws
  }
  for (.i in seq_along(.ws)) {
    warning(.ws[.i])
  }
  if (inherits(.ret, "nlmixrFitCore")) {
    if (save) {
      AIC(.ret) # Calculate SAEM AIC when saving...
      .env <- .ret$env
      .extra <- (proc.time() - .nlmixrTime)["elapsed"] - sum(.env$time)
      .env$time <- .data.frame(.env$time, "other" = .extra, check.names = FALSE)
      saveRDS(.ret, file = .saveFile)
    } else {
      .env <- .ret$env
      .extra <- (proc.time() - .nlmixrTime)["elapsed"] - sum(.env$time)
      .env$time <- .data.frame(.env$time, "other" = .extra, check.names = FALSE)
    }
  }
  return(.ret)
}

##' Control Options for SAEM
##'
##' @param seed Random Seed for SAEM step.  (Needs to be set for
##'     reproducibility.)  By default this is 99.
##'
##' @param nBurn Number of iterations in the Stochastic Approximation
##'     (SA), or burn-in step. This is equivalent to Monolix's \code{K_0} or \code{K_b}.
##'
##' @param nEm Number of iterations in the Expectation-Maximization
##'     (EM) Step. This is equivalent to Monolix's \code{K_1}.
##'
##' @param nmc Number of Markov Chains. By default this is 3.  When
##'     you increase the number of chains the numerical integration by
##'     MC method will be more accurate at the cost of more
##'     computation.  In Monolix this is equivalent to \code{L}
##'
##' @param nu This is a vector of 3 integers. They represent the
##'     numbers of transitions of the three different kernels used in
##'     the Hasting-Metropolis algorithm.  The default value is \code{c(2,2,2)},
##'     representing 40 for each transition initially (each value is
##'     multiplied by 20).
##'
##'     The first value represents the initial number of multi-variate
##'     Gibbs samples are taken from a normal distribution.
##'
##'     The second value represents the number of uni-variate, or multi-
##'     dimensional random walk Gibbs samples are taken.
##'
##'     The third value represents the number of bootstrap/reshuffling or
##'     uni-dimensional random samples are taken.
##'
##' @param print The number it iterations that are completed before
##'     anything is printed to the console.  By default, this is 1.
##'
##' @param covMethod  Method for calculating covariance.  In this
##'     discussion, R is the Hessian matrix of the objective
##'     function. The S matrix is the sum of each individual's
##'     gradient cross-product (evaluated at the individual empirical
##'     Bayes estimates).
##'
##'  "\code{linFim}" Use the Linearized Fisher Information Matrix to calculate the covariance.
##'
##'  "\code{fim}" Use the SAEM-calculated Fisher Information Matrix to calculate the covariance.
##'
##'  "\code{r,s}" Uses the sandwich matrix to calculate the covariance, that is: \eqn{R^-1 \times S \times R^-1}
##'
##'  "\code{r}" Uses the Hessian matrix to calculate the covariance as \eqn{2\times R^-1}
##'
##'  "\code{s}" Uses the crossproduct matrix to calculate the covariance as \eqn{4\times S^-1}
##'
##'  "" Does not calculate the covariance step.
##'
##' @param logLik boolean indicating that log-likelihood should be
##'     calculate by Gaussian quadrature.
##'
##' @param trace An integer indicating if you want to trace(1) the
##'     SAEM algorithm process.  Useful for debugging, but not for
##'     typical fitting.
##'
##' @param nnodes.gq number of nodes to use for the Gaussian
##'     quadrature when computing the likelihood with this method
##'     (defaults to 1, equivalent to the Laplaclian likelihood)
##'
##' @param nsd.gq span (in SD) over which to integrate when computing
##'     the likelihood by Gaussian quadrature. Defaults to 3 (eg 3
##'     times the SD)
##'
##' @param adjObf is a boolean to indicate if the objective function
##'     should be adjusted to be closer to NONMEM's default objective
##'     function.  By default this is \code{TRUE}
##'
##' @param ... Other arguments to control SAEM.
##'
##' @inheritParams RxODE::rxSolve
##' @inheritParams foceiControl
##'
##' @return List of options to be used in \code{\link{nlmixr}} fit for
##'     SAEM.
##' @author Wenping Wang & Matthew L. Fidler
##' @export
saemControl <- function(seed = 99,
                        nBurn = 200, nEm = 300,
                        nmc = 3,
                        nu = c(2, 2, 2),
                        atol = 1e-06,
                        rtol = 1e-04,
                        method = "lsoda",
                        transitAbs = FALSE,
                        print = 1,
                        trace = 0,
                        covMethod = c("linFim", "fim", "r,s", "r", "s", ""),
                        calcTables = TRUE,
                        logLik = FALSE,
                        nnodes.gq = 3,
                        nsd.gq = 1.6,
                        optExpression = TRUE,
                        maxsteps = 100000L,
                        adjObf = TRUE,
                        ...) {
  .xtra <- list(...)
  .rm <- c()
  if (missing(transitAbs) && !is.null(.xtra$transit_abs)) {
    transitAbs <- .xtra$transit_abs
    .rm <- c(.rm, "transit_abs")
  }
  if (missing(nBurn) && !is.null(.xtra$n.burn)) {
    nBurn <- .xtra$n.burn
    .rm <- c(.rm, "n.burn")
  }
  if (missing(nEm) && !is.null(.xtra$n.em)) {
    nEm <- .xtra$n.em
    .rm <- c(.rm, "n.em")
  }
  .ret <- list(
    mcmc = list(niter = c(nBurn, nEm), nmc = nmc, nu = nu),
    ODEopt = RxODE::rxControl(
      atol = atol, rtol = rtol, method = method,
      transitAbs = transitAbs, maxsteps = maxsteps, ...
    ),
    seed = seed,
    print = print,
    DEBUG = trace,
    optExpression = optExpression,
    nnodes.gq = nnodes.gq,
    nsd.gq = nsd.gq,
    adjObf = adjObf,
    ...
  )
  if (length(.rm) > 0) {
    .ret <- .ret[!(names(.ret) %in% .rm)]
  }
  .ret[["covMethod"]] <- match.arg(covMethod)
  .ret[["logLik"]] <- logLik
  .ret[["calcTables"]] <- calcTables
  class(.ret) <- "saemControl"
  .ret
}

##' Add CWRES
##'
##' This returns a new fit object with CWRES attached
##'
##' @param fit nlmixr fit without WRES/CWRES
##' @param updateObject Boolean indicating if the original fit object
##'     should be updated. By default this is true.
##' @param envir Environment that should be checked for object to
##'     update.  By default this is the global environment.
##' @return fit with CWRES
##' @author Matthew L. Fidler
##' @export
addCwres <- function(fit, updateObject = TRUE, envir = globalenv()) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  .pt <- proc.time()
  .oTime <- fit$env$time
  .objName <- substitute(fit)
  if (any(names(fit) == "CWRES")) {
    ## warning("Already contains CWRES")
    return(fit)
  }
  .uif <- fit$uif
  .saem <- fit$saem
  .od <- fit$origData
  if (!is.null(.saem)) {
    assign("saem", NULL, fit$env)
    on.exit({
      assign("saem", .saem, fit$env)
    })
    .newFit <- as.focei.saemFit(.saem, .uif,
      data = .nmGetData(fit), calcResid = TRUE, obf = NA,
      calcCov = fit$cov, covMethod = fit$covMethod,
      calcCovTime = as.vector(fit$time[["covariance"]])
    )
    .ob1 <- .newFit$objDf
    .ob2 <- fit$objDf
    if (any(names(.ob2) == "Condition Number")) {
      .cn <- unique(.ob2[["Condition Number"]])
      .cn <- .cn[!is.na(.cn)]
      if (length(.cn) == 1) {
        .ob1[, "Condition Number"] <- .cn
      } else {
        .ob1[, "Condition Number"] <- NA
      }
    }
    .ob1 <- rbind(.ob1, .ob2)
    .ob1 <- .ob1[order(row.names(.ob1)), ]
    .ob1 <- .ob1[!is.na(.ob1$OBJF), ]
    assign("objDf", .ob1, envir = .newFit$env)
    assign("saem", .saem, fit$env)
    assign("adjObj", fit$env$adjObj, .newFit$env)
    .df <- .newFit[, c("WRES", "CRES", "CWRES", "CPRED")]
    ## Add CWRES timing to fit.
    .new <- cbind(fit, .df)
  }
  .nlme <- fit$nlme
  if (!is.null(.nlme)) {
    .newFit <- as.focei(.nlme, .uif, data = .nmGetData(fit), calcResid = TRUE)
    assign("adjObj", fit$env$adjObj, .newFit$env)
    .df <- .newFit[, c("WRES", "CRES", "CWRES", "CPRED")]
    .new <- cbind(fit, .df)
  }
  class(.new) <- class(.newFit)
  if (!is.null(.saem)) {
    setOfv(.new, "FOCEi")
  }
  if (updateObject) {
    .parent <- envir
    .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
      if (.cur == .objName && identical(.parent[[.cur]]$env, fit$env)) {
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
  .env <- .new$env
  .env$time <- .data.frame(.oTime, cwres = (proc.time() - .pt)["elapsed"], check.names = FALSE)
  assign("origData", .od, .env)
  return(.new)
}
