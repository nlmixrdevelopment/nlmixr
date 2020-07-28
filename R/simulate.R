## Add RxODE THETA/ETA replacement mini DSL
.repSim <- function(x, theta = c(), eta = c(), lhs = c()) {
  ret <- eval(parse(text = sprintf("quote({%s})", x)))
  f <- function(x) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      return(x)
    } else if (is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`[`))) {
        type <- tolower(as.character(x[[2]]))
        if (type == "theta") {
          return(eval(parse(text = sprintf("quote(%s)", theta[as.numeric(x[[3]])]))))
        } else if (type == "eta") {
          return(eval(parse(text = sprintf("quote(%s)", eta[as.numeric(x[[3]])]))))
        }
        stop("Only theta/eta translation supported.")
      } else if (length(x[[2]]) == 1 &&
        ((identical(x[[1]], quote(`=`))) ||
          (identical(x[[1]], quote(`~`))))) {
        if (any(as.character(x[[2]]) == lhs)) {
          if (any(as.character(x[[2]]) == c("rx_pred_", "rx_r_"))) {
            x[[1]] <- quote(`~`)
            return(as.call(lapply(x, f)))
          } else {
            return(NULL)
          }
        }
        return(as.call(lapply(x, f)))
      } else {
        return(as.call(lapply(x, f)))
      }
    } else {
      stop("Don't know how to handle type ", typeof(x),
        call. = FALSE
      )
    }
  }
  ret <- deparse(f(ret))[-1]
  ret <- ret[regexpr("^ *NULL$", ret) == -1]
  ret <- paste(ret[-length(ret)], collapse = "\n")
  return(RxODE::rxNorm(RxODE::rxGetModel(ret)))
}

.simInfo <- function(object) {
  .mod <- RxODE::rxNorm(object$model$pred.only)
  ## .mod <- gsub(rex::rex("d/dt(", capture(except_any_of("\n;)")), ")", or("=", "~")), "d/dt(\\1)~", .mod);
  .lhs <- object$model$pred.only$lhs
  .lhs <- .lhs[.lhs != "rx_pred_"]
  .lhs <- .lhs[.lhs != "rx_r_"]
  .omega <- object$omega
  .etaN <- dimnames(.omega)[[1]]
  .params <- nlme::fixed.effects(object)
  .thetaN <- names(.params)
  .newMod <- paste0(
    .repSim(.mod, theta = .thetaN, eta = .etaN, c(.lhs, "rx_pred_", "rx_r_")),
    "ipred=rxTBSi(rx_pred_, rx_lambda_, rx_yj_);"
  )
  ## paste0(gsub("\n\n+", "\n", gsub(rex::rex(capture(or(.lhs)), or("=", "~"), except_any_of("\n;"),one_of("\n;")), "",
  ##                                 gsub(rex::rex(capture(or("rx_pred_", "rx_r_")), or("=", "~")), "\\1~",
  ##                                      .repThetaEta(.mod, theta=.thetaN, eta=.etaN)))),
  ##        "ipred=rxTBSi(rx_pred_, rx_lambda_, rx_yj_);"));
  .sim <- "\nsim=rxTBSi(rx_pred_+rx_r_"
  .err <- object$uif$err
  .w <- which(!is.na(object$uif$err))
  .mat <- diag(length(.w))
  .dimn <- character(length(.w))
  for (.i in seq_along(.w)) {
    .ntheta <- object$uif$ini$ntheta[.w[.i]]
    .cur <- .thetaN[.ntheta]
    .dimn[.i] <- .cur
    if (any(.err[.w[.i]] == c("add", "norm", "dnorm", "lnorm", "dlnorm", "logn"))) {
      ## .sim <- paste0(.sim, "+", .cur);
      .mat[.i, .i] <- .params[.ntheta]^2
      .params[.ntheta] <- NA_real_
    } else if (.err[.w[.i]] == "prop") {
      ## .sim <- paste0(.sim, "+ipred*", .cur);
      .mat[.i, .i] <- .params[.ntheta]^2
      .params[.ntheta] <- NA_real_
    } else if (.err[.w[.i]] == "pow") {
      ## .sim <- paste0(.sim, "+", .cur, "*ipred^(", object$uif$ini$name[which(object$uif$ini$err == "pow2")], ")");
      .mat[.i, .i] <- .params[.ntheta]^2
      .params[.ntheta] <- NA_real_
    }
  }
  .params <- .params[!is.na(.params)]
  dimnames(.mat) <- list(.dimn, .dimn)
  .w <- which(!(.dimn %in% names(.params)))
  .mat <- .mat[.w, .w, drop = FALSE]
  .sigma <- .mat
  .sigmaNames <- dimnames(.mat)[[1]]
  .newMod <- paste0(.newMod, .sim, ", rx_lambda_, rx_yj_);\n")
  .newMod <- strsplit(.newMod, "\n")[[1]]
  .w <- which(regexpr("rx_r_~", .newMod) != -1)
  .subs <- function(x, .what, .with) {
    if (all(as.character(x) == .what)) {
      return(eval(parse(text = sprintf("quote(%s)", .with))))
    } else if (is.call(x)) {
      as.call(lapply(x, .subs, .what = .what, .with = .with))
    } else if (is.pairlist(x)) {
      as.pairlist(lapply(x, .subs, .what = .what, .with = .with))
    } else {
      return(x)
    }
  }
  for (.i in .w) {
    .cur <- sub(";", "", sub("rx_r_~", "", .newMod[.i]))
    .cur <- eval(parse(text = sprintf("RxODE::rxSplitPlusQ(quote(%s))", .cur)))
    .cur <- paste(sapply(.cur, function(x) {
      .ret <- sprintf("sqrt(%s)", x)
      for (.what in .sigmaNames) {
        .with <- "1"
        .old <- paste(deparse(eval(parse(text = sprintf("quote(%s)", .ret)))), collapse = "")
        .new <- paste(deparse(eval(parse(text = sprintf(
          ".subs(quote(%s),.what=%s,.with=%s)", .ret,
          deparse(.what), deparse(.with)
        )))), collapse = "")
        if (.old != .new) {
          .ret <- sprintf("%s*%s", .what, .new)
        }
      }
      return(.ret)
    }), collapse = "+")
    .cur <- gsub("  +", " ", .cur)
    .cur <- gsub(rex::rex("*sqrt(Rx_pow_di(1, 2))"), "", .cur)
    .cur <- gsub(rex::rex("Rx_pow_di(1, 2)", any_spaces, "*", any_spaces), "", .cur)
    .newMod[.i] <- paste0("rx_r_~", .cur, ";")
  }
  .newMod <- paste(paste(.newMod, collapse = "\n"), "\n")
  .dfObs <- object$nobs
  .nlmixrData <- nlmixr::nlmixrData(nlme::getData(object))
  .dfSub <- length(unique(.nlmixrData$ID))
  .thetaMat <- nlme::getVarCov(object)
  if (all(is.na(object$uif$ini$neta1))) {
    .omega <- NULL
    .dfSub <- 0
  }
  return(list(
    rx = .newMod, params = .params, events = .nlmixrData,
    thetaMat = .thetaMat, omega = .omega, sigma = .sigma, dfObs = .dfObs, dfSub = .dfSub
  ))
}


##' Simulate a nlmixr solved system
##'
##' This takes the uncertainty in the model parameter estimates and to
##' simulate a number of theoretical studies.  Each study simulates a
##' realization of the parameters from the uncertainty in the fixed
##' parameter estimates.  In addition the omega and sigma matrices are
##' simulated from the uncertainty in the Omega/Sigma matrices based
##' on the number of subjects and observations the model was based on.
##'
##' @param object nlmixr object
##' @param ... Other arguments sent to \code{rxSolve}
##' @inheritParams RxODE::rxSolve
##' @export
nlmixrSim <- function(object, ...) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit({
    RxODE::.setWarnIdSort(TRUE)
  })
  save <- getOption("nlmixr.save", FALSE)
  .si <- .simInfo(object)
  .xtra <- list(...)
  if (any(names(.xtra) == "rx")) {
    .si$rx <- .xtra$rx
  }
  if (!is.null(.xtra$modelName)) {
    message(sprintf("Compiling %s model...", .xtra$modelName), appendLF = FALSE)
  } else {
    message("Compiling model...", appendLF = FALSE)
  }
  .newobj <- RxODE::RxODE(.si$rx)
  on.exit({
    RxODE::rxUnload(.newobj)
  })
  message("done")
  if ((any(names(.xtra) == "nStud") && .xtra$nStud <= 1) || !any(names(.xtra) == "nStud")) {
    .si$thetaMat <- NULL
    .si$dfSub <- NULL
    .si$dfObs <- NULL
  } else {
    if (RxODE::rxIs(.xtra$thetaMat, "matrix")) {
      .si$thetaMat <- .xtra$thetaMat
    } else if (!is.null(.xtra$thetaMat)) {
      if (is.na(.xtra$thetaMat)) {
        .si$thetaMat <- NULL
      }
    }
    if (any(names(.xtra) == "dfSub")) {
      .si$dfSub <- .xtra$dfSub
    }
    if (any(names(.xtra) == "dfObs")) {
      .si$dfObs <- .xtra$dfObs
    }
  }


  if (any(names(.xtra) == "omega")) {
    .si$omega <- .xtra$omega
    if (any(is.na(.xtra$omega))) {
      .si$omega <- NULL
    }
  }
  if (any(names(.xtra) == "sigma")) {
    .si$sigma <- .xtra$sigma
    if (any(is.na(.xtra$sigma))) {
      .si$sigma <- NULL
    }
  }
  if (any(names(.xtra) == "events") &&
    RxODE::rxIs(.xtra$events, "rx.event")) {
    .si$events <- .xtra$events
  }
  if (any(names(.xtra) == "params")) {
    .si$params <- .xtra$params
  }
  .xtra$object <- .newobj
  .xtra$params <- .si$params
  .xtra$events <- .si$events
  if (RxODE::rxIs(.xtra$thetaMat, "matrix")) {
    .xtra$thetaMat <- NULL
  } else {
    .xtra$thetaMat <- .si$thetaMat
  }
  .xtra$dfObs <- .si$dfObs
  .xtra$omega <- .si$omega
  .xtra$dfSub <- .si$dfSub
  .xtra$sigma <- .si$sigma
  if (save) {
    .modName <- ifelse(is.null(object$uif$model.name), "", paste0(object$uif$model.name, "-"))
    if (.modName == ".-") .modName <- ""
    .dataName <- ifelse(is.null(object$uif$data.name), "", paste0(object$uif$data.name, "-"))
    if (.dataName == ".-") .dataName <- ""
    .digest <- digest::digest(list(
      gsub("<-", "=", gsub(" +", "", object$uif$fun.txt)),
      as.data.frame(object$uif$ini),
      .xtra,
      as.character(utils::packageVersion("nlmixr")),
      as.character(utils::packageVersion("RxODE"))
    ))
    .saveFile <- file.path(
      getOption("nlmixr.save.dir", getwd()),
      paste0("nlmixr-nlmixrSim-", .modName, .dataName, "-", .digest, ".rds")
    )
    if (file.exists(.saveFile)) {
      message(sprintf("Loading nlmixrSim already run (%s)", .saveFile))
      .ret <- readRDS(.saveFile)
      return(.ret)
    }
  }
  .ret <- do.call(getFromNamespace("rxSolve", "RxODE"), .xtra, envir = parent.frame(2))
  if (inherits(.ret, "rxSolve")) {
    .rxEnv <- attr(class(.ret), ".RxODE.env")
    if (!is.null(.xtra$nsim)) {
      .rxEnv$nSub <- .xtra$nsim
    }
    if (!is.null(.xtra$nSub)) {
      .rxEnv$nSub <- .xtra$nSub
    }
    if (is.null(.xtra$nStud)) {
      .rxEnv$nStud <- 1
    } else {
      .rxEnv$nStud <- .xtra$nStud
    }
    .cls <- c("nlmixrSim", class(.ret))
    attr(.cls, ".RxODE.env") <- .rxEnv
    if (any(names(.ret) == "CMT") && any(names(object) == "CMT")) {
      if (is(object$CMT, "factor")) {
        levels(.ret$CMT) <- levels(object$CMT)
        class(.ret$CMT) <- "factor"
      }
    }
    class(.ret) <- .cls
  }
  if (save) {
    saveRDS(.ret, file = .saveFile)
  }
  return(.ret)
}

##' @export
plot.nlmixrSim <- function(x, y, ...) {
  p1 <- eff <- Percentile <- sim.id <- id <- p2 <- p50 <- p05 <- p95 <- . <- NULL
  .args <- list(...)
  save <- getOption("nlmixr.save", FALSE)
  RxODE::rxReq("dplyr")
  RxODE::rxReq("tidyr")
  if (is.null(.args$p)) {
    .p <- c(0.05, 0.5, 0.95)
  } else {
    .p <- .args$p
  }
  if (save) {
    .digest <- digest::digest(list(
      .args,
      as.character(utils::packageVersion("nlmixr")),
      as.character(utils::packageVersion("RxODE"))
    ))
    .saveFile <- file.path(
      getOption("nlmixr.save.dir", getwd()),
      paste0("nlmixrSimPlot-", .digest, ".rds")
    )
    if (file.exists(.saveFile)) {
      message(sprintf("Loading nlmixrSimPlot already summarized (%s)", .saveFile))
      .ret <- readRDS(.saveFile)
      return(.ret)
    }
  }
  if (x$env$nStud <= 1) {
    if (x$env$nSub < 2500) {
      warning("In order to put confidence bands around the intervals, you need at least 2500 simulations.")
      message("Summarizing data for plot")
      .ret <- x %>%
        dplyr::group_by(time) %>%
        dplyr::do(data.frame(p1 = .p, eff = quantile(.$sim, probs = .p))) %>%
        dplyr::mutate(Percentile = factor(sprintf("%02d%%", round(p1 * 100))))
      message("done.")
      .ret <- ggplot2::ggplot(.ret, aes(time, eff, col = Percentile, fill = Percentile)) +
        ggplot2::geom_line(size = 1.2)
      return(.ret)
    } else {
      .n <- round(sqrt(x$env$nSub))
    }
  } else {
    .n <- x$env$nStud
  }
  message("Summarizing data for plot")
  .ret <- x %>%
    dplyr::mutate(id = sim.id %% .n) %>%
    dplyr::group_by(id, time) %>%
    dplyr::do(data.frame(p1 = .p, eff = quantile(.$sim, probs = .p))) %>%
    dplyr::group_by(p1, time) %>%
    dplyr::do(data.frame(p2 = .p, eff = quantile(.$eff, probs = .p))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p2 = sprintf("p%02d", (p2 * 100))) %>%
    tidyr::spread(p2, eff) %>%
    dplyr::mutate(Percentile = factor(sprintf("%02d%%", round(p1 * 100))))
  message("done.")
  .ret <- ggplot2::ggplot(.ret, aes(time, p50, col = Percentile, fill = Percentile)) +
    ggplot2::geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.5) +
    ggplot2::geom_line(size = 1.2)
  if (save) {
    saveRDS(.ret, file = .saveFile)
  }
  return(.ret)
}

## Mini DSL to fix pred-only models
.predOnlyRx <- function(object) {
  .ret <- eval(parse(text = sprintf("quote({%s})", RxODE::rxNorm(object$model$pred.only))))
  f <- function(x) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      return(x)
    } else if (is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (length(x[[2]]) == 1 &&
        ((identical(x[[1]], quote(`=`))) ||
          (identical(x[[1]], quote(`~`))))) {
        x[[1]] <- quote(`~`)
        return(as.call(lapply(x, f)))
      } else {
        if (((identical(x[[1]], quote(`=`))) ||
          (identical(x[[1]], quote(`~`)))) &&
          length(x[[2]] == 3)) {
          if (identical(x[[2]][[1]], quote(`/`))) {
            x[[1]] <- quote(`~`)
            return(as.call(lapply(x, f)))
          }
        }
        return(as.call(lapply(x, f)))
      }
    } else {
      stop("Don't know how to handle type ", typeof(x),
        call. = FALSE
      )
    }
  }
  .ret <- deparse(f(.ret))[-1]
  .ret <- .ret[regexpr("^ *NULL$", .ret) == -1]
  .ret <- .ret[-length(.ret)]
  .w <- rev(which(regexpr("^ *cmt[(].*[)] *$", .ret) != -1))
  .cur <- 1
  ## Remove cmt(#) at the beginning
  while (any(.w == .cur)) {
    .w <- .w[.w != .cur]
    .cur <- .cur + 1
  }
  ## Put rx_pred_ at the end, but before cmt(name) statements
  if (length(.w) > 2) {
    while (.w[1] == .w[2] + 1) {
      .w <- .w[-1]
      if (length(.w) <= 2) break
    }
  }
  if (length(.w) > 0) {
    .w <- .w[1]
    .ret[.w] <- paste0("pred=rx_pred_\n", .ret[.w])
  } else {
    .ret <- c(.ret, "pred=rx_pred_")
  }
  .ret <- paste(.ret, collapse = "\n")
  return(RxODE::RxODE(.ret))
}

##' Predict a nlmixr solved system
##'
##' @param ipred Flag to calculate individual predictions. When
##'     \code{ipred} is \code{TRUE}, calculate individual predictions.
##'     When \code{ipred} is \code{FALSE}, set calculate typical population predations.
##'     When \code{ipred} is \code{NA}, calculate both individual and
##'     population predictions.
##'
##' @inheritParams RxODE::rxSolve
##'
##' @export
nlmixrPred <- function(object, ..., ipred = FALSE) {
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  lst <- as.list(match.call()[-1])
  if (RxODE::rxIs(lst$params, "rx.event")) {
    if (!is.null(lst$events)) {
      tmp <- lst$events
      lst$events <- lst$params
      lst$params <- tmp
    } else {
      lst$events <- lst$params
      lst$params <- NULL
    }
  }
  if (!RxODE::rxIs(lst$events, "rx.event")) {
    lst$events <- nlmixrData(getData(object))
  }
  message("Compiling model...", appendLF = FALSE)
  lst$object <- .predOnlyRx(object)
  message("done")
  params <- fixed.effects(object)
  names(params) <- sprintf("THETA[%d]", seq_along(params))
  do.ipred <- FALSE
  do.pred <- FALSE
  if (is.na(ipred)) {
    do.ipred <- TRUE
    do.pred <- TRUE
  } else if (ipred) {
    do.ipred <- TRUE
  } else {
    do.pred <- TRUE
  }
  if (do.ipred) {
    re <- random.effects(object)
    if (is.null(re)) {
      .tmp <- lst$events
      .w <- which(tolower(names(.tmp)) == "id")
      .nid <- length(unique(.tmp[[.w]]))
      ipred.par <- data.frame(t(params),
        rx_err_ = rep(0, .nid),
        check.names = FALSE
      )
    } else {
      re <- re[, -1]
      names(re) <- sprintf("ETA[%d]", seq_along(names(re)))
      ipred.par <- data.frame(re, t(params),
        rx_err_ = 0, check.names = FALSE
      )
    }
  }
  if (do.pred) {
    neta <- dim(object$omega)[1]
    if (neta == 0) {
      pred.par <- c(params, rx_err_ = 0)
    } else {
      pred.par <- c(params, setNames(rep(0, neta + 1), c(sprintf("ETA[%d]", seq(1, neta)), "rx_err_")))
    }
  }
  on.exit(
    {
      RxODE::rxUnload(lst$object)
    },
    add = TRUE
  )
  if (!is.na(ipred)) {
    if (do.pred) {
      lst$params <- pred.par
    } else {
      lst$params <- ipred.par
    }
    ret <- suppressWarnings(do.call(getFromNamespace("rxSolve", "RxODE"), lst, envir = parent.frame(2)))
    if (do.ipred) {
      names(ret) <- sub("pred", "ipred", names(ret))
    }
    return(ret)
  } else {
    lst$params <- pred.par
    ret.pred <- suppressWarnings(do.call(getFromNamespace("rxSolve", "RxODE"), lst, envir = parent.frame(2)))
    lst$params <- ipred.par
    ret.pred$ipred <- suppressWarnings(do.call(getFromNamespace("rxSolve", "RxODE"), lst, envir = parent.frame(2)))$pred
    return(ret.pred)
  }
}
##' @rdname nlmixrPred
##' @export
predict.nlmixrFitData <- function(object, ...) {
  nlmixrPred(object, ...)
}

##' Augmented Prediction for nlmixr fit
##'
##'
##' @param object Nlmixr fit object
##' @inheritParams nlme::augPred
##' @inheritParams RxODE::rxSolve
##' @return Stacked data.frame with observations, individual/population predictions.
##' @author Matthew L. Fidler
##' @export
nlmixrAugPred <- function(object, ..., covsInterpolation = c("locf", "linear", "nocb", "midpoint"),
                          primary = NULL, minimum = NULL, maximum = NULL, length.out = 51L) {
  force(object)
  if (!inherits(object, "nlmixrFitData")) {
    stop("Need a nlmixr fit object")
  }
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  save <- getOption("nlmixr.save", FALSE)
  if (save) {
    .modName <- ifelse(is.null(object$uif$model.name), "", paste0(object$uif$model.name, "-"))
    if (.modName == ".-") .modName <- ""
    .dataName <- ifelse(is.null(object$uif$data.name), "", paste0(object$uif$data.name, "-"))
    if (.dataName == ".-") .dataName <- ""
    .digest <- digest::digest(list(
      gsub("<-", "=", gsub(" +", "", object$uif$fun.txt)),
      as.data.frame(object$uif$ini),
      covsInterpolation,
      primary, minimum, maximum, length.out,
      as.character(utils::packageVersion("nlmixr")),
      as.character(utils::packageVersion("RxODE"))
    ))
    .saveFile <- file.path(
      getOption("nlmixr.save.dir", getwd()),
      paste0("nlmixr-augPred-", .modName, .dataName, "-", .digest, ".rds")
    )
    if (file.exists(.saveFile)) {
      message(sprintf("Loading augPred already run (%s)", .saveFile))
      .ret <- readRDS(.saveFile)
      return(.ret)
    }
  }
  uif <- object$uif
  ## Using the model will drop the subjects that were dropped by the fit
  dat <- suppressWarnings(nlmixrData(getData(object), object$model$pred.only))
  dat <- as.data.frame(dat)
  names(dat) <- toupper(names(dat))
  attr(dat$ID, "levels") <- attr(object$ID, "levels")
  attr(dat$ID, "class") <- "factor"
  .isMulti <- (any(names(object) == "DVID") || any(names(object) == "CMT"))
  up.covs <- toupper(uif$all.covs)
  up.names <- toupper(names(dat))
  for (i in seq_along(up.covs)) {
    w <- which(up.covs[i] == up.names)
    if (length(w) == 1) {
      names(dat)[w] <- uif$all.covs[i]
    }
  }
  ids <- unique(dat$ID)
  .multiType <- NULL
  if (.isMulti) {
    stop("multiple endpoint augPred not supported yet.")
  }
  if (.isMulti) {
    if (any(names(dat) == "DVID")) {
      new.pts <- lapply(unique(object$uif$predDf$dvid), function(dvid) {
        .dat <- dat[dat$DVID == dvid, , drop = FALSE]
        r <- range(.dat$TIME, na.rm = TRUE, finite = TRUE)
        if (is.null(minimum) || is.infinite(minimum)) {
          minimum <- r[1]
        }
        if (is.null(maximum) || is.infinite(maximum)) {
          maximum <- r[2]
        }
        new.time <- seq(minimum, maximum, length.out = length.out)
        new.pts <- expand.grid(TIME = new.time, ID = ids, DVID = dvid)
      })
      new.pts <- do.call("rbind", new.pts)
      .multiType <- "DVID"
    } else {
      new.pts <- lapply(unique(object$uif$predDf$cmt), function(cmt) {
        .dat <- dat[dat$CMT == cmt, , drop = FALSE]
        r <- range(.dat$TIME, na.rm = TRUE, finite = TRUE)
        if (is.null(minimum) || is.infinite(minimum)) {
          minimum <- r[1]
        }
        if (is.null(maximum) || is.infinite(maximum)) {
          maximum <- r[2]
        }
        new.time <- seq(minimum, maximum, length.out = length.out)
        new.pts <- expand.grid(TIME = new.time, ID = ids, CMT = cmt)
      })
      new.pts <- do.call("rbind", new.pts)
      .multiType <- "CMT"
    }
  } else {
    r <- range(dat$TIME, na.rm = TRUE, finite = TRUE)
    if (is.null(minimum) || is.infinite(minimum)) {
      minimum <- r[1]
    }
    if (is.null(maximum) || is.infinite(maximum)) {
      maximum <- r[2]
    }
    new.time <- sort(unique(c(seq(minimum, maximum, length.out = length.out), dat$TIME)))
    new.pts <- expand.grid(TIME = new.time, ID = ids)
  }
  ## Add covariates in the augmented prediction
  covsi <- match.arg(covsInterpolation)
  all.covs <- uif$all.covs
  if (length(all.covs) > 0) {
    fs <- c(locf = 0, nocb = 1, midpoint = 0.5, linear = 0)
    new.cov <- lapply(all.covs, function(cov) {
      unlist(lapply(ids, function(id) {
        dat.id <- dat[dat$ID == id, ]
        fun <- stats::approxfun(dat.id$TIME, dat.id[[cov]],
          method = ifelse(covsi == "linear", "linear", "constant"),
          rule = 2,
          f = fs[covsi]
        )
        return(fun(new.time))
      }))
    })
    names(new.cov) <- all.covs
    new.pts <- cbind(new.pts, new.cov)
  }
  new.pts$EVID <- 0
  new.pts$AMT <- 0
  dat.old <- dat
  dat <- rbind(dat[, names(new.pts), drop = FALSE], new.pts)
  dat <- dat[order(dat$ID, dat$TIME), ]
  dat <- dat[!duplicated(paste(dat$ID, dat$TIME)), ]
  lst <- as.list(match.call()[-1])
  lst <- lst[!(names(lst) %in% c("primary", "minimum", "maximum", "length.out"))]
  lst$object <- object
  if (all(is.na(uif$ini$neta1))) {
    lst$ipred <- FALSE
  } else {
    lst$ipred <- NA
  }
  lst$events <- dat
  lst$params <- NULL
  if (.isMulti) {
    lst$keep <- .multiType
  }
  lst$returnType <- "data.frame.TBS"
  dat.new <- do.call("nlmixrPred", lst, envir = parent.frame(2))
  .Call(
    `_nlmixr_augPredTrans`, dat.new$pred, dat.new$ipred,
    dat.new$rxLambda, dat.new$rxYj
  )
  dat.new <- dat.new[, !(names(dat.new) %in% c("rxLambda", "rxYj"))]
  dat.new$id <- factor(dat.new$id)
  levels(dat.new$id) <- levels(object$ID)
  if (.isMulti) {
    if (.multiType == "DVID") {
      .tmp <- object$uif$nmodel$predDf[, c("cond", "dvid")]
      names(.tmp) <- c("Endpoint", "DVID")
      dat.new <- merge(dat.new, .tmp, by = "DVID")
      dat.new <- dat.new[, names(dat.new) != "DVID", drop = FALSE]
      .endpoint <- dat.new[, "Endpoint"]
      dat.new <- dat.new[, names(dat.new) != "Endpoint", drop = FALSE]
      dat.new <- data.frame(dat.new[, 1:2], Endpoint = .endpoint, stack(dat.new[, -(1:2), drop = FALSE]))
      levels(dat.new$ind) <- gsub("pred", "Population", gsub("ipred", "Individual", levels(dat.new$ind)))
      dat.new$Endpoint <- factor(dat.new$Endpoint)
    } else {
      .tmp <- object$uif$nmodel$predDf[, c("cond", "cmt")]
      names(.tmp) <- c("Endpoint", "CMT")
      dat.new <- merge(dat.new, .tmp, by = "CMT")
      dat.new <- dat.new[, names(dat.new) != "CMT", drop = FALSE]
      .endpoint <- dat.new[, "Endpoint"]
      dat.new <- dat.new[, names(dat.new) != "Endpoint", drop = FALSE]
      dat.new <- data.frame(dat.new[, 1:2], Endpoint = .endpoint, stack(dat.new[, -(1:2), drop = FALSE]))
      levels(dat.new$ind) <- gsub("pred", "Population", gsub("ipred", "Individual", levels(dat.new$ind)))
      dat.new$Endpoint <- factor(dat.new$Endpoint)
    }
  } else {
    dat.new <- data.frame(dat.new[, 1:2], stack(dat.new[, -(1:2), drop = FALSE]))
    levels(dat.new$ind) <- gsub("pred", "Population", gsub("ipred", "Individual", levels(dat.new$ind)))
  }
  names(dat.old) <- tolower(names(dat.old))
  dat.old <- dat.old[dat.old$evid == 0, ]
  if (.isMulti) {
    if (.multiType == "DVID") {
      ## FIXME
      if (inherits(dat.old$dvid, "factor") || inherits(dat.old$dvid, "character")) {
        dat.old <- data.frame(id = dat.old$id, time = dat.old$time, Endpoint = dat.old$dvid, values = dat.old$dv, ind = "Observed")
      } else {
        .tmp <- object$uif$nmodel$predDf[, c("cond", "dvid")]
        names(.tmp) <- c("Endpoint", "DVID")
        dat.old <- merge(dat.old, .tmp, by = "DVID")
        dat.old <- data.frame(
          id = dat.old$id, time = dat.old$time, Endpoint = dat.old$Endpoint,
          values = dat.old$dv, ind = "Observed"
        )
      }
    } else {
      .tmp <- object$uif$nmodel$predDf[, c("cond", "cmt")]
      names(.tmp) <- c("Endpoint", "CMT")
      .low <- tolower(names(dat.old))
      .w <- which(.low == "cmt")
      if (length(.w) == 1) {
        names(dat.old)[.w] <- "CMT"
      }
      dat.old <- merge(dat.old, .tmp, by = "CMT")
      dat.old <- data.frame(
        id = dat.old$id, time = dat.old$time, Endpoint = dat.old$Endpoint,
        values = dat.old$dv, ind = "Observed"
      )
    }
  } else {
    dat.old <- data.frame(id = dat.old$id, time = dat.old$time, values = dat.old$dv, ind = "Observed")
  }
  .ret <- rbind(dat.new, dat.old)
  if (save) {
    saveRDS(.ret, file = .saveFile)
  }
  return(.ret)
}

##' @rdname nlmixrAugPred
##' @export
augPred.nlmixrFitData <- memoise::memoise(function(object, primary = NULL, minimum = NULL, maximum = NULL,
                                                   length.out = 51, ...) {
  .ret <- nlmixrAugPred(
    object = object, primary = primary,
    minimum = minimum, maximum = maximum,
    length.out = length.out, ...
  )
  class(.ret) <- c("nlmixrAugPred", "data.frame")
  return(.ret)
})

##' @export
plot.nlmixrAugPred <- function(x, y, ...) {
  if (any(names(x) == "Endpoint")) {
    for (.tmp in unique(x$Endpoint)) {
      .x <- x[x$Endpoint == .tmp, names(x) != "Endpoint"]
      plot.nlmixrAugPred(.x)
    }
  } else {
    ids <- unique(x$id)
    time <- values <- ind <- id <- NULL # Rcheck fix
    for (i in seq(1, length(ids), by = 16)) {
      tmp <- ids[seq(i, i + 15)]
      tmp <- tmp[!is.na(tmp)]
      d1 <- x[x$id %in% tmp, ]
      dobs <- d1[d1$ind == "Observed", ]
      dpred <- d1[d1$ind != "Observed", ]
      p3 <- ggplot(d1, aes(time, values, col = ind)) +
        geom_line(data = dpred, size = 1.2) +
        geom_point(data = dobs) +
        facet_wrap(~id)
      print(p3)
    }
  }
}

##' @rdname nlmixrSim
##' @export
rxSolve.nlmixrFitData <- function(object, params = NULL, events = NULL, inits = NULL, scale = NULL,
                                  method = c("liblsoda", "lsoda", "dop853"),
                                  transitAbs = NULL, atol = 1.0e-6, rtol = 1.0e-4,
                                  maxsteps = 5000L, hmin = 0L, hmax = NULL, hini = 0L, maxordn = 12L, maxords = 5L, ...,
                                  cores, covsInterpolation = c("linear", "locf", "nocb", "midpoint"),
                                  addCov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL,
                                  nCoresRV = 1L, sigmaIsChol = FALSE, nDisplayProgress = 10000L,
                                  amountUnits = NA_character_, timeUnits = "hours", stiff,
                                  theta = NULL, eta = NULL, addDosing = FALSE, updateObject = FALSE,
                                  omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,
                                  nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,
                                  nStud = 1L, dfSub = 0.0, dfObs = 0.0, returnType = c("rxSolve", "matrix", "data.frame"),
                                  seed = NULL, nsim = NULL) {
  do.call("nlmixrSim", as.list(match.call()[-1]), envir = parent.frame(2))
}

##' @rdname nlmixrSim
##' @export
simulate.nlmixrFitData <- function(object, nsim = 1, seed = NULL, ...) {
  nlmixr::nlmixrSim(object, ..., nsim = nsim, seed = seed)
}

##' @rdname nlmixrSim
##' @export
solve.nlmixrFitData <- function(a, b, ...) {
  lst <- as.list(match.call()[-1])
  n <- names(lst)
  if (!missing(a)) {
    n[n == "a"] <- ""
  }
  if (!missing(b)) {
    n[n == "b"] <- ""
  }
  names(lst) <- n
  do.call("nlmixrSim", lst, envir = parent.frame(2))
}

##' @importFrom RxODE rxSolve
##' @export
RxODE::rxSolve
