##' VPC based on ui model
##'
##' @param fit nlmixr fit object
##' @param data this is the data to use to augment the VPC fit.  By
##'     default is the fitted data, (can be retrieved by
##'     \code{\link[nlme]{getData}}), but it can be changed by specifying
##'     this argument.
##' @param n Number of VPC simulations.  By default 100
##' @inheritParams vpc::vpc
##' @inheritParams RxODE::rxSolve
##' @param ... Args sent to \code{\link[RxODE]{rxSolve}}
##' @return Simulated dataset (invisibly)
##' @author Matthew L. Fidler
##' @export
vpc_ui <- function(fit, data = NULL, n = 100, bins = "jenks",
                   n_bins = "auto", bin_mid = "mean",
                   show = NULL, stratify = NULL, pred_corr = FALSE,
                   pred_corr_lower_bnd = 0, pi = c(0.05, 0.95), ci = c(0.05, 0.95),
                   uloq = NULL, lloq = NULL, log_y = FALSE, log_y_min = 0.001,
                   xlab = NULL, ylab = NULL, title = NULL, smooth = TRUE, vpc_theme = NULL,
                   facet = "wrap", labeller = NULL, vpcdb = FALSE, verbose = FALSE, ...) {
  if (!requireNamespace("vpc", quietly = TRUE)) {
    stop("'vpc' is required; Download from CRAN or github https://github.com/ronkeizer/vpc")
  }
  RxODE::.setWarnIdSort(FALSE)
  on.exit(RxODE::.setWarnIdSort(TRUE))
  save <- getOption("nlmixr.save", FALSE)
  if (save) {
    .modName <- ifelse(is.null(fit$uif$model.name), "", paste0(fit$uif$model.name, "-"))
    if (.modName == ".-") .modName <- ""
    .dataName <- ifelse(is.null(fit$uif$data.name), "", paste0(fit$uif$data.name, "-"))
    if (.dataName == ".-") .dataName <- ""
    .digest <- digest::digest(list(
      gsub("<-", "=", gsub(" +", "", fit$uif$fun.txt)),
      .as.data.frame(fit$uif$ini),
      data, n, bins,
      n_bins,
      bin_mid,
      show, stratify, pred_corr,
      pred_corr_lower_bnd,
      pi, ci,
      uloq, lloq, log_y, log_y_min,
      xlab, ylab, title, smooth, vpc_theme,
      facet, labeller, vpcdb, verbose,
      as.character(utils::packageVersion("nlmixr")),
      as.character(utils::packageVersion("RxODE"))
    ))
    .saveFile <- file.path(
      getOption("nlmixr.save.dir", getwd()),
      paste0("nlmixr-vpc-", .modName, .dataName, "-", .digest, ".rds")
    )
    if (file.exists(.saveFile)) {
      message(sprintf("Loading vpc already run (%s)", .saveFile))
      .ret <- readRDS(.saveFile)
      return(.ret)
    }
  }
  if (is(data, "numeric") | is(data, "integer")) {
    if (missing(n)) {
      n <- data
      data <- NULL
    } else {
      stop("Data needs to be a data frame instead of a numeric value.")
    }
  }
  .xtra <- list(...)
  if (inherits(fit, "nlmixrVpc")) {
    sim <- attr(class(fit), "nlmixrVpc")
  } else {
    .xtra$object <- fit
    ## .xtra$returnType <- "data.frame";
    .xtra$returnType <- "rxSolve"
    .xtra$modelName <- "VPC"
    pt <- proc.time()
    .dt <- fit$origData
    .si <- fit$simInfo
    if (is.null(data)) {
      if (!is.null(stratify)) {
        .rx <- paste0(.si$rx, "\nrx_dummy_var~", paste(stratify, collapse = "+"), "\n")
      } else {
        .rx <- .si$rx
      }
      dat <- nlmixrData(.dt, .rx)
    } else {
      dat <- data
    }
    ## Get stratify columns
    if (!is.null(stratify)) {
      cols <- c(stratify, "dv")
      stratify <- tolower(stratify)
    } else if (inherits(fit, "nlmixrFitCore")) {
      .uif <- fit$uif
      if (length(.uif$predDf$cmt) > 1) {
        cols <- c("cmt", "dv")
        stratify <- "cmt"
      } else {
        cols <- c("dv")
      }
    } else {
      cols <- c("dv")
    }
    .nd <- names(dat)
    .xtra$rx <- paste(c(do.call("c", lapply(cols, function(x) {
      if (any(tolower(x) == c("dv"))) {
        return(NULL)
      }
      .w <- which(tolower(.nd) == tolower(x))
      if (length(.w) == 1L) {
        if (tolower(x) == "cmt") {
          names(dat)[.w] <<- paste0("nlmixr_", tolower(.nd[.w]))
          return(paste0("cmt0=0+nlmixr_", tolower(.nd[.w])))
        } else if (tolower(x) == "dvid") {
          names(dat)[.w] <<- paste0("nlmixr_", tolower(.nd[.w]))
          return(paste0("dvid0=0+nlmixr_", tolower(.nd[.w])))
        } else {
          names(dat)[.w] <<- paste0("nlmixr_", tolower(.nd[.w]))
          return(paste0(.nd[.w], "=0+nlmixr_", tolower(.nd[.w])))
        }
      }
      return(NULL)
    })), .si$rx), collapse = "\n")
    .xtra$nStud <- n
    if (!is.null(.xtra$nsim)) {
      .xtra$nStud <- .xtra$nsim
      .xtra$nsim <- NULL
    }
    .xtra$dfObs <- 0
    .xtra$dfSub <- 0
    .xtra$thetaMat <- NA
    names(dat) <- sapply(names(dat), function(x) {
      if (any(tolower(x) == cols)) {
        return(tolower(x))
      }
      return(x)
    })
    if (any(names(dat) == "nlmixr_cmt")) dat$cmt <- dat$nlmixr_cmt
    if (any(names(dat) == "nlmixr_dvid")) dat$dvid <- dat$nlmixr_dvid
    .xtra$events <- dat
    sim <- do.call("nlmixrSim", .xtra)
    max.resim <- 10
    while (any(is.na(sim$sim))) {
      sim <- sim[!is.na(sim$sim), ]
      sim$id <- as.integer(factor(paste(sim$id)))
      print(summary(sim))
    }
    if (any(names(sim) == "cmt0")) names(sim)[names(sim) == "cmt0"] <- "cmt"
    if (any(names(sim) == "dvid0")) names(sim)[names(sim) == "dvid0"] <- "dvid"
    sim0 <- sim
    if (any(names(dat) == "nlmixr_cmt")) dat <- dat[, names(dat) != "nlmixr_cmt"]
    if (any(names(dat) == "nlmixr_dvid")) dat <- dat[, names(dat) != "nlmixr_dvid"]
    names(dat) <- gsub("nlmixr_", "", names(dat))
    onames <- names(dat)
    names(dat) <- tolower(onames)
    w <- which(duplicated(names(dat)))
    if (length(w) > 0) {
      warning(sprintf("Dropping duplicate columns (case insensitive): %s", paste(onames, collapse = ", ")))
      dat <- dat[, -w, drop = FALSE]
    }
    names(sim) <- gsub("^sim$", "dv", tolower(names(sim)))
    sim0 <- sim
    ##
    if (pred_corr) {
      .xtra.prd <- .xtra
      .xtra.prd$modelName <- "Pred (for pcVpc)"
      .xtra.prd$params <- c(
        .si$params, setNames(rep(0, dim(.si$omega)[1]), dimnames(.si$omega)[[2]]),
        setNames(rep(0, dim(.si$sigma)[1]), dimnames(.si$sigma)[[2]])
      )
      .xtra.prd$omega <- NA
      .xtra.prd$sigma <- NA
      .xtra.prd$returnType <- "data.frame"
      .xtra.prd$nStud <- 1
      .xtra.prd$nsim <- NULL
      sim2 <- do.call("nlmixrSim", .xtra.prd)
      sim$pred <- sim2$sim
    }
    diff <- proc.time() - pt
    message(sprintf("done (%.2f sec)", diff["elapsed"]))

    ##
    ## Assume this is in the observed dataset. Add it to the current dataset
    dat <- dat[dat$evid == 0, ]
    if (any(names(sim) == "cmt") && any(names(fit) == "CMT")) {
      if (is(fit$CMT, "factor")) {
        sim$cmt <- factor(sim$cmt, sort(unique(sim$cmt)), labels = levels(fit$CMT))
      }
    }
    if (any(names(dat) == "cmt") && any(names(fit) == "CMT")) {
      if (is(fit$CMT, "factor")) {
        dat$cmt <- factor(dat$cmt, sort(unique(dat$cmt)), labels = levels(fit$CMT))
      }
    }
    sim <- list(rxsim = sim0, sim = .as.data.frame(sim), obs = dat)
    class(sim) <- "rxHidden"
    attr(sim, "nsim") <- .xtra$nsim
    class(sim) <- "nlmixrVpc"
  }
  ns <- loadNamespace("vpc")
  if (exists("vpc_vpc", ns)) {
    vpcn <- "vpc_vpc"
  } else {
    vpcn <- "vpc"
  }
  call <- as.list(match.call(expand.dots = TRUE))[-1]
  call <- call[names(call) %in% methods::formalArgs(getFromNamespace(vpcn, "vpc"))]
  call$obs_cols <- list(id = "id", dv = "dv", idv = "time")
  call$sim_cols <- list(id = "id", dv = "dv", idv = "time")
  call$stratify <- stratify
  p <- do.call(getFromNamespace(vpcn, "vpc"), c(sim, call), envir = parent.frame(1))
  cls <- c("nlmixrVpc", class(p))
  attr(cls, "nlmixrVpc") <- sim
  class(p) <- cls
  if (save) {
    saveRDS(p, file = .saveFile)
  }
  return(p)
}

##' @export
`$.nlmixrVpc` <- function(obj, arg, exact = TRUE) {
  if (arg == "gg") {
    .x <- obj
    .cls <- class(.x)[class(.x) != "nlmixrVpc"]
    attr(.cls, "nlmixrVpc") <- NULL
    class(.x) <- .cls
    return(.x)
  } else if (any(arg == c("rxsim", "sim", "obs"))) {
    .info <- attr(class(obj), "nlmixrVpc")
    return(.info[[arg]])
  }
  NextMethod()
}

##' @export
print.nlmixrVpc <- function(x, ...) {
  cat(sprintf(
    "nlmixr vpc object of %d simulations.\n",
    attr(attr(class(x), "nlmixrVpc"), "nsim")
  ))
  cat("  $rxsim = original simulated data\n")
  cat("  $sim = merge simulated data\n")
  cat("  $obs = observed data\n")
  cat("  $gg = vpc ggplot\n")
  cat("use vpc(...) to change plot options\n")
  cat("plotting the object now\n")
  NextMethod()
}

##' @rdname vpc_ui
##' @export
vpc.nlmixrFitData <- function(sim, ...) {
  vpc_ui(fit = sim, ...)
}

##' @export
vpc.nlmixrFOCEi <- vpc.nlmixrFitData

##' @export
vpc.nlmixrSaem <- vpc.nlmixrFitData

##' @export
vpc.nlmixrNlme <- vpc.nlmixrFitData

##' @export
vpc.nlmixrPosthoc <- vpc.nlmixrFitData

##' @rdname vpc_ui
##' @export
vpc.nlmixrVpc <- function(sim, ...) {
  vpc_ui(fit = sim, ...)
}

##' @rdname vpc_ui
##' @export vpc.ui
##' @rawNamespace S3method(vpc,ui)
vpc.ui <- function(sim, ...) {
  vpc_ui(fit = sim, ...)
}
