.setupPlotData <- function(data, cmt){
  .dat <- as.data.frame(data)
  .doCmt <- FALSE
  if (any(names(.dat) == "CMT")) {
    if (length(levels(.dat$CMT)) > 1) {
      .doCmt <- TRUE
    }
  }
  if (!.doCmt) {
    .dat$CMT <- factor(rep(1, length(.dat[, 1])), 1, "All Data")
  } else {
    levels(.dat$CMT) <- paste("Endpoint: ", levels(.dat$CMT))
  }
  if (any(names(.dat) == "CENS")) {
    .censLeft <- any(.dat$CENS == 1)
    .censRight <- any(.dat$CENS == -1)
    if (.censLeft & .censRight){
      .dat$CENS <- factor(.dat$CENS, c(-1, 0, 1), c("Right censored data", "Observed data", "left censored data"))
    } else if (.censLeft){
      .dat$CENS <- factor(.dat$CENS, c(0, 1), c("Observed data", "Censored data"))
    } else if (.censRight) {
      .dat$CENS <- factor(.dat$CENS, c(0, -1), c("Observed data", "Censored data"))
    } else {
      .dat <- .dat[, names(.dat) != "CENS"]
    }
  }
  return(.dat)
}

.dvPlot <- function(.dat0, vars){
  if (any(names(.dat0) == "CENS")){
    .data <- data.frame(DV = .dat0$DV, CENS=.dat0$CENS, stack(.dat0[, vars]))
    .aes <- ggplot2::aes(.data$values, .data$DV, color = .data$CENS)
  } else {
    .data <- data.frame(DV = .dat0$DV, stack(.dat0[, vars]))
    .aes <- ggplot2::aes(.data$values, .data$DV)
  }
  ggplot2::ggplot(.data, .aes) +
      ggplot2::facet_wrap(~ind) +
      ggplot2::geom_abline(slope = 1, intercept = 0, col = "red", size = 1.2) +
      ## ggplot2::geom_smooth(col="blue", lty=2, formula=DV ~ values + 0, size=1.2) +
      ggplot2::geom_point() +
      xlab("Predictions") +
      RxODE::rxTheme()
}


#' Plot a nlmixr data object
#'
#' Plot some standard goodness of fit plots for the focei fitted object
#'
#' @param x a focei fit object
#' @param ... additional arguments
#' @return NULL
#' @author Wenping Wang & Matthew Fidler
#' @export
plot.nlmixrFitData <- function(x, ...) {
  .lst <- list()
  object <- x
  IWRES <- NULL
  .tp <- traceplot(x)
  if (!is.null(.tp)) .lst[[length(.lst) + 1]] <- .tp
  .dat <- .setupPlotData(x)
  for (.cmt in levels(.dat$CMT)) {
    .dat0 <- .dat[.dat$CMT == .cmt, ]
    .hasCwres <- any(names(.dat0) == "CWRES")
    .hasNpde <- any(names(.dat0) == "NPDE")
    .p1 <- .dvPlot(.dat0, c("PRED", "IPRED")) +
      ggplot2::ggtitle(.cmt, "DV vs PRED/IPRED")
    .lst[[length(.lst) + 1]] <- .p1

    if (.hasCwres) {
      .d1 <- data.frame(DV = .dat0$DV, stack(.dat0[, c("CPRED", "IPRED")]))
      .p1 <- ggplot2::ggplot(.d1, ggplot2::aes_string("values", "DV")) +
        ggplot2::facet_wrap(~ind) +
        ggplot2::geom_abline(slope = 1, intercept = 0, col = "red", size = 1.2) +
        ## ggplot2::geom_smooth(col="blue", lty=2, formula=DV ~ values + 0, size=1.2) +
        ggplot2::geom_point() +
        xlab("Predictions") +
        ggplot2::ggtitle(.cmt, "DV vs CPRED/IPRED") +
        RxODE::rxTheme()
      .lst[[length(.lst) + 1]] <- .p1
    }

    if (.hasNpde) {
      .d1 <- data.frame(DV = .dat0$DV, stack(.dat0[, c("EPRED", "IPRED")]))
      .p1 <- ggplot2::ggplot(.d1, ggplot2::aes_string("values", "DV")) +
        ggplot2::facet_wrap(~ind) +
        ggplot2::geom_abline(slope = 1, intercept = 0, col = "red", size = 1.2) +
        ## ggplot2::geom_smooth(col="blue", lty=2, formula=DV ~ values + 0, size=1.2) +
        ggplot2::geom_point() +
        xlab("Predictions") +
        ggplot2::ggtitle(.cmt, "DV vs EPRED/IPRED") +
        RxODE::rxTheme()
      .lst[[length(.lst) + 1]] <- .p1
    }

    .p2 <- ggplot2::ggplot(.dat0, ggplot2::aes_string(x = "IPRED", y = "IRES")) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
      ggplot2::ggtitle(.cmt, "IRES vs IPRED") +
      RxODE::rxTheme()
    .lst[[length(.lst) + 1]] <- .p2

    .p2 <- ggplot2::ggplot(.dat0, ggplot2::aes_string(x = "TIME", y = "IRES")) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
      ggplot2::ggtitle(.cmt, "IRES vs TIME") +
      RxODE::rxTheme()
    .lst[[length(.lst) + 1]] <- .p2

    .p2 <- ggplot2::ggplot(.dat0, ggplot2::aes_string(x = "IPRED", y = "IWRES")) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
      ggplot2::ggtitle(.cmt, "IWRES vs IPRED") +
      RxODE::rxTheme()
    .lst[[length(.lst) + 1]] <- .p2

    .p2 <- ggplot2::ggplot(.dat0, ggplot2::aes_string(x = "TIME", y = "IWRES")) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
      ggplot2::ggtitle(.cmt, "IWRES vs IPRED")+
      RxODE::rxTheme()
    .lst[[length(.lst) + 1]] <- .p2

    if (.hasCwres) {
      .p2 <- ggplot2::ggplot(.dat0, ggplot2::aes_string(x = "CPRED", y = "CWRES")) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
        ggplot2::ggtitle(.cmt, "CWRES vs CPRED") +
        RxODE::rxTheme()
      .lst[[length(.lst) + 1]] <- .p2

      .p2 <- ggplot2::ggplot(.dat0, ggplot2::aes_string(x = "TIME", y = "CWRES")) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
        ggplot2::ggtitle(.cmt, "CWRES vs CPRED")+
        RxODE::rxTheme()
      .lst[[length(.lst) + 1]] <- .p2
    }
    if (.hasNpde) {
      .p2 <- ggplot2::ggplot(.dat0, ggplot2::aes_string(x = "EPRED", y = "NPDE")) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
        ggplot2::ggtitle(.cmt, "NPDE vs EPRED") +
        RxODE::rxTheme()
      .lst[[length(.lst) + 1]] <- .p2

      .p2 <- ggplot2::ggplot(.dat0, ggplot2::aes_string(x = "TIME", y = "NPDE")) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
        ggplot2::ggtitle(.cmt, "NPDE vs EPRED") +
        RxODE::rxTheme()
      .lst[[length(.lst) + 1]] <- .p2
    }
    ## .idPlot <- try(plot.nlmixrAugPred(nlmixrAugPred(object)));
    ## if (inherits(.idPlot, "try-error")){
    .ids <- unique(.dat0$ID)
    .s <- seq(1, length(.ids), by = 16)
    .j <- 0
    for (i in .s) {
      .j <- .j + 1
      .tmp <- .ids[seq(i, i + 15)]
      .tmp <- .tmp[!is.na(.tmp)]
      .d1 <- .dat0[.dat0$ID %in% .tmp, ]

      .p3 <- ggplot2::ggplot(.d1, aes(x = TIME, y = DV)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(aes(x = TIME, y = IPRED), col = "red", size = 1.2) +
        ggplot2::geom_line(aes(x = TIME, y = PRED), col = "blue", size = 1.2) +
        ggplot2::facet_wrap(~ID) +
        ggplot2::ggtitle(.cmt, sprintf("Individual Plots (%s of %s)", .j, length(.s))) +
        RxODE::rxTheme()
      .lst[[length(.lst) + 1]] <- .p3
    }
  }
  ## if (grDevices::dev.cur() != 1){
  ##     .x  <- .lst
  ##     for (.i in seq_along(.x)){
  ##         plot(.x[[.i]])
  ##     }
  ## }
  class(.lst) <- "nlmixrPlotList"
  return(.lst)
}

##' @export
plot.nlmixrPlotList <- function(x, y, ...) {
  .x <- x
  class(.x) <- NULL
  for (.i in seq_along(.x)) {
    plot(.x[[.i]])
  }
}

##' @export
plot.nlmixrFitCore <- function(x, ...) {
  stop("This is not a nlmixr data frame and cannot be plotted")
}

##' @export
plot.nlmixrFitCoreSilent <- plot.nlmixrFitCore


#' @title Produce trace-plot for fit if applicable
#'
#' @param x fit object
#' @param ... other parameters
#' @return Fit traceplot or nothing.
#' @author Rik Schoemaker, Wenping Wang & Matthew L. Fidler
#' @export
traceplot <- function(x, ...) {
  UseMethod("traceplot")
}

##' @rdname traceplot
##' @export
traceplot.nlmixrFitCore <- function(x, ...) {
  .m <- x$parHistStacked
  if (!is.null(.m)) {
    .p0 <- ggplot(.m, aes(iter, val)) +
      geom_line() +
      facet_wrap(~par, scales = "free_y")
    if (!is.null(x$mcmc)) {
      .p0 <- .p0 + ggplot2::geom_vline(xintercept = x$mcmc$niter[1], col = "blue", size = 1.2)
    }
    .p0 <- .p0 + RxODE::rxTheme()
    return(.p0)
  } else {
    return(invisible(NULL))
  }
}

##' @export
traceplot.nlmixrFitCoreSilent <- traceplot.nlmixrFitCore
