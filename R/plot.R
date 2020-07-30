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

.dvPlot <- function(.dat0, vars, log=FALSE){
  .xgxr <- getOption("RxODE.xgxr", TRUE) &&
    requireNamespace("xgxr", quietly = TRUE)
  if (any(names(.dat0) == "CENS")){
    .data <- data.frame(DV = .dat0$DV, CENS=.dat0$CENS, stack(.dat0[, vars]))
    .aes <- ggplot2::aes(.data$values, .data$DV, color = .data$CENS)
    if (length(levels(.data$CENS)) == 3){
      .color <- ggplot2::scale_color_manual(values=c("blue", "black", "red"))
    } else {
      .color <- ggplot2::scale_color_manual(values=c("black", "red"))
    }
    .legendPos <- ggplot2::theme(legend.position="bottom",legend.box="horizontal",
                                 legend.title = ggplot2::element_blank())
  } else {
    .data <- data.frame(DV = .dat0$DV, stack(.dat0[, vars]))
    .aes <- ggplot2::aes(.data$values, .data$DV)
    .color <- NULL
    .legendPos <- NULL
  }
  .logx <- NULL
  .logy <- NULL
  if (log){
    if (.xgxr) {
      .logx <- xgxr::xgx_scale_x_log10()
      .logy <- xgxr::xgx_scale_y_log10()
    } else {
      .logx <- ggplot2::scale_x_log10()
      .logy <- ggplot2::scale_y_log10()
    }
  }
  ggplot2::ggplot(.data, .aes) +
    ggplot2::facet_wrap(~ind) +
    ggplot2::geom_abline(slope = 1, intercept = 0, col = "red", size = 1.2) +
    .logx + .logy +
    ## ggplot2::geom_smooth(col="blue", lty=2, formula=DV ~ values + 0, size=1.2) +
    ggplot2::geom_point(alpha=0.5) +
    xlab("Predictions") +
    RxODE::rxTheme() + .color + .legendPos
}

.scatterPlot <- function(.dat0, vars, .cmt, log=FALSE){
  .data <- .dat0
  .data$x <- .data[[vars[1]]]
  .data$y <- .data[[vars[2]]]
  if (any(names(.dat0) == "CENS")){
    .aes <- ggplot2::aes(.data$x, .data$y, color = .data$CENS)
    if (length(levels(.data$CENS)) == 3){
      .color <- ggplot2::scale_color_manual(values=c("blue", "black", "red"))
    } else {
      .color <- ggplot2::scale_color_manual(values=c("black", "red"))
    }
    .legendPos <- ggplot2::theme(legend.position="bottom",legend.box="horizontal",
                                 legend.title = ggplot2::element_blank())
  } else {
    .aes <- ggplot2::aes(.data$x, .data$y)
    .aes <- ggplot2::aes(.data$values, .data$DV)
    .color <- NULL
    .legendPos <- NULL
  }
  .xgxr <- getOption("RxODE.xgxr", TRUE) &&
    requireNamespace("xgxr", quietly = TRUE)
  .logx <- NULL
  if (log){
    if (.xgxr) {
      .logx <- xgxr::xgx_scale_x_log10()
    } else {
      .logx <- ggplot2::scale_x_log10()
    }
  }
  ggplot2::ggplot(.data, .aes) +
    ggplot2::geom_point(alpha=0.5) +
    ggplot2::geom_abline(slope = 0, intercept = 0, col = "red") +
    ggplot2::ggtitle(.cmt, paste0(vars[1], " vs ", vars[2])) +
    ggplot2::xlab(vars[1]) + ggplot2::ylab(vars[2]) +
    RxODE::rxTheme() + .color + .legendPos + .logx
}

#' @title Produce trace-plot for fit if applicable
#'
#' @param x fit object
#' @param ... other parameters
#' @return Fit traceplot or nothing.
#' @author Vipul Mann,  Matthew L. Fidler
#' @export
bootplot <- function(x, ...){
  UseMethod("bootplot")
}
##' @rdname traceplot
##' @export
bootplot.nlmixrFitCore <- function(x, ...) {
  if (inherits(x, "nlmixrFitCore")) {
    if (exists(".bootPlotData", x$env)){
      .chisq <- x$env$.bootPlotData$chisq
      .dfD <- x$env$.bootPlotData$dfD
      .deltaN <- x$env$.bootPlotData$deltaN
      .df2 <- x$env$.bootPlotData$df2
      .plot <- ggplot2::ggplot(.chisq, aes(quantiles, deltaofv, color=Distribution)) +
        ggplot2::geom_line() + ggplot2::ylab("\u0394 objective function") +
        ggplot2::geom_text(data=.dfD, aes(label=label), hjust=0) +
        ggplot2::xlab("Distribution quantiles") +
        ggplot2::scale_color_manual(values=c("red", "blue")) +
        RxODE::rxTheme() +
        ggplot2::theme(legend.position="bottom",legend.box="horizontal")

      if (requireNamespace("ggtext", quietly = TRUE)) {
        .plot <- .plot +
          ggplot2::theme(plot.title = ggtext::element_markdown(),
                         legend.position="none") +
          ggplot2::labs(
            title = paste0(
              'Bootstrap <span style="color:blue; opacity: 0.2;">\u0394 objective function (', .deltaN,
              ' models, df\u2248', .df2, ')</span> vs <span style="color:red; opacity: 0.2;">reference \u03C7\u00B2(df=',
              length(fit$ini$est), ")</style>"
            ),
            caption = "\u0394 objective function curve should be on or below the reference distribution curve"
          )
      } else {
        .plot <- ggplot2::labs(
          title = paste0("Distribution of \u0394 objective function values for ", .deltaN, ' df=', .df2, " models"),
          caption = "\u0394 objective function curve should be on or below the reference distribution curve"
        )
      }
      .plot
    } else {
      stop("this nlmixr object does not include boostrap distribution statics for comparison",
           call.=FALSE)
    }
  } else {
    stop("this is not a nlmixr object",
         call.=FALSE)
  }
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
  if (exists(".bootPlotData", object$env)){
    .bp <- bootplot(x)
    .lst[[length(.lst) + 1]] <- .bp
  }
  .dat <- .setupPlotData(x)
  for (.cmt in levels(.dat$CMT)) {
    .dat0 <- .dat[.dat$CMT == .cmt, ]
    .hasCwres <- any(names(.dat0) == "CWRES")
    .hasNpde <- any(names(.dat0) == "NPDE")
    .p1 <- .dvPlot(.dat0, c("PRED", "IPRED")) +
      ggplot2::ggtitle(.cmt, "DV vs PRED/IPRED")
    .lst[[length(.lst) + 1]] <- .p1

    .p1 <- .dvPlot(.dat0, c("PRED", "IPRED"), TRUE) +
      ggplot2::ggtitle(.cmt, "log-scale DV vs PRED/IPRED")
    .lst[[length(.lst) + 1]] <- .p1

    if (.hasCwres) {
      .p1 <- .dvPlot(.dat0, c("CPRED", "IPRED")) +
        ggplot2::ggtitle(.cmt, "DV vs CPRED/IPRED")
      .lst[[length(.lst) + 1]] <- .p1

      .p1 <- .dvPlot(.dat0, c("CPRED", "IPRED"), TRUE) +
        ggplot2::ggtitle(.cmt, "log-scale DV vs CPRED/IPRED")
      .lst[[length(.lst) + 1]] <- .p1
    }

    if (.hasNpde) {
      .p1 <- .dvPlot(.dat0, c("EPRED", "IPRED"))  +
        ggplot2::ggtitle(.cmt, "DV vs EPRED/IPRED")
      .lst[[length(.lst) + 1]] <- .p1

      .p1 <- .dvPlot(.dat0, c("EPRED", "IPRED"), TRUE)  +
        ggplot2::ggtitle(.cmt, "log-scale DV vs EPRED/IPRED")
      .lst[[length(.lst) + 1]] <- .p1
    }

    for (x in c("IPRED", "PRED", "CPRED", "EPRED", "TIME")) {
      if (any(names(.dat0) == x)) {
        for (y in c("IWRES", "IRES", "RES", "CWRES", "NPDE")) {
          if (any(names(.dat0) == y)) {
            if (y == "CWRES" && x %in% c("TIME", "CPRED")){
              .doIt <- TRUE
            } else if (y == "NPDE" && x %in% c("TIME", "EPRED")) {
              .doIt <- TRUE
            } else if (!(y %in% c("CWRES", "NPDE"))) {
              .doIt <- TRUE
            }
            if (.doIt){
              .p2 <- .scatterPlot(.dat0, c(x, y), .cmt, log=FALSE)
              .lst[[length(.lst) + 1]] <- .p2
              .p2 <- .scatterPlot(.dat0, c(x, y), .cmt, log=TRUE)
              .lst[[length(.lst) + 1]] <- .p2
            }

          }
        }
      }
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
      if (any(names(.d1) == "lowerLim")) {
        .p3 <- .p3 + geom_cens(aes(lower=lowerLim, upper=upperLim), fill="purple")
      }
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
