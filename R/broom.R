.fixNames <- function(.df) {
  names(.df) <- gsub("Back-.*", "estimate", names(.df))
  names(.df) <- gsub("Estimate", "model.est", names(.df))
  names(.df) <- gsub(".*RSE.*", "rse", names(.df))
  names(.df) <- gsub("SE", "std.error", names(.df))
  names(.df) <- gsub("CI Lower", "conf.low", names(.df))
  names(.df) <- gsub("CI Upper", "conf.high", names(.df))
  names(.df) <- gsub("BSV.*", "bsv", names(.df))
  names(.df) <- gsub("Shrink.*", "shrink", names(.df))
  return(.df)
}

##' @export
confint.nlmixrFitCore <- function(object, parm, level = 0.95, ...) {
  .extra <- list(...)
  .exponentiate <- ifelse(any(names(.extra) == "exponentiate"), .extra$exponentiate, FALSE)
  .ciNames <- ifelse(any(names(.extra) == "ciNames"), .extra$ciNames, TRUE)
  .df <- .fixNames(object$parFixedDf)
  if (!missing(parm)) .df <- .df[parm, ]
  .zv <- qnorm(1 - (1 - level) / 2)
  .low <- .df$model.est - .df$std.error * .zv
  .hi <- .df$model.est + .df$std.error * .zv
  if (is.na(.exponentiate)) {
    .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
    .hi[.exp] <- exp(.hi[.exp])
    .low[.exp] <- exp(.low[.exp])
  } else if (.exponentiate) {
    .hi <- exp(.hi)
    .low <- exp(.low)
  }
  .df <- data.frame(model.est = .df$model.est, estimate = .df$estimate, conf.low = .low, conf.high = .hi)
  if (.ciNames) names(.df)[3:4] <- paste(c((1 - level) / 2, (1 - (1 - level) / 2)) * 100, "%")
  .df
}

##' @export
confint.nlmixrFitCoreSilent <- confint.nlmixrFitCore

#' @importFrom tibble as_tibble
.nlmixrTidyFixed <- function(x, ..., .ranpar = FALSE) {
  .extra <- list(...)
  .conf.int <- ifelse(any(names(.extra) == "conf.int"), .extra$conf.int, ifelse(any(names(.extra) == "conf.level"), TRUE, FALSE))
  .conf.level <- ifelse(any(names(.extra) == "conf.level"), .extra$conf.level, 0.95)
  .exponentiate <- ifelse(any(names(.extra) == "exponentiate"), .extra$exponentiate, FALSE)
  .quick <- ifelse(any(names(.extra) == "quick"), .extra$quick, FALSE)
  .rse <- ifelse(any(names(.extra) == "rse"), .extra$rse, FALSE)
  .bsv <- ifelse(any(names(.extra) == "bsv"), .extra$bsv, FALSE)
  .shrink <- ifelse(any(names(.extra) == "shrink"), .extra$shrink, FALSE)
  if (.quick) warning("quick does not do anything for nlmixr fit objects")
  .df <- .fixNames(x$parFixedDf)
  .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
  if (is.na(.exponentiate)) {
    ## se(exp(x)) ~= exp(mu_x)*se_x
    ##
    ## Norris N. The Standard Errors of the Geometric and Harmonic Means and Their Application to Index Numbers
    ## Ann. Math. Statist. Volume 11, Number 4 (1940), 445-448.
    .df$std.error[.exp] <- exp(.df$model.est[.exp]) * .df$std.error[.exp]
  } else if (.exponentiate) {
    .df$std.error <- exp(.df$model.est) * .df$std.error
    .df$estimate <- exp(.df$model.est)
  } else if (!.exponentiate) {
    .df$estimate <- .df$model.est
  }
  .df$statistic <- .df$estimate / .df$std.error
  .df$p.value <- stats::pt(.df$statistic, nobs(x) - attr(logLik(x), "df"), lower.tail = FALSE)
  if (.conf.int) {
    .ci <- confint.nlmixrFitCore(x, level = .conf.level, ciNames = FALSE, exponentiate = .exponentiate)
    .df$conf.low <- .ci$conf.low
    .df$conf.high <- .ci$conf.high
  } else {
    .df <- .df[, regexpr("^conf[.]", names(.df)) == -1]
  }
  if (!.rse) {
    .df <- .df[, names(.df) != "rse"]
  }
  if (!.bsv) {
    .df <- .df[, names(.df) != "bsv"]
  }
  if (!.shrink) {
    .df <- .df[, names(.df) != "shrink"]
  }
  .ini <- x$uif$ini$err[!is.na(x$uif$ini$ntheta)]
  ## effect   group   term            estimate std.error statistic
  .df <- data.frame(
    effect = "fixed",
    term = row.names(.df), .df, stringsAsFactors = FALSE
  )
  if (!.ranpar) {
    .df <- .df[is.na(x$uif$ini$err[!is.na(x$uif$ini$ntheta)]), ]
  } else {
    .tmp <- x$uif$ini$err[!is.na(x$uif$ini$ntheta)]
    .df <- .df[!is.na(.tmp), ]
    .tmp <- .tmp[!is.na(.tmp)]
    .df$group <- paste0("Residual(", .tmp, ")")
    .df$effect <- "ran_pars"
  }
  tibble::as_tibble(.df)
}

#' @importFrom tibble as_tibble
.nlmixrTidyRandom <- function(x, ...) {
  .d <- dim(x$omegaR)
  if (.d[1] > 0) {
    .tmp <- stack(x$eta[, -1])
    .df <- data.frame(group = "ID", level = x$eta$ID, term = .tmp$ind, estimate = .tmp$values)
    return(tibble::as_tibble(.df))
  } else {
    return(NULL)
  }
}

#' @importFrom tibble as_tibble
.nlmixrTidyRandomPar <- function(x, ...) {
  .pars <- .getR(x$omegaR, TRUE)
  if (length(.pars) > 0) {
    .p1 <- data.frame(
      effect = "ran_pars", group = "ID", term = names(.pars), estimate = .pars, std.error = NA_real_,
      statistic = NA_real_, p.value = NA_real_, stringsAsFactors = FALSE
    ) %>%
      .reorderCols()
    .p2 <- data.frame(.nlmixrTidyFixed(x, .ranpar = TRUE), stringsAsFactors = FALSE) %>%
      .reorderCols()
    .df <- rbind(.p1, .p2)
    for (.v in c("statistic", "p.value", "std.error")) {
      if (all(is.na(.df[, .v]))) {
        .df <- .df[, names(.df) != .v]
      }
    }
    return(tibble::as_tibble(.df))
  } else {
    return(NULL)
  }
}
##   effect   group    term                  estimate std.error statistic
##   <chr>    <chr>    <chr>                    <dbl>     <dbl>     <dbl>
## 1 fixed    NA       (Intercept)           251.          6.82     36.8
## 2 fixed    NA       Days                   10.5         1.55      6.77
## 3 ran_pars Subject  sd__(Intercept)        24.7        NA        NA
## 4 ran_pars Subject  sd__Days                5.92       NA        NA
## 5 ran_pars Subject  cor__(Intercept).Days   0.0656     NA        NA
## 6 ran_pars Residual sd__Observation        25.6        NA        NA

## Row names and order taken & adapted from
## https://github.com/bbolker/broom.mixed/blob/master/R/utilities.R#L238-L248
.reorderCols <- function(x) {
  allCols <- c(
    "response", "effect",
    "component", ## glmmTMB, brms
    "group", "level", "term", "index", "estimate",
    "std.error", "statistic",
    "df", "p.value",
    "conf.low", "conf.high", "rhat", "ess"
  )
  return(x[, intersect(allCols, names(x))])
}

#' @importFrom tibble as_tibble
.coefPar <- function(x, exponentiate = FALSE, ...) {
  .d <- dim(x$omegaR)
  if (.d[1] == 0) {
    return(NULL)
  }
  .muRef <- x$mu.ref
  .theta <- x$theta
  .df <- .fixNames(x$parFixedDf)
  if (is.na(exponentiate)) {
    .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
  } else if (exponentiate) {
    .exp <- setNames(rep(TRUE, length(.theta)), names(.theta))
  } else {
    .exp <- setNames(rep(FALSE, length(.theta)), names(.theta))
  }
  .eta <- x$eta
  .noMuRef <- c()
  .x <- setNames(
    data.frame(lapply(names(.eta), function(n) {
      if (any(n == names(.muRef))) {
        .n <- .muRef[[n]]
        .ret <- .eta[[n]] + .theta[.n]
        if (any(.n == names(.exp))) {
          if (.exp[.n]) {
            .ret <- exp(.ret)
          }
        }
      } else {
        .ret <- .eta[[n]]
        if (n != "ID") {
          warning(sprintf("The parameter '%s' is not mu-referenced and the coef will not be returned", n))
          .noMuRef <<- c(.noMuRef, n)
        }
      }
      return(.ret)
    })), sapply(names(.eta), function(n) {
      if (any(n == names(.muRef))) {
        .n <- .muRef[[n]]
        return(.n)
      }
      return(n)
    })
  )
  .tmp <- stack(.x[, -1])
  .df <- data.frame(group = "ID", level = .x$ID, term = .tmp$ind, estimate = .tmp$values)
  return(tibble::as_tibble(.df))
}

#' @export
#' @importFrom tibble as_tibble
tidy.nlmixrFitCore <- function(x, ...) {
  .extra <- list(...)
  if (any(names(.extra) == "effects")) {
    .effects <- .extra$effects
  } else {
    if (any(names(.extra) == "effect")) {
      .effects <- .extra$effect
    } else {
      .effects <- c("fixed", "ran_pars")
    }
  }
  .effects <- match.arg(.effects, c("fixed", "random", "ran_vals", "ran_pars", "ran_coef"),
    several.ok = TRUE
  )
  .ret <- list()
  if (any(.effects == "fixed")) {
    .ret$fixed <- .nlmixrTidyFixed(x, ...)
  }
  if (any(.effects == "random") || any(.effects == "ran_vals")) {
    .ret$ran_vals <- .nlmixrTidyRandom(x, ...)
  }
  if (any(.effects == "ran_pars")) {
    .ret$ran_pars <- .nlmixrTidyRandomPar(x, ...)
  }
  if (any(.effects == "ran_coef")) {
    .ret$ran_coef <- .coefPar(x, ...)
  }
  if (all(unlist(lapply(seq_along(.ret), is.null)))) {
    return(NULL)
  }
  return(dplyr::bind_rows(.ret, .id = "effect") %>%
    tibble::as_tibble() %>%
    .reorderCols())
}

##' @export
tidy.nlmixrFitCoreSilent <- tidy.nlmixrFitCore

#' @export
#' @importFrom tibble as_tibble
glance.nlmixrFitCore <- function(x, ...) {
  .lst <- list(...)
  if (any(names(.lst) == "type")) {
    setOfv(x, type = .lst$type)
  }
  .aic <- AIC(x) ## To calculate AIC if needed
  .df <- x$objDf[x$objDf$AIC == .aic, ]
  names(.df) <- gsub("Log-likelihood", "logLik", names(.df))
  names(.df) <- gsub("Condition Number", "conditionNumber", names(.df))
  tibble::as_tibble(.df)
}

##' @export
glance.nlmixrFitCoreSilent <- glance.nlmixrFitCore

##' @export
augment.nlmixrFitCore <- function(x, ...) {
  stop("augment is not yet implemented for nlmixr models")
}

##' @export
augment.nlmixrFitCoreSilent <- augment.nlmixrFitCore
