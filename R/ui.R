
nlmixrfindLhs <- function(x) {
  ## Modified from http://adv-r.had.co.nz/Expressions.html find_assign4
  if (is.atomic(x) || is.name(x)) {
    character()
  } else if (is.call(x)) {
    if ((identical(x[[1]], quote(`<-`)) ||
      identical(x[[1]], quote(`=`))) &&
      is.name(x[[2]])) {
      lhs <- as.character(x[[2]])
    } else {
      lhs <- character()
    }
    unique(c(lhs, unlist(lapply(x, nlmixrfindLhs))))
  } else if (is.pairlist(x)) {
    unique(unlist(lapply(x, nlmixrfindLhs)))
  } else {
    stop("Don't know how to handle type ", typeof(x),
      call. = FALSE
    )
  }
}

.deparse <- function(expr) {
  deparse(expr, width.cutoff = 500, control = "useSource")
}

.bodyDewrap <- function(ret) {
  .ret <- ret
  if (length(.ret) > 1) {
    .ret[1] <- sub(rex::rex(
      start, or(group("function", any_spaces, "(", any_spaces, ")", any_spaces), ""),
      any_spaces, any_of("("), any_spaces, "{", any_spaces
    ), "", .ret[1], perl = TRUE)
    if (.ret[1] == "") .ret <- .ret[-1]
  }
  if (length(.ret) > 1) {
    .len <- length(.ret)
    .ret[.len] <- sub(rex::rex(any_spaces, "}", any_spaces, any_of(")"), any_spaces, end), "", .ret[.len])
    if (.ret[.len] == "") .ret <- .ret[-.len]
  }
  return(.ret)
}

.deparse1 <- function(expr) {
  return(.bodyDewrap(deparse(expr, width.cutoff = 500)))
}

.pipedData <- NULL
.clearPipedData <- function() {
  assignInMyNamespace(".pipedData", NULL)
}

.getPipedData <- function() {
  .pipedData
}

.getUif <- function(what) {
  if (inherits(what, "nlmixrFitCore")) {
    .uif <- what$uif
    assignInMyNamespace(".pipedData", getData(what))
  } else if (inherits(what, "nlmixrUI")) {
    .uif <- what
  } else if (inherits(what, "function")) {
    .uif <- nlmixrUI(what)
  } else {
    stop("Do not know how to handle object")
  }
  return(.uif)
}

.processLotri <- function(.code, .uif) {
  .mat <- try(eval(parse(text = sprintf("lotri::lotri(%s)", .code))), silent = TRUE)
  if (inherits(.mat, "matrix")) {
    .d <- dimnames(.mat)[[1]]
    .ini2 <- .as.data.frame(.uif$ini)
    .ini1 <- .ini2[is.na(.ini2$neta1), ]
    .ini2 <- .ini2[!is.na(.ini2$neta1), ]
    .ini4 <- .ini2[!is.na(.ini2$err), ]
    .ini2 <- .ini2[is.na(.ini2$err), ]
    .m2 <- as.vector(.mat)
    .l2 <- length(.d) * 2
    .df <- .data.frame(n1 = rep(.d, each = length(.d)), n2 = rep(.d, length(.d)), val = .m2)
    .dfI <- .as.data.frame(.uif$ini[, c("neta1", "neta2", "name")])
    .dfI <- .dfI[!is.na(.dfI$neta1), ]
    .dfI <- .dfI[which(.dfI$neta1 == .dfI$neta2), c("neta1", "name")]
    .d2 <- paste(.dfI$name)
    .diff <- setdiff(.d, .d2)
    if (length(.diff) > 0) {
      stop(sprintf(
        "trying to provide an estimate for non-existant eta: %s",
        paste(.diff, collapse = ", ")
      ))
    }
    names(.dfI)[2] <- "n1"
    .df$n1 <- paste(.df$n1)
    .dfI$n1 <- paste(.dfI$n1)
    .df <- merge(.df, .dfI, all.x = TRUE)
    names(.dfI) <- c("neta2", "n2")
    .df <- merge(.df, .dfI, all.x = TRUE)
    .df <- .df[order(.df$neta1, .df$neta2), ]
    ## Name becomes (eta.cl,eta.ka)
    .df$name <- ifelse(.df$neta1 == .df$neta2, .df$n1, paste0("(", .df$n1, ",", .df$n2, ")"))
    .df$name2 <- ifelse(.df$neta1 == .df$neta2, .df$n1, paste0("(", .df$n2, ",", .df$n1, ")"))
    .ini3 <- do.call("rbind", lapply(seq_along(.df$name), function(i) {
      .name1 <- paste(.df$name[i])
      .name2 <- paste(.df$name2[i])
      .w <- which(.name1 == .ini2$name)
      if (length(.w) == 1) {
        .ini2$est[.w] <<- .df$val[i]
        return(NULL)
      }
      .w <- which(.ini2$name == .name2)
      if (length(.w) == 1) {
        .ini2$est[.w] <<- .df$val[i]
        return(NULL)
      }
      if (.df$neta1[i] < .df$neta2[i]) {
        return(NULL)
      }
      .df2 <- nlmixrBoundsTemplate
      .df2$neta1 <- .df$neta1[i]
      .df2$neta2 <- .df$neta2[i]
      .df$name <- .df$name[i]
      .df2$est <- .df$val[i]
      .df2$lower <- -Inf
      .df2$upper <- Inf
      .df2$condition <- "ID"
      return(.df2)
    }))
    .ini2 <- rbind(.ini2, .ini3)
    .ini2 <- .ini2[order(.ini2$neta1, .ini2$neta2), ]
    .ini2 <- rbind(.ini1, .ini2, .ini4)
    class(.ini2) <- c("nlmixrBounds", "data.frame")
    .uif$ini <- .ini2
    return(.uif)
  } else {
    stop(sprintf("invalid syntax: %s", .code))
  }
}

#' nlmixr ini block handling
#'
#' The ini block controls initial conditions for 'theta' (fixed effects),
#' 'omega' (random effects), and 'sigma' (residual error) elements of the model.
#'
#' 'theta' and 'sigma' can be set using either \code{<-} or \code{=} such as
#' \code{tvCL <- 1} or equivalently \code{tvCL = 1}.  'omega' can be set with a
#' \code{~}.
#'
#' Parameters can be named or unnamed (though named parameters are preferred).
#' A named parameter is set using the name on the left of the assignment while
#' unnamed parameters are set without an assignment operator.  \code{tvCL <- 1}
#' would set a named parameter of \code{tvCL} to \code{1}.  Unnamed parameters
#' are set using just the value, such as \code{1}.
#'
#' For some estimation methods, lower and upper bounds can be set for 'theta'
#' and 'sigma' values.  To set a lower and/or upper bound, use a vector of
#' values.  The vector is \code{c(lower, estimate, upper)}.  The vector may be
#' given with just the estimate (\code{c(estimate)}), the lower bound and
#' estimate (\code{c(lower, estimate)}), or all three (\code{c(lower, estimate,
#' upper)}).  To set an estimate and upper bound without a lower bound, set the
#' lower bound to \code{-Inf}, \code{c(-Inf, estimate, upper)}.  When an
#' estimation method does not support bounds, the bounds will be ignored with a
#' warning.
#'
#' 'omega' values can be set as a single value or as the values of a
#' lower-triangular matrix.  The values may be set as either a
#' variance-covariance matrix (the default) or as a correlation matrix for the
#' off-diagonals with the standard deviations on the diagonals.  Names may be
#' set on the left side of the \code{~}.  To set a variance-covariance matrix
#' with variance values of 2 and 3 and a covariance of -2.5 use \code{~c(2, 2.5,
#' 3)}.  To set the same matrix with names of \code{iivKa} and \code{iivCL}, use
#' \code{iivKa + iivCL~c(2, 2.5, 3)}.  To set a correlation matrix with standard
#' deviations on the diagonal, use \code{cor()} like \code{iivKa + iivCL~cor(2,
#' -0.5, 3)}.
#'
#' Values may be fixed (and therefore not estimated) using either the name
#' \code{fixed} at the end of the assignment or by calling \code{fixed()} as a
#' function for the value to fix.  For 'theta' and 'sigma', either the estimate
#' or the full definition (including lower and upper bounds) may be included in
#' the fixed setting.  For example, the following are all effectively equivalent
#' to set a 'theta' or 'sigma' to a fixed value (because the lower and upper
#' bounds are ignored for a fixed value): \code{tvCL <- fixed(1)}, \code{tvCL <-
#' fixed(0, 1)}, \code{tvCL <- fixed(0, 1, 2)}, \code{tvCL <- c(0, fixed(1),
#' 2)}, or \code{tvCL <- c(0, 1, fixed)}.  For 'omega' assignment, the full
#' block or none of the block must be set as \code{fixed}.  Examples of setting
#' an 'omega' value as fixed are: \code{iivKa~fixed(1)}, \code{iivKa +
#' iivCL~fixed(1, 2, 3)}, or \code{iivKa + iivCL~c(1, 2, 3, fixed)}.  Anywhere
#' that \code{fixed} is used, \code{FIX}, \code{FIXED}, or \code{fix} may be
#' used equivalently.
#'
#' For any value, standard mathematical operators or functions may be used to
#' define the value.  For example, \code{exp(2)} and \code{24*30} may be used to
#' define a value anywhere that a number can be used (e.g. lower bound,
#' estimate, upper bound, variance, etc.).
#'
#' Values may be labeled using the \code{label()} function after the assignment.
#' Labels are are used to make reporting easier by giving a human-readable
#' description of the parameter, but the labels do not have any effect on
#' estimation.  The typical way to set a label so that the parameter \code{tvCL}
#' has a label of "Typical Value of Clearance (L/hr)" is \code{tvCL <- 1;
#' label("Typical Value of Clearance (L/hr)")}.
#'
#' \code{nlmixr} will attempt to determine some back-transformations for the
#' user.  For example, \code{CL <- exp(tvCL)} will detect that \code{tvCL} must
#' be back-transformed by \code{exp()} for easier interpretation.  When you want
#' to control the back-transformation, you can specify the back-transformation
#' using \code{backTransform()} after the assignment.  For example, to set the
#' back-transformation to \code{exp()}, you can use \code{tvCL <- 1;
#' backTransform(exp())}.
#'
#' @param ini Ini block or nlmixr model object
#' @param ... Other arguments parsed by nlmixr
#' @return bounds expression or parsed ui object
#' @author Matthew L. Fidler
#' @export
ini <- function(ini, ...) {
  if (is(substitute(ini), "{")) {
    .f <- eval(parse(text = sprintf("function() %s", paste(.deparse(substitute(ini)), collapse = "\n"))))
    .fb <- nlmixrBounds(.f)
    if (!exists(".ini", parent.frame(1))) {
      assign(".ini", .f, parent.frame(1))
    } else {
      .ini <- get(".ini", parent.frame(1))
      if (is(.ini, "function")) {
        assign(".ini", .f, parent.frame(1))
      }
    }
    return(.fb)
  } else {
    .uif <- .getUif(ini)
    .ini <- .uif$ini
    .call <- match.call(expand.dots = TRUE)[-1]
    .call <- .call[-1]
    .ns <- names(.call)
    if (!is.null(.ns)) {
      .ini <- .uif$ini
      for (.n in .ns) {
        .w <- which(.n == .ini$name)
        if (length(.w) == 1) {
          if (any(.deparse(.call[[.n]]) == c("fix", "fixed", "FIX", "FIXED"))) {
            if (.uif$ini$fix[.w]) {
              warning(sprintf("Trying to fix '%s', but already fixed.", .n))
            } else {
              .uif$ini$fix[.w] <- TRUE
            }
          } else if (any(.deparse(.call[[.n]]) == c("unfix", "unfixed", "UNFIX", "UNFIXED"))) {
            if (.uif$ini$fix[.w]) {
              .uif$ini$fix[.w] <- FALSE
            } else {
              warning(sprintf("Trying to unfix '%s', but not fixed.", .n))
            }
          } else if (regexpr(rex::rex(or(c("fix", "fixed", "FIX", "FIXED")), "(", anything, ")"), .deparse(.call[[.n]])) != -1) {
            .val <- eval(.call[[.n]][[2]])
            if (.uif$ini$fix[.w]) {
              warning(sprintf("Trying to fix '%s', but already fixed.  Still assigned to '%s'", .n, .val))
            } else {
              .uif$ini$fix[.w] <- TRUE
            }
            .uif$ini$est[.w] <- .val
          } else if (regexpr(rex::rex(or(c("unfix", "unfixed", "UNFIX", "UNFIXED")), "(", anything, ")"), .deparse(.call[[.n]])) != -1) {
            .val <- eval(.call[[.n]][[2]])
            if (.uif$ini$fix[.w]) {
              .uif$ini$fix[.w] <- FALSE
            } else {
              warning(sprintf("Trying to unfix '%s', but not fixed.  Still assigned to '%s'", .n, .val))
            }
            .uif$ini$est[.w] <- .val
          } else {
            .val <- try(eval(.call[[.n]]), silent = TRUE)
            .found <- TRUE
            if (inherits(.val, "try-error")) {
              .found <- FALSE
              .val <- as.character(.call[[.n]])
              .val2 <- .val
              .frames <- seq(1, sys.nframe())
              .frames <- .frames[.frames != 0]
              for (.f in .frames) {
                .env <- parent.frame(.f)
                if (exists(.val, envir = .env)) {
                  .val2 <- try(get(.val, envir = .env), silent = TRUE)
                  if (!inherits(.val2, "try-error")) {
                    .val <- .val2
                    .found <- TRUE
                    break
                  }
                }
              }
            }
            if (!.found) {
              stop(sprintf("object '%s' not found", .val))
            }
            if (length(.val) == 1) {
              .uif$ini$est[.w] <- .val
            } else if (length(.val) == 2) {
              .uif$ini$lower[.w] <- .val[1]
              .uif$ini$est[.w] <- .val[2]
              ## Warning here? The upper should be inf?
              .uif$ini$est[.w] <- Inf
            } else if (length(.val) == 3) {
              .uif$ini$lower[.w] <- .val[1]
              .uif$ini$est[.w] <- .val[2]
              ## Warning here? The upper should be inf?
              .uif$ini$upper[.w] <- .val[3]
            } else {
              stop(sprintf("Cannot figure out what you are trying to do to the '%s' estimate.", .n))
            }
          }
        } else {
          warning(sprintf("The model does not have a parameter named '%s', modification ignored.", .n))
        }
      }
    } else if (length(.call) == 1) {
      .lst <- eval(.call[[1]])
      .ns <- names(.lst)
      if (inherits(.lst, "list") && !is.null(.ns)) {
        for (.n in .ns) {
          .w <- which(.n == .ini$name)
          if (length(.w) == 1) {
            .val <- .lst[[.n]]
            if (length(.val) == 1) {
              .uif$ini$est[.w] <- .val
            } else if (length(.val) == 2) {
              .uif$ini$lower[.w] <- .val[1]
              .uif$ini$est[.w] <- .val[2]
              .uif$ini$upper[.w] <- Inf
            } else if (length(.val) == 3) {
              .uif$ini$lower[.w] <- .val[1]
              .uif$ini$est[.w] <- .val[2]
              .uif$ini$upper[.w] <- .val[3]
            } else {
              stop(sprintf("Cannot figure out what to do with '%s'", .n))
            }
          }
        }
      } else if (inherits(.lst, "numeric") && !is.null(.ns)) {
        for (.n in .ns) {
          .w <- which(.n == .ini$name)
          if (length(.w) == 1) {
            .uif$ini$est[.w] <- .lst[.n]
          }
        }
      } else {
        ## If this is a+b~c(...) parse using lotri
        ## FIXME handle conditional changes like a+b~c(...) | occ
        .code <- paste(deparse(.lst), collapse = " ")
        .uif <- .processLotri(.code, .uif)
      }
    } else {
      stop("Do not know what to do with the model's ini call.")
    }
    return(.uif)
  }
}


.thetaModelReg <- rex::rex(or(
  group(start, or("tv", "t", "pop", "POP", "Pop", "TV", "T", "cov", "err", "eff")),
  group(or("tv", "pop", "POP", "pop", "TV", "cov", "err", "eff"), end)
))
.etaParts <- c(
  "eta", "ETA", "Eta", "ppv", "PPV", "Ppv", "iiv", "Iiv", "bsv", "Bsv", "BSV",
  "bpv", "Bpv", "BPV", "psv", "PSV", "Psv"
)
.etaModelReg <- rex::rex(or(group(start, or(.etaParts)), group(or(.etaParts), end)))

##' nlmixr model block
##'
##' @param model Model specification
##' @param ... Other arguments to model object parsed by nlmixr
##' @param .lines This is an internal argument when code{model} is
##'     being called recursively and should not be used.
##' @return Parsed UI object
##' @author Matthew L. Fidler
##' @export
model <- function(model, ..., .lines = NULL) {
  if (is(substitute(model), "{")) {
    .f <- eval(parse(text = sprintf("function() %s", paste(.deparse(substitute(model)), collapse = "\n"))))
    if (!exists(".model", parent.frame(1))) {
      assign(".model", .f, parent.frame(1))
    } else {
      .model <- get(".model", parent.frame(1))
      if (is(.model, "function")) {
        assign(".model", .f, parent.frame(1))
      }
    }
    if (exists(".ini", parent.frame(1))) {
      .ini <- get(".model", parent.frame(1))
      .model <- get(".model", parent.frame(1))
      if (is(.ini, "function") & is(.model, "function")) {
        return(nlmixrUI(parent.frame(1)))
      }
    }
    return(.f)
  } else {
    .uif <- .getUif(model)
    .ini <- .as.data.frame(.uif$ini)
    .call <- match.call(expand.dots = TRUE)[-(1:2)]
    if (is.null(.lines)) {
      .lines <- .deparse(.call[[1]])
    }
    .f <- try(eval(parse(text = sprintf("function() %s", paste(.lines, collapse = "\n")))), silent = TRUE)
    .origLines <- strsplit(.uif$fun.txt, "\n")[[1]]
    .fNew <- c("function(){", .origLines, "}")
    .origLines <- .fNew
    if (!inherits(.f, "try-error")) {
      .fOrig <- eval(parse(text = paste(.fNew, collapse = "\n")))
      .lhsOrig <- nlmixrfindLhs(body(.fOrig))
      .lhsNew <- nlmixrfindLhs(body(.f))
      .shared <- intersect(.lhsNew, .lhsOrig)
      .distOrig <- .findDist(body(.fOrig))
      .distNew <- .findDist(body(.f))
      .sharedDist <- intersect(.distNew, .distOrig)
      if (length(.shared) == 0 && length(.sharedDist) == 0) {
        .code <- paste(.lines, collapse = " ")
        .processLotri(.code, .uif)
        .uif <- try(.processLotri(.code, .uif), silent = TRUE)
        if (inherits(.uif, "try-error")) {
          stop("Could not find a part of the model to modify.")
        }
        return(.uif)
      }
      if (length(.sharedDist) > 0) {
        for (.v in .sharedDist) {
          .regShared <- rex::rex(start, any_spaces, .v, any_spaces, "~")
          .newLine <- .lines[regexpr(.regShared, .lines) != -1]
          .oldLine <- .origLines[regexpr(.regShared, .origLines) != -1]
          .w <- which(regexpr(.regShared, .fNew) != -1)
          if (length(.w) > 1) stop("Cannot modify a multi-line definition")
          .fNew[.w] <- .newLine
          .oldVarsL <- allVars(body(eval(parse(text = sprintf("function(){%s}", .oldLine)))))
          .newVarsL <- allVars(body(eval(parse(text = sprintf("function(){%s}", .newLine)))))
          .rmVars <- setdiff(.oldVarsL, .newVarsL)
          .addVars <- setdiff(.newVarsL, .oldVarsL)
          for (.rm in .rmVars) {
            .wh <- .ini[.ini$name == .rm, ]
            if (length(.wh$name) == 1) {
              if (is.na(.wh$ntheta)) {
                ## removing eta
                .curEta <- .wh$neta1
                .maxEta <- max(.ini$neta1, na.rm = TRUE)
                .s <- seq(.curEta, .maxEta)
                .s <- .s[-length(.s)]
                .w <- unique(sort(c(
                  which(is.na(.ini$neta1)),
                  which(.ini$neta1 != .curEta),
                  which(.ini$neta2 != .curEta)
                )))
                .ini <- .ini[.w, ]
                for (.rmI in .s) {
                  .ini$neta1[.ini$neta1 == .rmI + 1] <- .rmI
                  .ini$neta2[.ini$neta2 == .rmI + 1] <- .rmI
                }
              } else {
                .curTheta <- .wh$ntheta
                .maxTheta <- max(.ini$ntheta, na.rm = TRUE)
                .s <- seq(.curTheta, .maxTheta)
                .s <- .s[-length(.s)]
                .w <- unique(sort(c(
                  which(is.na(.ini$ntheta)),
                  which(.ini$ntheta != .curTheta)
                )))
                .ini <- .ini[.w, ]
                for (.rmI in .s) {
                  .ini$ntheta[.ini$ntheta == .rmI + 1] <- .rmI
                }
              }
            }
          }
          ## Now add variables; Currently assume thetas; Perhaps something better?
          ## Something other than thetas are not currently supported by UI currently.
          for (.new in .addVars) {
            if (!any(.ini$name == .new)) {
              .maxTheta <- max(.ini$ntheta, na.rm = TRUE)
              .d2 <- nlmixrBoundsTemplate
              .d2$name <- .new
              .d2$ntheta <- .maxTheta + 1
              .d2$est <- 1
              .d2$lower <- -Inf
              .d2$upper <- Inf
              .ini <- rbind(.ini, .d2)
            }
          }
        }
      }
      if (length(.shared) > 0) {
        for (.v in .shared) {
          .regShared <- rex::rex(start, any_spaces, .v, any_spaces, or("=", "<-"))
          .newLine <- .lines[regexpr(.regShared, .lines) != -1]
          .oldLine <- .origLines[regexpr(.regShared, .origLines) != -1]
          .w <- which(regexpr(.regShared, .fNew) != -1)
          if (length(.w) > 1) stop("Cannot modify a multi-line definition")
          .fNew[.w] <- .newLine
          .oldVarsL <- allVars(body(eval(parse(text = sprintf("function(){%s}", .oldLine)))))
          .newVarsL <- allVars(body(eval(parse(text = sprintf("function(){%s}", .newLine)))))
          .rmVars <- setdiff(.oldVarsL, .newVarsL)
          .addVars <- setdiff(.newVarsL, .oldVarsL)
          for (.rm in .rmVars) {
            .wh <- .ini[.ini$name == .rm, ]
            if (length(.wh$name) == 1) {
              if (is.na(.wh$ntheta)) {
                ## removing eta
                .curEta <- .wh$neta1
                .maxEta <- max(.ini$neta1, na.rm = TRUE)
                .s <- seq(.curEta, .maxEta)
                .s <- .s[-length(.s)]
                .w <- unique(sort(c(
                  which(is.na(.ini$neta1)),
                  which(.ini$neta1 != .curEta),
                  which(.ini$neta2 != .curEta)
                )))
                .ini <- .ini[.w, ]
                for (.rmI in .s) {
                  .ini$neta1[.ini$neta1 == .rmI + 1] <- .rmI
                  .ini$neta2[.ini$neta2 == .rmI + 1] <- .rmI
                }
              } else {
                .curTheta <- .wh$ntheta
                .maxTheta <- max(.ini$ntheta, na.rm = TRUE)
                .s <- seq(.curTheta, .maxTheta)
                .s <- .s[-length(.s)]
                .w <- unique(sort(c(
                  which(is.na(.ini$ntheta)),
                  which(.ini$ntheta != .curTheta)
                )))
                .ini <- .ini[.w, ]
                for (.rmI in .s) {
                  .ini$ntheta[.ini$ntheta == .rmI + 1] <- .rmI
                }
              }
            }
          }
          ## Now add variables; Currently parsed based on variable name; Perhaps something better?
          for (.new in .addVars) {
            if (!any(.ini$name == .new)) {
              if (regexpr(.etaModelReg, .new) != -1) {
                .maxEta <- suppressWarnings(max(.ini$neta1, na.rm = TRUE))
                if (is.infinite(.maxEta)) .maxEta <- 0
                .maxTheta <- max(.ini$ntheta, na.rm = TRUE)
                .d2 <- nlmixrBoundsTemplate
                .d2$neta1 <- .maxEta + 1
                .d2$neta2 <- .maxEta + 1
                .d2$est <- 1
                .d2$lower <- -Inf
                .d2$upper <- Inf
                .d2$condition <- "ID"
                .ini <- rbind(.ini, .d2)
              } else if (regexpr(.thetaModelReg, .new) != -1) {
                .maxTheta <- max(.ini$ntheta, na.rm = TRUE)
                .d2 <- nlmixrBoundsTemplate
                .d2$ntheta <- .maxTheta + 1
                .d2$est <- 1
                .d2$lower <- -Inf
                .d2$upper <- Inf
                .d2$name <- .new
                .ini <- rbind(.ini, .d2)
              }
            }
          }
        }
      }

      class(.ini) <- c("nlmixrBounds", "data.frame")
      .model <- eval(parse(text = paste(.fNew, collapse = "\n"), keep.source = TRUE))
      return(.finalizeUiModel(
        nlmixrUIModel(.model, .ini, NULL),
        as.list(.uif$meta)
      ))
    }
  }
}

.nlmixrUpdate <- function(object, ...) {
  .uif <- .getUif(object)
  .iniNames <- paste(.uif$ini$name)
  .call <- match.call(expand.dots = TRUE)[-c(1:2)]
  .callNames <- names(.call)
  if (is.null(.callNames)) {
    .callNames <- rep("", length(.call))
  }
  if (any(.callNames %in% c("data", "est", "control", "table", "save"))) {
    stop("Cannot update these arguments: 'data', 'est', 'control', 'table', or 'save'\nUse these options with nlmixr.")
  }
  ## Any time the parameters are named, use ini call
  .iniArgs <- which(.callNames %in% .iniNames)
  .iniNames2 <- .callNames[.iniArgs]
  if (length(.iniArgs) > 0) {
    .uif <- do.call(ini, c(list(.uif), setNames(as.list(.call[.iniArgs]), .iniNames2)))
    .iniI <- seq_along(.call)[-.iniArgs]
  } else {
    .iniI <- seq_along(.call)
  }
  for (.i in .iniI) {
    ## Call model() multiple times if needed.
    if (.callNames[.i] != "") {
      ## construct model(.uif,{name()<-expr()})
      .cn <- eval(parse(text = paste0("quote(", .callNames[.i], ")")))
      .uif <- eval(as.call(list(
        `model`, quote(.uif),
        as.call(list(
          quote(`{`),
          as.call(list(
            quote(`<-`), .cn,
            .call[[.i]]
          ))
        ))
      )))
    } else {
      .x <- .call[[.i]]
      .didCall <- FALSE
      if (is.call(.x) || is.pairlist(.x)) {
        if (identical(.x[[1]], quote(`c`)) ||
          identical(.x[[1]], quote(`list`))) {
          .uif <- eval(as.call(list(`ini`, quote(.uif), .x)))
          .didCall <- TRUE
        }
      } else if (is.name(.x)) {
        .cur <- 1
        .tmp <- eval(.x, sys.parent(1))
        if (is(.tmp, "list") || is(.tmp, "numeric") || is(.tmp, "integer")) {
          .uif <- eval(bquote(ini(.(.uif), .(.tmp))))
        } else {
          .uif <- model(.uif, .lines = .deparse(.tmp))
        }
        .didCall <- TRUE
      }
      if (!.didCall) {
        .uif <- model(.uif, .lines = .deparse(.x))
      }
    }
  }
  return(.uif)
}

##' @export
update.nlmixrUI <- .nlmixrUpdate

##' @export
update.nlmixrFitCore <- .nlmixrUpdate

##' @export
update.function <- .nlmixrUpdate

.finalizeUiModel <- function(fun2, meta) {
  class(fun2) <- "nlmixrUI"
  var.ini <- c(fun2$focei.names, fun2$eta.names)
  var.def <- fun2$all.vars
  diff <- setdiff(var.ini, var.def)
  if (length(diff) > 0) {
    stop(sprintf("Model error: initial estimates provided without variables being used: %s", paste(diff, collapse = ", ")))
  }
  ns <- fun2$name[which(!is.na(fun2$neta1) & !is.na(fun2$err))]
  if (length(ns) > 0) {
    stop(sprintf("Residual error component(s) need to be defined with assignment ('=' or '<-') in ini block (not '~'): %s", paste(ns, collapse = ", ")))
  }
  ns <- fun2$name[is.na(fun2$est)]
  if (length(ns) > 0) {
    stop(sprintf("The following parameters initial estimates are NA: %s", paste(ns, collapse = ", ")))
  }
  fun2$meta <- list2env(meta, parent = emptyenv())
  .mv <- RxODE::rxModelVars(fun2$rxode)
  .tmp <- paste(fun2$nmodel$predDf$var)
  .testVars <- c(.mv$lhs, .mv$state)
  .tmp <- .tmp[!(.tmp %in% .testVars)]
  if (length(.tmp > 0)) {
    .predDf <- fun2$predDf
    if (length(.predDf$var) > 1) {
      .predDf <- .predDf[.predDf$var %in% .tmp, ]
      .predDf <- .predDf[.predDf$var == .predDf$cond, , drop = FALSE]
      if (length(.predDf$var) > 0) {
        stop(sprintf("Multiple compartment models with expressions need to be conditioned by `|`\n ie log(cp) ~ add(err) | cmt\n The following endpoints need to be corrected: %s", paste(.predDf$var, collapse = ", ")))
      }
    }
    .tmp <- lapply(.tmp, function(x) {
      .newPars <- RxODE::rxModelVars(paste0("nlmixr_tmp=", x))$params
      if (all(.newPars %in% .testVars)) {
        return(NA_character_)
      }
      return(.newPars)
    })
    .tmp <- stats::na.omit(unlist(.tmp))
    if (length(.tmp) > 0) {
      stop(sprintf(
        "Modeled responses need to be defined in the model; Add definition for: %s",
        paste(.tmp, collapse = ", ")
      ))
    }
  }
  return(fun2)
}

##' Prepares the UI function and returns a list.
##'
##' @param fun UI function
##' @return nlmixr UI function
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nlmixrUI <- function(fun) {
  if (is(fun, "function")) {
    lhs0 <- nlmixrfindLhs(body(fun))
    dum.fun <- function() {
      return(TRUE)
    }
    env.here <- environment(dum.fun)
    env <- new.env(parent = .GlobalEnv)
    assign("fun", fun, env)
    fun2 <- attr(fun, "srcref")
    if (is.null(fun2)) {
      message(sprintf("Cannot run this way in non-interactive mode.\nTry running:\n\nR -e 'source(\"script\", keep.source=TRUE)'\n\nFor Batch-mode type of running, you can use:\n\nR -e \"source('script.R', keep.source=TRUE, echo=TRUE)\" %s> script.Rout 2>&1", ifelse((.Platform$OS.type == "unix"), "", "1")))
      stop("option \"keep.source\" must be TRUE for nlmixr models.")
      return(eval(fun(), parent.frame(1)))
    }
    fun2 <- as.character(fun2, useSource = TRUE)
    rg <- rex::rex("function", any_spaces, "(", anything, ")")
    w <- which(regexpr(rg, fun2) != -1)
    if (length(w) > 0) {
      w <- w[1]
      fun2[w] <- sub(rg, "", fun2[w])
    }
    fun2 <- gsub(rex::rex(boundary, "ini(", any_spaces, "{"), ".ini <- function()({", fun2)
    fun2 <- gsub(rex::rex(boundary, "model(", any_spaces, "{"), ".model <- function()({", fun2)
    if (fun2[length(fun2)] != "}") {
      fun2[length(fun2)] <- sub(rex::rex("}", end), "", fun2[length(fun2)])
      fun2[length(fun2) + 1] <- "}"
    }
    w <- which(regexpr(rex::rex(start, any_spaces, "#", anything), fun2) != -1)
    if (length(w) > 0 && all(lhs0 != "desc")) {
      w2 <- w[1]
      if (length(w) > 1) {
        for (i in 2:length(w)) {
          if (w[i] - 1 == w[i - 1]) {
            w2[i] <- w[i]
          } else {
            break
          }
        }
      }
      desc <- paste(gsub(
        rex::rex(any_spaces, end), "",
        gsub(rex::rex(start, any_spaces, any_of("#"), any_spaces), "", fun2[w2])
      ),
      collapse = " "
      )
      lhs0 <- c(lhs0, "desc")
    }
    .model <- .ini <- NULL ## Rcheck hack
    eval(parse(text = paste(fun2, collapse = "\n"), keep.source = TRUE))
  } else if (is(fun, "environment")) {
    env.here <- fun
    .model <- fun$.model
    .ini <- fun$.ini
    env <- new.env(parent = .GlobalEnv)
    lhs0 <- c("data", "desc", "ref", "imp", "est", "control", "table")
    warning("Some information (like parameter labels) is lost by evaluating a nlmixr function")
    fun <- paste0(
      "function(){\n",
      paste(sapply(lhs0, function(var) {
        if (exists(var, envir = env.here)) {
          if (!inherits(get(var, envir = env.here), "function")) {
            return(sprintf("\n%s <- %s;", var, .deparse(get(var, envir = env.here))))
          } else {
            return("")
          }
        } else {
          return("")
        }
      }), collapse = ""),
      sprintf("\nini(%s)", paste(.deparse(body(.ini)), collapse = "\n")),
      sprintf("\nmodel(%s)", paste(.deparse(body(.model)), collapse = "\n")),
      "\n}"
    )
    fun <- eval(parse(text = fun))
    assign("fun", fun, env)
  }
  ini <- nlmixrBounds(.ini)
  meta <- list()
  for (var in lhs0) {
    if (!any(var == c(".ini", ".model")) &&
      exists(var, envir = env.here)) {
      meta[[var]] <- get(var, envir = env.here)
    }
  }
  if (inherits(.ini, "try-error")) {
    stop("Error parsing initial estimates.")
  }

  fun2 <- nlmixrUIModel(.model, ini, fun)
  ## if (inherits(fun2, "try-error")){
  ##     stop("Error parsing model.")
  ## }
  return(.finalizeUiModel(fun2, meta))
}
nlmixrUI.multipleEndpoint <- function(x) {
  if (length(x$nmodel$predDf$cond) > 1) {
    .info <- x$nmodel$predDf
    if (getOption("RxODE.combine.dvid", TRUE)) {
      .info <- .info[order(.info$dvid), ]
    }
    .info <- with(.info, .data.frame(
      variable = paste(var, "~", ifelse(use.utf(), "\u2026", "...")),
      cmt = paste0("cmt='", cond, "' or cmt=", cmt),
      "dvid*" = ifelse(is.na(dvid), "",
        paste0("dvid='", cond, "' or dvid=", dvid)
      ),
      check.names = FALSE
    ))
    if (!getOption("RxODE.combine.dvid", TRUE)) {
      .info <- .info[, names(.info) != "dvid*"]
    }
    if (requireNamespace("huxtable", quietly = TRUE)) {
      .hux <- huxtable::hux(.info) %>%
        huxtable::add_colnames() %>%
        huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
        huxtable::set_position("center") %>%
        huxtable::set_all_borders(TRUE)
      if (getOption("RxODE.combine.dvid", TRUE)) {
        .hux <- .hux %>%
          huxtable::add_footnote("* If dvids are outside this range, all dvids are re-numered sequentially, ie 1,7, 10 becomes 1,2,3 etc")
      }
    } else {
      .hux <- .info
    }
    return(.hux)
  } else {
    return(NULL)
  }
}
##' Print UI function
##'
##' @param x  UI function
##' @param ... other arguments
##' @author Matthew L. Fidler
##' @export
print.nlmixrUI <- function(x, ...) {
  cat(cli::rule(x$model.desc, line = "bar2"), "\n")
  cat(cli::rule(crayon::bold("Initialization:")), "\n")
  print(x$ini)
  if (length(x$all.covs) > 0) {
    cat("\n Covariates or Uninitialized Parameters ($all.covs)\n")
    print(x$all.covs)
  }
  if (length(x$predDf$cond) > 1) {
    cat(cli::rule(paste0(crayon::bold("Multiple Endpoint Model"), " (", crayon::bold$blue("$multipleEndpoint"), "):")), "\n")

    if (requireNamespace("huxtable", quietly = TRUE)) {
      x$multipleEndpoint %>%
        huxtable::print_screen(colnames = FALSE)
    } else {
      print(x$multipleEndpoint)
    }
    cat("\n")
  }
  .mu <- x$muRefTable
  if (length(.mu) > 0) {
    cat(cli::rule(paste0(
      crayon::bold(paste0(ifelse(use.utf(), "\u03bc", "mu"), "-referencing")),
      " (", crayon::bold$blue("$muRefTable"), "):"
    )), "\n")
    if (requireNamespace("huxtable", quietly = TRUE)) {
      .mu %>%
        huxtable::print_screen(colnames = FALSE)
    } else {
      print(.mu)
    }
    cat("\n")
  }
  cat(cli::rule(crayon::bold(sprintf("Model%s:", ifelse(class(x$rxode) == "RxODE", " (RxODE)", "")))), "\n")
  cat(x$fun.txt, "\n")
  cat(cli::rule(line = "bar2"), "\n")
}

## This is a list of supported distributions with the number of arguments they currently support.
dists <- list(
  "dpois" = 0,
  "dbinom" = 0:1,
  "dbern" = 0,
  "bern" = 0,
  "dbeta" = 2:3,
  ##
  ## "dnbinom"=2:3,  ## dnbinom is in R; FIXME: how does ot compare to dneg_binomial
  ## "dneg_binomial", ## not in base R (but in glnmm2)
  ##
  ## Available as external package http://ugrad.stat.ubc.ca/R/library/rmutil/html/BetaBinom.html
  ## "dbetabinomial", ## not in base R (but in glnmm2)
  "dt" = 1:2,
  "pois" = 0,
  "binom" = 0:1,
  "beta" = 2:3,
  "t" = 1:2,
  "add" = 1,
  "norm" = 1,
  "dnorm" = 1,
  "prop" = 1,
  "pow" = 2,
  "tbs" = 1,
  "boxCox" = 1,
  "tbsYj" = 1,
  "yeoJohnson" = 1,
  "logn" = 1,
  "dlogn" = 1,
  "lnorm" = 1,
  "dlnorm" = 1
)

distsPositive <- c("add", "norm", "dnorm", "prop", "pow", "logn", "dlogn", "lnorm", "dlnorm")

allVars <- function(x) {
  defined <- character()
  this.env <- environment()
  f <- function(x) {
    if (is.atomic(x)) {
      character()
    } else if (is.name(x)) {
      as.character(x)
    } else if (is.call(x) || is.pairlist(x)) {
      if (identical(x[[1]], quote(`~`)) ||
        identical(x[[1]], quote(`=`)) ||
        identical(x[[1]], quote(`<-`))) {
        if (is.call(x[[3]])) {
          ret <- unique(unlist(lapply(x[[3]][-1], f)))
        } else {
          ret <- unique(unlist(lapply(x[[3]], f)))
        }
        ret <- ret[!(ret %in% defined)]
        assign("defined", unique(c(defined, x[[2]])), this.env)
        return(ret)
      } else {
        children <- lapply(x[-1], f)
        unique(unlist(children))
      }
    } else {
      stop("Don't know how to handle type ", typeof(x),
        call. = FALSE
      )
    }
  }
  f(x)
}

allNames <- function(x) {
  if (is.atomic(x)) {
    character()
  } else if (is.name(x)) {
    as.character(x)
  } else if (is.call(x) || is.pairlist(x)) {
    children <- lapply(x[-1], allNames)
    unique(unlist(children))
  } else {
    stop("Don't know how to handle type ", typeof(x),
      call. = FALSE
    )
  }
}

allCalls <- function(x) {
  if (is.atomic(x) || is.name(x)) {
    character()
  } else if (is.call(x)) {
    fname <- as.character(x[[1]])
    children <- lapply(x[-1], allCalls)
    unique(c(fname, unlist(children)))
  } else if (is.pairlist(x)) {
    unique(unlist(lapply(x[-1], allCalls), use.names = FALSE))
  } else {
    stop("Don't know how to handle type ", typeof(x), call. = FALSE)
  }
}

##' Find the assignments in R expression
##'
##' @param x R expression
##' @return list of assigned parameters
##' @author Hadley Wickham and Matthew L. Fidler
##' @keywords internal
##' @export
nlmixrfindLhs <- function(x) {
  ## Modified from http://adv-r.had.co.nz/Expressions.html find_assign4
  if (is.atomic(x) || is.name(x)) {
    character()
  } else if (is.call(x)) {
    if ((identical(x[[1]], quote(`<-`)) ||
      identical(x[[1]], quote(`=`))) &&
      is.name(x[[2]])) {
      lhs <- as.character(x[[2]])
    } else {
      lhs <- character()
    }
    unique(c(lhs, unlist(lapply(x, nlmixrfindLhs))))
  } else if (is.pairlist(x)) {
    unique(unlist(lapply(x, nlmixrfindLhs)))
  } else {
    stop("Don't know how to handle type ", typeof(x),
      call. = FALSE
    )
  }
}

.findDist <- function(x) {
  ## Modified from http://adv-r.had.co.nz/Expressions.html find_assign4
  if (is.atomic(x) || is.name(x)) {
    character()
  } else if (is.call(x)) {
    if (identical(x[[1]], quote(`~`))) {
      if (is.name(x[[2]])) {
        lhs <- as.character(x[[2]])
      } else {
        if (is.call(x[[2]]) && identical(x[[2]][[1]], quote(`linCmt`))) {
          lhs <- "linCmt()"
        } else {
          lhs <- character()
        }
      }
    } else {
      lhs <- character()
    }
    unique(c(lhs, unlist(lapply(x, .findDist))))
  } else if (is.pairlist(x)) {
    unique(unlist(lapply(x, .findDist)))
  } else {
    stop("Don't know how to handle type ", typeof(x),
      call. = FALSE
    )
  }
}

## cauchy = t w/df=1; Support?
## dgeom is a special negative binomial; Support?
## fixme
unsupported.dists <- c(
  "dchisq", "chisq", "dexp", "df", "f", "dgeom", "geom",
  "dhyper", "hyper", "dunif", "unif",
  "dweibull", "weibull",
  ## for testing...
  "nlmixrDist"
)
add.dists <- c("add", "prop", "norm", "pow", "dnorm", "logn", "lnorm", "dlnorm", "tbs", "tbsYj", "boxCox", "yeoJohnson")
nlmixrUIModel <- function(fun, ini = NULL, bigmodel = NULL) {
  .md5 <- digest::digest(list(
    .deparse1(body(fun)), .deparse(ini), .deparse(bigmodel),
    ## Should give different models for different nlmixr versions
    sessionInfo()$otherPkgs$nlmixr$Version
  ))
  .wd <- RxODE::rxTempDir()
  if (.wd == "") {
    warning("rxTempDir did not work.")
    .wd <- tempfile()
    dir.create(.wd, recursive = TRUE)
  } else {
    ## This makes all the parsing files, cpp and so in their own
    ## directory.  No collisions.
    suppressWarnings({
      dir.create(.wd, recursive = TRUE)
    })
  }
  .uiFile <- file.path(.wd, paste0("ui-", .md5))
  .uiLock <- paste0(.uiFile, ".lock")
  .uiBad <- paste0(.uiFile, ".bad")
  .uiFile <- paste0(.uiFile, ".uid")
  if (file.exists(.uiLock)) {
    message(sprintf("Waiting for UI to parse on another thread (%s)", .uiLock), appendLF = FALSE)
    while (file.exists(.uiLock)) {
      Sys.sleep(0.5)
      message(".", appendLF = FALSE)
    }
    message("")
    if (file.exists(.uiBad)) {
      load(.uiBad)
      for (.w in .lst$ws) {
        warning(.w)
      }
      stop(paste(.lst$em, "(another thread)"))
    }
    load(file = .uiFile)
    return(ret)
  } else if (file.exists(.uiBad)) {
    load(.uiBad)
    for (.w in .lst$ws) {
      warning(.w)
    }
    stop(sprintf("%s\nBad parsed model (cached %s).", .lst$em, .uiBad))
  } else if (file.exists(.uiFile)) {
    load(file = .uiFile)
    return(ret)
  } else {
    sink(.uiLock)
    cat("")
    sink()
    on.exit(
      {
        unlink(.uiLock)
      },
      add = TRUE
    )
    sink(.uiBad)
    cat("")
    sink()
    .thisEnv <- environment()
    .thisEnv$ws <- character(0)
    ret <- tryCatch(suppressWarnings(withCallingHandlers(.nlmixrUIModel(fun, ini, bigmodel),
      warning = function(w) {
        assign("ws", unique(c(w$message, .thisEnv$ws)), .thisEnv)
      }
    )), error = function(e) {
      assign("em", e$message, .thisEnv)
    })
    .lst <- list()
    if (exists("ws", envir = .thisEnv)) {
      .lst$ws <- .thisEnv$ws
    }
    if (exists("em", envir = .thisEnv)) {
      .lst$em <- .thisEnv$em
      save(.lst, file = .uiBad)
      ## stop(.lst$em);
      ## Let it error out on its own for the backtrace...
      ret <- .nlmixrUIModel(fun, ini, bigmodel)
    }
    save(ret, .lst, file = .uiFile)
    unlink(.uiBad)
    return(ret)
  }
}
.nlmixrUIModel <- function(fun, ini = NULL, bigmodel = NULL) {
  ## Parses the UI function to extract predictions and errors, and the other model specification.
  .fun000 <- fun
  rxode <- FALSE
  all.names <- allNames(body(fun))
  .diff <- setdiff(paste(ini$name), all.names)
  if (length(.diff) > 0) {
    .diff <- .diff[regexpr("[(]", .diff) == -1]
    if (length(.diff) > 0) {
      stop(sprintf("The following parameter(s) were in the ini block but not in the model block: %s", paste(.diff, collapse = ", ")))
    }
  }
  if (any(regexpr(rex::rex(start, or("rx_", "nlmixr_")), paste(all.names)) != -1)) {
    stop("Parameters/States/Variables cannot start with `rx_` or `nlmixr_`")
  }
  all.vars <- allVars(body(fun))
  all.funs <- allCalls(body(fun))
  all.lhs <- nlmixrfindLhs(body(fun))
  errs.specified <- c()
  add.prop.errs <- .data.frame(y = character(), add = logical(), prop = logical())
  bounds <- ini
  theta.names <- c()
  theta.ord <- c()
  eta.names <- c()
  .mu.ref <- list()
  .oneTheta <- c()
  cov.ref <- list()
  cov.theta <- c()
  log.theta <- c()
  log.eta <- c()
  this.env <- environment()
  if (!is.null(ini)) {
    unnamed.thetas <- ini$ntheta[(!is.na(ini$ntheta) & is.na(ini$name))]
    if (length(unnamed.thetas) > 0) {
      stop(sprintf("The following THETAs are unnamed: %s", paste(sprintf("THETA[%d]", unnamed.thetas), collapse = ", ")))
    }
    unnamed.etas <- ini$neta1[!is.na(ini$neta1) & (ini$neta1 == ini$neta2) & is.na(ini$name)]
    if (length(unnamed.etas) > 0) {
      stop(sprintf("The following ETAs are unnamed: %s", paste(sprintf("ETA[%d]", unnamed.etas), collapse = ", ")))
    }
    theta.names <- ini$theta.names
    eta.names <- ini$eta.names
  }
  errn <- 0

  any.theta.names <- function(what, th.names) {
    return(any(what[1] == th.names))
  }

  find.theta <- function(x) {
    if (length(x) == 0) {
    } else if (is.name(x) && length(x) == 1 && any.theta.names(as.character(x), theta.names)) {
      return(as.character(x))
    } else if (is.call(x) || is.pairlist(x)) {
      if (length(x) == 0) {
      } else if (length(x) == 1 && any.theta.names(as.character(x), theta.names)) {
        return(as.character(x))
      } else if (identical(x[[1]], quote(`+`))) {
        th <- c()
        if (length(x) >= 3) {
          if (any.theta.names(as.character(x[[3]]), theta.names)) {
            th <- as.character(x[[3]])
          }
        }
        if (length(x) >= 2) {
          if (length(x[[2]]) > 1) {
            return(c(th, find.theta(x[[2]])))
          } else {
            if (any.theta.names(as.character(x[[2]]), theta.names)) {
              th <- c(th, as.character(x[[2]]))
            }
            return(th)
          }
        }
      }
    }
  }
  .regPar <- rex::rex("nlmixr_", capture(anything), "_par")
  .doDist <- function(distName, distArgs, curCond = NULL) {
    if (any(regexpr(.regPar, distArgs) != -1)) {
      .tmp <- distArgs[regexpr(.regPar, distArgs) != -1]
      distArgs <- gsub(.regPar, "\\1", distArgs)
    }
    if (!any(names(dists) == distName)) {
      stop(sprintf("The %s distribution is currently unsupported.", distName))
    }
    .nargs <- dists[[distName]]
    if (length(.nargs) == 1) {
      if (length(distArgs) != .nargs) {
        stop(sprintf("The %s distribution requires %s argument%s.", distName, .nargs, ifelse(.nargs == 1, "", "s")))
      }
    } else {
      .minNargs <- min(.nargs)
      .maxNargs <- max(.nargs)
      if (length(distArgs) < .minNargs | length(distArgs) > .maxNargs) {
        stop(sprintf("The %s distribution requires %s-%s arguments.", distName, min(.nargs), max(.nargs)))
      }
    }
    if (length(distArgs) == 0L) {
      .tmp <- .as.data.frame(bounds)
      .tmp1 <- .tmp[1, ]
      .tmp1[1, ] <- NA
      .tmp1[, "err"] <- distName
      .tmp1[, "est"] <- -Inf
      if (!is.null(curCond)) {
        .tmp1[, "condition"] <- .deparse(curCond)
      } else {
        .tmp1[, "condition"] <- ""
      }
      .tmp <- rbind(.tmp, .tmp1)
      class(.tmp) <- c("nlmixrBounds", "data.frame")
      assign("bounds", .tmp, this.env)
      return(errn)
    }
    for (.i in seq_along(distArgs)) {
      .tmp <- suppressWarnings(as.numeric(distArgs[.i]))
      errn <- errn + 1
      if (!is.na(.tmp)) {
        ## FIXME: allow numeric estimates...?
        stop("Distribution parameters cannot be numeric, but need to be estimated")
      }
      .w <- which(bounds$name == distArgs[.i])
      if (length(.w) == 0) {
        stop("Residual distribution parameter(s) estimates were not found in ini block")
      }
      .tmp <- .as.data.frame(bounds)
      .tmp$err[.w] <- ifelse(.i == 1, distName, paste0(distName, .i))
      if (any(distName == distsPositive)) {
        if (any(is.na(.tmp$lower[.w])) || any(is.infinite(.tmp$lower[.w]))) {
          .tmp$lower[.w] <- 0
        }
        if ((.tmp$lower[.w] < 0) || .tmp$est[.w] < 0 || .tmp$upper[.w] < 0) {
          stop(sprintf("The distribution '%s' must have positive parameter estimates", distName))
        }
      }
      if (!is.null(curCond)) {
        .tmp$condition[.w] <- sub(rex::rex(or("cmt", "CMT"), any_spaces, "==", any_spaces), "", .deparse(curCond))
      } else {
        .tmp$condition[.w] <- ""
      }
      class(.tmp) <- c("nlmixrBounds", "data.frame")
      assign("bounds", .tmp, this.env)
    }
    return(errn)
  }
  .predDf <- NULL
  .doDist1 <- function(err1, err1.v, err1.args, x2, x3, curCond = NULL) {
    .bCond <- is.null(curCond)
    if (is.null(curCond)) curCond <- x2
    if (!is.na(suppressWarnings(as.numeric(err1.v)))) {
      stop("Distribution parameters cannot be numeric, but need to be estimated.")
    }
    assign("errs.specified", unique(errs.specified, err1), this.env)
    if (any(do.pred == c(2, 4, 5))) {
      return(quote(nlmixrIgnore()))
    }
    else if (any(do.pred == c(1, 4, 5))) {
      assign(".predDf", rbind(
        .predDf,
        .data.frame(
          cond = ifelse(.bCond, "", sub(rex::rex(or("cmt", "CMT"), any_spaces, "==", any_spaces), "", .deparse(curCond))),
          var = .deparse(x2)
        )
      ), this.env)
      return(bquote(if (CMT == .(curCond)) {
        nlmixr_pred <- .(x2)
      }))
    } else if (do.pred == 3) {
      .doDist(err1, err1.args, curCond)
      tmp <- bounds
      if ((any(paste(tmp$err) == "add") || any(paste(tmp$err) == "norm") || any(paste(tmp$err) == "dnorm") ||
        any(paste(tmp$err) == "lnorm") || any(paste(tmp$err) == "dlnorm") || any(paste(tmp$err) == "logn") ||
        any(paste(tmp$err) == "dlogn")) &&
        any(paste(tmp$err) == "prop")) {
        assign("errn", errn + 1, this.env)
        assign("add.prop.errs", rbind(
          add.prop.errs,
          .data.frame(y = sprintf("Y%02d", errn), add = TRUE, prop = TRUE)
        ), this.env)
      } else if (any(paste(tmp$err) == "prop")) {
        assign("errn", errn + 1, this.env)
        assign("add.prop.errs", rbind(
          add.prop.errs,
          .data.frame(y = sprintf("Y%02d", errn), add = FALSE, prop = TRUE)
        ), this.env)
      }
      return(bquote(return(.(sprintf("Y%02d", errn)))))
    } else {
      return(bquote(if (CMT == .(curCond)) {
        return(.(x3))
      }))
    }
  }
  .doDist2 <- function(err1, err1.v, err1.args, err2, err2.v, err2.args, x2, x3, curCond = NULL) {
    .bCond <- is.null(curCond)
    if (is.null(curCond)) curCond <- x2
    if (!is.na(suppressWarnings(as.numeric(err1.v)))) {
      stop("Distribution parameters cannot be numeric, but need to be estimated.")
    }
    if (!is.na(suppressWarnings(as.numeric(err2.v)))) {
      stop("Distribution parameters cannot be numeric, but need to be estimated.")
    }
    if (any(err1 == add.dists) &&
      any(err2 == add.dists)) {
      tmp <- paste(sort(c(err1, err2)), collapse = "+")
      assign("errs.specified", unique(errs.specified, tmp), this.env)
      if (any(do.pred == c(2, 4, 5))) {
        return(quote(nlmixrIgnore()))
      }
      else if (any(do.pred == c(1, 4, 5))) {
        assign(".predDf", rbind(
          .predDf,
          .data.frame(
            cond = ifelse(.bCond, "", sub(rex::rex(or("cmt", "CMT"), any_spaces, "==", any_spaces), "", .deparse(curCond))),
            var = .deparse(x2)
          )
        ), this.env)
        return(bquote(if (CMT == .(curCond)) {
          nlmixr_pred <- .(x2)
        }))
      } else if (do.pred == 3) {
        .doDist(err1, err1.args, curCond)
        .doDist(err2, err2.args, curCond)
        tmp <- bounds
        if ((any(paste(tmp$err) == "add") ||
          any(paste(tmp$err) == "norm") ||
          any(paste(tmp$err) == "dnorm") ||
          any(paste(tmp$err) == "lnorm") ||
          any(paste(tmp$err) == "dlnorm") ||
          any(paste(tmp$err) == "logn") ||
          any(paste(tmp$err) == "dlogn")
        ) && any(paste(tmp$err) == "prop")) {
          assign("errn", errn + 1, this.env)
          assign("add.prop.errs", rbind(
            add.prop.errs,
            .data.frame(y = sprintf("Y%02d", errn), add = TRUE, prop = TRUE)
          ), this.env)
        }
        return(bquote(return(.(sprintf("Y%02d", errn)))))
      } else {
        return(bquote(if (CMT == .(curCond)) {
          return(.(x3))
        }))
      }
    } else {
      stop(sprintf(
        "The %s and %s distributions cannot be combined\nCurrently can combine: %s",
        as.character(x3[[2]][[1]]), as.character(x3[[3]][[1]]),
        paste(add.dists, collapse = ", ")
      ))
    }
  }
  .doDist3 <- function(err1, err1.v, err1.args, err2, err2.v, err2.args, err3, err3.v, err3.args, x2, x3, curCond = NULL) {
    .bCond <- is.null(curCond)
    if (is.null(curCond)) curCond <- x2
    if (!is.na(suppressWarnings(as.numeric(err1.v)))) {
      stop("Distribution parameters cannot be numeric, but need to be estimated.")
    }
    if (!is.na(suppressWarnings(as.numeric(err2.v)))) {
      stop("Distribution parameters cannot be numeric, but need to be estimated.")
    }
    if (!is.na(suppressWarnings(as.numeric(err3.v)))) {
      stop("Distribution parameters cannot be numeric, but need to be estimated.")
    }
    if (any(err1 == add.dists) &&
      any(err2 == add.dists) &&
      any(err3 == add.dists)) {
      tmp <- paste(sort(c(err1, err2, err3)), collapse = "+")
      assign("errs.specified", unique(errs.specified, tmp), this.env)
      if (any(do.pred == c(2, 4, 5))) {
        return(quote(nlmixrIgnore()))
      }
      else if (any(do.pred == c(1, 4, 5))) {
        assign(".predDf", rbind(
          .predDf,
          .data.frame(
            cond = ifelse(.bCond, "", sub(rex::rex(or("cmt", "CMT"), any_spaces, "==", any_spaces), "", .deparse(curCond))),
            var = .deparse(x2)
          )
        ), this.env)
        return(bquote(if (CMT == .(curCond)) {
          nlmixr_pred <- .(x2)
        }))
      } else if (do.pred == 3) {
        .doDist(err1, err1.args, curCond)
        .doDist(err2, err2.args, curCond)
        .doDist(err3, err3.args, curCond)
        tmp <- bounds
        if ((any(paste(tmp$err) == "add") ||
          any(paste(tmp$err) == "dnorm") ||
          any(paste(tmp$err) == "norm") ||
          any(paste(tmp$err) == "lnorm") ||
          any(paste(tmp$err) == "dlnorm") ||
          any(paste(tmp$err) == "logn") ||
          any(paste(tmp$err) == "dlogn")
        ) && any(paste(tmp$err) == "prop")) {
          assign("errn", errn + 1, this.env)
          assign("add.prop.errs", rbind(
            add.prop.errs,
            .data.frame(y = sprintf("Y%02d", errn), add = TRUE, prop = TRUE)
          ), this.env)
        }
        return(bquote(return(.(sprintf("Y%02d", errn)))))
      } else {
        return(bquote(if (CMT == .(curCond)) {
          return(.(x3))
        }))
      }
    } else {
      stop(sprintf(
        "The %s, %s and %s distributions cannot be combined\nCurrently can combine: %s",
        err1, err2, err3, paste(add.dists, collapse = ", ")
      ))
    }
  }
  f <- function(x) {
    if (is.name(x)) {
      if (any.theta.names(as.character(x), theta.names)) {
        assign("theta.ord", unique(c(theta.ord, as.character(x))), this.env)
      }
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`~`)) &&
        as.character(x[[3]][[1]]) == "|" &&
        any(as.character(x[[3]][[2]][[1]]) == c(names(dists), unsupported.dists))) {
        ch.dist <- as.character(x[[3]][[2]])
        if (length(x[[3]][[3]]) == 1) {
          curCond <- sprintf("%s", as.character(x[[3]][[3]]))
        } else {
          curCond <- .deparse(x[[3]][[3]])
        }
        curCond <- eval(parse(text = sprintf("quote(%s)", curCond)))
        if (length(ch.dist) == 1) {
          err1 <- as.character(x[[3]][[2]][[1]])
          err1.v <- err1
          err1.args <- character(0)
          .doDist1(err1, err1.v, err1.args, x[[2]], x[[3]][[2]], curCond = curCond)
        } else {
          err1 <- as.character(x[[3]][[2]][[1]])
          err1.v <- as.character(x[[3]][[2]][[2]])
          err1.args <- as.character(x[[3]][[2]][-1])
          .doDist1(err1, err1.v, err1.args, x[[2]], x[[3]][[2]], curCond = curCond)
        }
      } else if (identical(x[[1]], quote(`~`)) &&
        as.character(x[[3]][[1]]) == "|" &&
        identical(x[[3]][[2]][[1]], quote(`+`)) &&
        length(as.character(x[[3]][[2]])) == 3 &&
        any(as.character(x[[3]][[2]][[2]][[1]]) == c(names(dists), unsupported.dists)) &&
        any(as.character(x[[3]][[2]][[3]][[1]]) == c(names(dists), unsupported.dists))) {
        err1 <- as.character(x[[3]][[2]][[2]][[1]])
        err1.v <- as.character(x[[3]][[2]][[2]][[2]])
        err1.args <- as.character(x[[3]][[2]][[2]][-1])
        err2 <- as.character(x[[3]][[2]][[3]][[1]])
        err2.v <- as.character(x[[3]][[2]][[3]][[2]])
        err2.args <- as.character(x[[3]][[2]][[3]][-1])
        if (length(x[[3]][[3]]) == 1) {
          curCond <- sprintf("%s", as.character(x[[3]][[3]]))
        } else {
          curCond <- .deparse(x[[3]][[3]])
        }
        curCond <- eval(parse(text = sprintf("quote(%s)", curCond)))
        .doDist2(err1, err1.v, err1.args, err2, err2.v, err2.args, x[[2]], x[[3]][[2]], curCond = curCond)
      } else if (identical(x[[1]], quote(`~`)) &&
        as.character(x[[3]][[1]]) == "|" &&
        identical(x[[3]][[2]][[1]], quote(`+`)) &&
        length(as.character(x[[3]][[2]])) == 3 &&
        any(as.character(x[[3]][[2]][[2]][[2]][[1]]) == c(names(dists), unsupported.dists))) {
        err1 <- as.character(x[[3]][[2]][[2]][[2]][[1]])
        err1.v <- as.character(x[[3]][[2]][[2]][[2]][[2]])
        err1.args <- as.character(x[[3]][[2]][[2]][[2]][-1])
        err2 <- as.character(x[[3]][[2]][[3]][[1]])
        err2.v <- as.character(x[[3]][[2]][[3]][[2]])
        err2.args <- as.character(x[[3]][[2]][[3]][-1])
        err3 <- as.character(x[[3]][[2]][[2]][[3]][[1]])
        err3.v <- as.character(x[[3]][[2]][[2]][[3]][[2]])
        err3.args <- as.character(x[[3]][[2]][[2]][[3]][-1])
        if (length(x[[3]][[3]]) == 1) {
          curCond <- sprintf("%s", as.character(x[[3]][[3]]))
        } else {
          curCond <- .deparse(x[[3]][[3]])
        }
        curCond <- eval(parse(text = sprintf("quote(%s)", curCond)))
        .doDist3(err1, err1.v, err1.args, err2, err2.v, err2.args, err3, err3.v, err3.args, x[[2]], x[[3]][[2]], curCond = curCond)
      } else if (identical(x[[1]], quote(`~`)) &&
        any(as.character(x[[3]][[1]]) == c(names(dists), unsupported.dists))) {
        ch.dist <- as.character(x[[3]])
        if (length(ch.dist) == 1) {
          err1 <- as.character(x[[3]][[1]])
          err1.v <- err1
          err1.args <- character(0)
          .doDist1(err1, err1.v, err1.args, x[[2]], x[[3]], curCond = NULL)
        } else {
          err1 <- as.character(x[[3]][[1]])
          err1.v <- as.character(x[[3]][[2]])
          err1.args <- as.character(x[[3]][-1])
          .doDist1(err1, err1.v, err1.args, x[[2]], x[[3]], curCond = NULL)
        }
      } else if (identical(x[[1]], quote(`~`)) && ## Arg parsing 4 should be the last....
        identical(x[[3]][[1]], quote(`+`)) &&
        identical(x[[3]][[2]][[1]], quote(`+`)) &&
        identical(x[[3]][[2]][[2]][[1]], quote(`+`)) &&
        any(as.character(x[[3]][[2]][[2]][[2]][[1]]) == c(names(dists), unsupported.dists))) {
        stop(sprintf(
          "Only 3 distributions can be combined.\nCurrently can combine: %s",
          paste(add.dists, collapse = ", ")
        ))
      } else if (identical(x[[1]], quote(`~`)) &&
        identical(x[[3]][[1]], quote(`+`)) &&
        identical(x[[3]][[2]][[1]], quote(`+`)) &&
        any(as.character(x[[3]][[2]][[2]][[1]]) == c(names(dists), unsupported.dists))) {
        err1 <- as.character(x[[3]][[2]][[2]][[1]])
        err1.v <- as.character(x[[3]][[2]][[2]][[2]])
        err1.args <- as.character(x[[3]][[2]][[2]][-1])
        err2 <- as.character(x[[3]][[3]][[1]])
        err2.v <- as.character(x[[3]][[3]][[2]])
        err2.args <- as.character(x[[3]][[3]][-1])
        err3 <- as.character(x[[3]][[2]][[3]][[1]])
        err3.v <- as.character(x[[3]][[2]][[3]][[2]])
        err3.args <- as.character(x[[3]][[2]][[3]][-1])
        .doDist3(err1, err1.v, err1.args, err2, err2.v, err2.args, err3, err3.v, err3.args, x[[2]], x[[3]], curCond = NULL)
      } else if (identical(x[[1]], quote(`~`)) &&
        identical(x[[3]][[1]], quote(`+`)) &&
        length(as.character(x[[3]])) == 3 &&
        any(as.character(x[[3]][[2]][[1]]) == c(names(dists), unsupported.dists)) &&
        any(as.character(x[[3]][[3]][[1]]) == c(names(dists), unsupported.dists))) {
        err1 <- as.character(x[[3]][[2]][[1]])
        err1.v <- as.character(x[[3]][[2]][[2]])
        err1.args <- as.character(x[[3]][[2]][-1])
        err2 <- as.character(x[[3]][[3]][[1]])
        err2.v <- as.character(x[[3]][[3]][[2]])
        err2.args <- as.character(x[[3]][[3]][-1])
        .doDist2(err1, err1.v, err1.args, err2, err2.v, err2.args, x[[2]], x[[3]])
      } else if (identical(x[[1]], quote(`~`)) && (do.pred != 2)) {
        return(quote(nlmixrIgnore()))
      } else if (identical(x[[1]], quote(`<-`)) && !any(do.pred == c(2, 4, 5))) {
        return(quote(nlmixrIgnore()))
      } else if (identical(x[[1]], quote(`=`)) && !any(do.pred == c(2, 4, 5))) {
        return(quote(nlmixrIgnore()))
      } else if (identical(x[[1]], quote(`<-`)) && do.pred == 4) {
        ## SAEM requires = instead of <-
        x[[1]] <- quote(`=`)
        return(as.call(lapply(x, f)))
      } else if (identical(x[[1]], quote(`exp`)) && any(do.pred == c(4, 5))) {
        ## Need traverse the parsing tree to get log theta/eta
        ## parameters.
        find.log <- function(x) {
          if (is.atomic(x) || is.name(x)) {
            if (any.theta.names(as.character(x), theta.names)) {
              assign("log.theta", unique(c(log.theta, as.character(x))), this.env)
            } else if (any.theta.names(as.character(x), eta.names)) {
              assign("log.eta", unique(c(log.eta, as.character(x))), this.env)
            }
            return(x)
          } else if (is.pairlist(x)) {
            return(lapply(x, find.log))
          } else if (is.call(x)) {
            return(lapply(x, find.log))
          } else {
            stop("Don't know how to handle type ", typeof(x),
              call. = FALSE
            )
          }
        }
        find.log(x[[2]])
        if (length(x[[2]]) == 1 && any.theta.names(as.character(x[[2]]), theta.names)) {
          tmp <- as.character(x[[2]])
          .oneTheta <<- unique(c(.oneTheta, tmp))
        }
        return(as.call(lapply(x, f)))
      } else if (identical(x[[1]], quote(`+`)) ||
        identical(x[[1]], quote(`-`))) {
        isMinus <- identical(x[[1]], quote(`-`))
        ## print(as.character(x))
        if (any(do.pred == c(4, 5))) {
          if (length(x) >= 3) {
            ## message("---")
            ## print(x[[1]])
            ## print(x[[2]]);
            ## print(x[[3]]);
            wm <- NULL
            muRef <- FALSE
            if (length(x[[3]]) == 3) {
              if (identical(x[[3]][[1]], quote(`*`))) {
                wm <- 3
                wm0 <- 2
                muRef <- TRUE
              } else if (identical(x[[3]][[1]], quote(`/`)) ||
                identical(x[[3]][[1]], quote(`^`)) ||
                identical(x[[3]][[1]], quote(`**`))) {
                wm <- 3
                wm0 <- 2
              }
            }
            if (length(x[[2]]) == 3) {
              if (identical(x[[2]][[1]], quote(`*`))) {
                wm <- 2
                wm0 <- 3
                muRef <- TRUE
              } else if (identical(x[[2]][[1]], quote(`/`)) ||
                identical(x[[2]][[1]], quote(`^`)) ||
                identical(x[[2]][[1]], quote(`**`))) {
                wm <- 3
                wm0 <- 2
              }
            }
            if (!is.null(wm)) {
              cur <- 2
              th <- 3
              w <- try(which(all.covs == as.character(x[[wm]][[2]])[1]), silent = TRUE)
              if (inherits(w, "try-error")) w <- integer(0)
              if (length(w) == 0) {
                cur <- 3
                th <- 2
                w <- try(which(all.covs == as.character(x[[wm]][[3]])[1]), silent = TRUE)
                if (inherits(w, "try-error")) w <- integer(0)
              }
              if (length(w) == 1) {
                cov <- all.covs[w]
                th <- as.character(x[[wm]][[th]])
                th0 <- find.theta(x[[wm0]])
                if (muRef && !isMinus && length(th0) == 1 && do.pred == 4) {
                  tmp <- get("cov.ref", this.env)
                  tmp[[cov]] <- c(tmp[[cov]], structure(th0, .Names = th))
                  assign("cov.ref", tmp, this.env)
                  return(f(x[[wm0]]))
                }
                tmp <- get("cov.theta", this.env)
                assign("cov.theta", unique(c(tmp, th)), this.env)
              }
            }
            if (any.theta.names(as.character(x[[2]]), eta.names) &&
              any.theta.names(as.character(x[[3]]), theta.names)) {
              ## Found ETA+THETA
              tmp <- .mu.ref
              tmp[[as.character(x[[2]])]] <- as.character(x[[3]])
              ## assign("mu.ref", tmp, this.env);
              .mu.ref <<- tmp
              tmp <- as.character(x[[3]])
              .oneTheta <<- unique(c(.oneTheta, tmp))
              ## Collapse to THETA
              return(x[[3]])
            } else if (any.theta.names(as.character(x[[3]]), eta.names) &&
              any.theta.names(as.character(x[[2]]), theta.names)) {
              ## Found THETA+ETA
              tmp <- .mu.ref
              tmp[[as.character(x[[3]])]] <- as.character(x[[2]])
              ## assign(".mu.ref", tmp, this.env)
              .mu.ref <<- tmp
              tmp <- as.character(x[[2]])
              .oneTheta <<- unique(c(.oneTheta, tmp))
              ## Collapse to THETA
              ## model$omega=diag(c(1,1,0))
              ## 0 is not estimated.
              ## inits$omega has the initial estimate
              ## mod$res.mod = 1 = additive or poisson
              ## mod$res.mod = 2 = proportional
              ## mod$res.mod = 3 = additive + proportional
              ## a+b*f
              ## mod$ares = initial estimate of res
              ## mod$bres = initial estimate of
              return(x[[2]])
            } else if (any.theta.names(as.character(x[[3]]), eta.names) &&
              length(x[[2]]) > 1) {
              ## This allows 123 + Cl + 123 + eta.Cl + 123
              ## And collapses to 123 + Cl + 123 + 123
              ## Useful for covariates...
              eta <- as.character(x[[3]])
              th <- find.theta(x[[2]])
              if (length(th) == 1) {
                tmp <- .mu.ref
                tmp[[eta]] <- th
                ## assign("tmp", .mu.ref, this.env)
                .mu.ref <<- tmp
                tmp <- as.character(th)
                .oneTheta <<- unique(c(.oneTheta, tmp))
                return(f(as.call(x[[2]])))
              }
            } else if (length(x) < 3) {
            } else if (any.theta.names(as.character(x[[3]]), theta.names) &&
              length(x[[2]]) > 1) {
              ## This allows 123 + eta.Cl + 123 + Cl + 123
              ## And collapses to 123  + 123 + Cl + 123
              ## Useful for covariates...
              theta <- as.character(x[[3]])
              .etas <- c()
              find.etas <- function(x) {
                if (is.atomic(x) || is.name(x)) {
                  return(x)
                } else if (is.pairlist(x)) {
                  return(lapply(x, find.etas))
                } else if (is.call(x)) {
                  if (identical(x[[1]], quote(`+`)) &&
                    any.theta.names(as.character(x[[3]]), eta.names)) {
                    .etas <<- c(.etas, as.character(x[[3]]))
                    return(x[[2]])
                  }
                  return(as.call(lapply(x, find.etas)))
                } else {
                  stop("Don't know how to handle type ", typeof(x),
                    call. = FALSE
                  )
                }
              }
              new <- find.etas(x[[2]])
              if (length(.etas) == 1) {
                tmp <- .mu.ref
                tmp[[.etas]] <- theta
                ## assign("tmp", .mu.ref, this.env);
                .mu.ref <<- tmp
                tmp <- as.character(theta)
                .oneTheta <<- unique(c(.oneTheta, tmp))
                x[[2]] <- new
              }
            }
          }
        }
        return(as.call(lapply(x, f)))
      } else {
        return(as.call(lapply(x, f)))
      }
    } else if (is.pairlist(x)) {
      as.pairlist(lapply(x, f))
    } else if (is.atomic(x)) {
      return(x)
    } else {
      stop("Don't know how to handle type ", typeof(x),
        call. = FALSE
      )
    }
  }
  rm.empty <- function(x) {
    ## empty if/else
    ## First remove if () followed by nlmixrIgnore()
    ignoreit <- rex::rex(any_spaces, "nlmixrIgnore()", any_spaces)
    w1 <- which(regexpr(rex::rex(start, any_spaces, or(group("if", any_spaces, "(", anything, ")"), "else"), any_spaces, end), x) != -1)
    if (length(w1) > 0) {
      w2 <- which(regexpr(ignoreit, x[w1 + 1]) != -1)
      if (length(w2) > 0) {
        x <- x[-w1[w2]]
      }
    }
    x <- x[regexpr(ignoreit, x, perl = TRUE) == -1]
    w1 <- which(regexpr(rex::rex(start, any_spaces, or("if", "else"), anything, "{", end), x) != -1)
    if (length(w1) > 0) {
      w2 <- w1 + 1
      w3 <- which(regexpr(rex::rex(start, any_spaces, "}", end), x[w2]) != -1)
      if (length(w3) > 0) {
        w1 <- w1[w3]
        w2 <- w2[w3]
        return(x[-c(w1, w2)])
      } else {
        return(x)
      }
    } else {
      return(x)
    }
  }
  new.fn <- function(x) {
    x <- rm.empty(x)
    if (do.pred == 2) {
      assign("rxode", any(regexpr(rex::rex(start, any_spaces, "d/dt(", anything, ")", any_spaces, or("=", "<-")), x) != -1), this.env)
    }
    x <- c("function(){", .bodyDewrap(x), "}")
    x <- eval(parse(text = paste(x, collapse = "\n")))
    return(x)
  }
  .fun00 <- function(funTxt) {
    .fun0 <- function(reg0 = rex::rex(
                        capture(any_spaces),
                        capture(except_any_of("()\n; ")),
                        capture(any_spaces)
                      ),
                      reg00 = "", reg01 = "", repE = "expr") {
      if (any(regexpr(.regRx, funTxt[1:w]) != -1)) {
        .rxBegin <- funTxt[1:w]
        .re <- rex::rex(
          start, capture(any_spaces), reg0,
          capture(any_spaces, or("=", "~", "<-")),
          capture(anything)
        )
        .w <- which(regexpr(.re, .rxBegin) != -1)
        .lines <- funTxt[.w]
        ## Extract any estimated parameter only expressions and promote them.
        .w2 <- which(regexpr(rex::rex(anything, or("=", "~", "<-"), .lhsReg), .lines) != -1)
        .doIt <- TRUE
        if (length(.w2) > 0) {
          ## Remove any parameters that depend on prior values and not covariates/expressions.
          .w <- .w[-.w2]
          if (length(.w) == 0) {
            .doIt <- FALSE
          } else {
            .lines <- funTxt[.w]
          }
        }
        if (.doIt) {
          ## Now remove any variables that are duplicated.  This often happens when
          ## if (x){y=a;}else{y=b;} blocks.
          .rxBegin <- funTxt[.w]
          .dups <- funTxt[regexpr(rex::rex(or("=", "~", "<-")), funTxt) != -1]
          .dups <- gsub(rex::rex(
            start, any_spaces, reg0, any_spaces,
            or("=", "~", "<-"), anything
          ), "\\2", .dups)
          .dups <- unique(.dups[duplicated(.dups)])
          .w2 <- which(regexpr(rex::rex(
            start, any_spaces,
            reg00, or(.dups), reg01, any_spaces,
            or("=", "~", "<-"), anything
          ), .lines) != -1)
          if (length(.w2) > 0) {
            .w <- .w[-.w2]
            if (length(.w) == 0) {
              .doIt <- FALSE
            } else {
              .lines <- funTxt[.w]
            }
          }
          if (.doIt) {
            .rxBegin <- funTxt
            .rxBegin[.w] <- gsub(.re, paste0("\\1\\2\\3\\4\\5nlmixr_\\3_", repE),
              .rxBegin[.w],
              perl = TRUE
            )
            .lines <- gsub(
              rex::rex(.re, anything),
              paste0("\\1nlmixr_\\3_", repE, "\\5\\6"), .lines
            )
            ## .w <- min(which(regexpr(.regRx, .rxBegin, perl=TRUE) != -1))-1;
            ## return(c(.rxBegin[1:.w],.lines,.rxBegin[-(1:.w)]));
            return(c(.lines, .rxBegin))
          }
        }
      }
      return(funTxt)
    }
    w <- which(regexpr(reg, funTxt, perl = TRUE) != -1)
    w <- max(w)
    funTxt <- .fun0(
      reg0 = rex::rex(capture(or("f(", "F(")), capture(except_any_of("()\n; ")), capture(")")),
      reg00 = rex::rex(or("f(", "F(")), reg01 = ")", repE = "F"
    )
    w <- which(regexpr(reg, funTxt, perl = TRUE) != -1)
    w <- max(w)
    funTxt <- .fun0(
      reg0 = rex::rex(capture(or("dur(", "d(")), capture(except_any_of("()\n; ")), capture(")")),
      reg00 = rex::rex(or("dur(", "d(")), reg01 = ")", repE = "dur"
    )
    w <- which(regexpr(reg, funTxt, perl = TRUE) != -1)
    w <- max(w)
    funTxt <- .fun0(
      reg0 = rex::rex(capture(or("lag(", "alag(")), capture(except_any_of("()\n; ")), capture(")")),
      reg00 = rex::rex(or("lag(", "alag(")), reg01 = ")", repE = "lag"
    )
    w <- which(regexpr(reg, funTxt, perl = TRUE) != -1)
    w <- max(w)
    funTxt <- .fun0(
      reg0 = rex::rex(capture(or("r(", "rate(")), capture(except_any_of("()\n; ")), capture(")")),
      reg00 = rex::rex(or("r(", "rate(")), reg01 = ")", repE = "rate"
    )
    w <- which(regexpr(reg, funTxt, perl = TRUE) != -1)
    w <- max(w)

    funTxt <- .fun0(
      reg0 = rex::rex(capture(any_spaces), capture(except_any_of("()\n; ")), capture("(0)")),
      reg00 = "", reg01 = "(0)", repE = "ini"
    )
    w <- which(regexpr(reg, funTxt, perl = TRUE) != -1)
    w <- max(w)

    funTxt <- .fun0()
    w <- which(regexpr(reg, funTxt, perl = TRUE) != -1)
    w <- max(w)

    .finalFix <- function() {
      if (any(regexpr(.regRx, funTxt[1:w],
        perl = TRUE
      ) != -1)) {
        ## There are still mixed PK parameters and ODEs
        w <- which(regexpr(.regRx, funTxt, perl = TRUE) != -1)
        w <- min(w) - 1
        if (w > 0) {
          .env <- new.env(parent = emptyenv())
          assign("extra", NULL, .env)
          .subs <- function(x) {
            if (is.atomic(x)) {
              x
            } else if (is.name(x)) {
              if (any(as.character(x) == ini$name)) {
                assign("extra", unique(c(.env$extra, as.character(x))), .env)
                return(eval(parse(text = sprintf("quote(nlmixr_%s_par)", as.character(x)))))
              } else {
                return(x)
              }
            } else if (is.call(x)) {
              as.call(lapply(x, .subs))
            } else if (is.pairlist(x)) {
              as.pairlist(lapply(x, .subs))
            } else {
              return(x)
            }
          }
          .x <- .deparse(.subs(body(eval(parse(text = paste("function(){\n", paste(funTxt[-(1:w)], collapse = "\n"), "\n}"))))))
          .x <- .x[-1]
          .x <- .x[-length(.x)]
          return(c(funTxt[1:w], paste0("nlmixr_", .env$extra, "_par <- ", .env$extra), .x))
        }
        return(funTxt)
      }
      return(funTxt)
    }
    funTxt <- .finalFix()
    return(funTxt)
  }
  .tmp <- .deparse1(body(fun))
  ## Assume ~ is boundaries
  reg <- rex::rex(or("=", "<-"), anything, boundary, or(ini$name), boundary)
  .regRx <- rex::rex(start, any_spaces, or(
    "d/dt(", "f(", "F(", "dur(", "d(",
    "lag(", "alag(", "r(", "rate(",
    group(anything, any_spaces, "(0)", any_spaces, or("=", "<-"))
  ))
  .lhs <- nlmixrfindLhs(body(
    eval(parse(text = paste(
      "function(){",
      paste(.tmp, collapse = "\n"),
      "}"
    )))
  ))
  .lhsReg <- rex::rex(boundary, or(.lhs), boundary)
  fun <- eval(parse(text = paste(c("function()({", .fun00(.tmp), "})"), collapse = "\n")))
  all.covs <- character()
  do.pred <- 1
  pred.txt <- .deparse(f(body(fun)))
  .reg <- rex::rex(or(
    group(any_spaces, "(", any_spaces, "{", any_spaces),
    group(any_spaces, "}", any_spaces, ")", any_spaces),
    group(any_spaces, "nlmixrIgnore()", any_spaces),
    group(any_spaces, "(", any_spaces, "{", any_spaces, "nlmixrIgnore()", any_spaces),
    group("nlmixrIgnore()", any_spaces, "}", any_spaces, ")")
  ))
  .pred <- sapply(pred.txt, function(x) {
    regexpr(.reg, x) != -1
  })
  if (all(.pred)) {
    stop("There must be at least one prediction in the model({}) block.  Use `~` for predictions")
  }
  pred <- new.fn(pred.txt)
  do.pred <- 0
  err <- new.fn(.deparse(f(body(fun))))
  do.pred <- 2
  rest.txt <- .deparse(f(body(fun)))
  rest <- new.fn(rest.txt)
  rest.funs <- allCalls(body(rest))
  rest.vars <- allVars(body(rest))
  all.covs <- setdiff(rest.vars, paste0(bounds$name))
  do.pred <- 3
  grp.fn <- new.fn(.deparse(f(body(fun))))
  do.pred <- 4
  saem.pars <- try(.deparse(f(body(fun))), silent = TRUE)
  nlme.mu.fun2 <- NULL
  if (inherits(saem.pars, "try-error")) {
    saem.pars <- NULL
  }
  do.pred <- 5
  nlme.mu.fun <- try(.deparse(f(body(fun))), silent = TRUE)
  if (inherits(nlme.mu.fun, "try-error")) {
    nlme.mu.fun <- NULL
  }
  .pred <- FALSE
  if (!rxode) {
    rxode <- TRUE
    .pred <- TRUE
  }
  .linCmt <- FALSE
  if (any(regexpr(rex::rex("limCmt("), .deparse(body(fun))) != -1)) {
    stop("You used `limCmt`, did you mean `linCmt`?")
  }
  if (any(regexpr(rex::rex("linCmt("), .deparse(body(fun))) != -1)) {
    .linCmt <- TRUE
    .hasLinCmt <- any(regexpr(rex::rex("linCmt("), .deparse(body(rest))) == -1)
    rx.txt <- .deparse1(body(rest))
    .regLin <- rex::rex(
      start, any_spaces,
      capture(
        or(
          group(one_of("Kk"), some_of("AaEe0123456789")),
          group(
            one_of("V", "v"),
            any_of("c", "C", "P", "p", "T", "t", "S", "s", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
          ),
          group(one_of("Qq"), any_of("0":"9", "p", "c")),
          group(one_of("Cc"), one_of("Ll"), any_of("Dd2")),
          "aob", "AOB", "alpha", "ALPHA", "Alpha",
          "beta", "BETA", "Beta"
        )
      ), any_spaces, or("=", "<-"),
      capture(anything)
    )
    .pred <- gsub(
      .regLin,
      "nlmixr_lincmt_\\1 <- \\2", rx.txt
    )
    .pred <- .pred[regexpr(rex::rex(
      start, any_spaces, or("dur", "d", "rate", "r", "lag", "alag", "f", "F"),
      any_spaces, "(",
      any_spaces, or("depot", "central"), any_spaces,
      ")", any_spaces, or("=", "<-"), any_spaces, capture(anything)
    ), .pred) == -1]
    .vars <- gsub(.regLin, "\\1", rx.txt[regexpr(.regLin, rx.txt) != 0])
    .rx <- gsub(
      rex::rex(
        start, any_spaces, capture(or(.vars)), any_spaces, or("=", "<-"),
        capture(anything)
      ),
      "\\1 <- nlmixr_lincmt_\\1", rx.txt
    )
    .pred <- .pred[regexpr(rex::rex("linCmt("), .pred) == -1]
    if (any(regexpr(rex::rex("linCmt("), .rx) != -1)) .hasLinCmt <- FALSE
    rest <- eval(parse(text = paste(c(
      "function(){",
      .pred, .rx,
      ifelse(.hasLinCmt, "nlmixr_lincmt_pred <- linCmt()\n}", "}")
    ),
    collapse = "\n"
    )))
    if (!is.null(saem.pars)) {
      saem.pars <- gsub(
        .regLin,
        "nlmixr_lincmt_\\1 = \\2", saem.pars
      )
    }
    if (!is.null(nlme.mu.fun)) {
      nlme.mu.fun <- gsub(
        .regLin,
        "nlmixr_lincmt_\\1 = \\2", nlme.mu.fun
      )
    }
    pred.txt <- gsub(
      rex::rex("linCmt(", any_spaces, ")"),
      "nlmixr_lincmt_pred", pred.txt
    )
    rest.txt <- gsub(.regLin, "nlmixr_lincmt_\\1 = \\2", rest.txt)
    .tmp <- .bodyDewrap(as.character(attr(fun, "srcref"), useSource = TRUE))
    .tmp <- .tmp[regexpr("~", .tmp) != -1]
    .tmp <- gsub(rex::rex("linCmt(", any_spaces, ")"), "nlmixr_lincmt_pred", .tmp)
    .tmp <- c("function()({", .pred, .rx, ifelse(.hasLinCmt, "nlmixr_lincmt_pred <- linCmt()", ""), "})")
    fun <- eval(parse(text = paste(.tmp, collapse = "\n"), keep.source = TRUE))
    rxode <- TRUE
    .pred <- FALSE
  }
  if (rxode) {
    rx.txt <- .deparse1(body(rest))
    w <- which(regexpr(reg, rx.txt, perl = TRUE) != -1)
    if (length(w) == 0) {
      stop("Error parsing model -- no parameters found.")
    }
    w <- which(regexpr(reg, rx.txt, perl = TRUE) != -1)
    w <- max(w)
    if (any(regexpr(.regRx, rx.txt[1:w], perl = TRUE) != -1)) {
      if (length(rx.txt) == w) {
        tmp <- RxODE::rxGetModel(paste(rx.txt, collapse = "\n"))
        tmp <- paste(tmp$params, "=", tmp$params)
        w <- length(tmp)
        rx.txt <- c(tmp, rx.txt)
      } else {
        stop("Mixed estimation types and ODEs.")
      }
    }
    rx.ode <- rx.txt[-(1:w)]
    rx.pred <- try(eval(parse(text = paste(c("function() {", rx.txt[1:w], "}"), collapse = "\n"))), silent = TRUE)
    if (inherits(rx.pred, "try-error")) {
      .w <- which(regexpr(rex::rex("}"), rx.ode) != -1)
      if (length(.w) > 0) {
        .i <- 1
        .w0 <- w + .w[.i]
        rx.pred <- rx.txt[1:.w0]
        rx.ode <- rx.txt[-(1:.w0)]
        if (regexpr(rex::rex("else"), rx.ode[1]) != -1) {
          .i <- .i + 1
          .w0 <- w + .w[.i]
          rx.pred <- rx.txt[1:.w0]
          rx.ode <- rx.txt[-(1:.w0)]
          w <- .w0
          if (regexpr(rex::rex("else"), rx.ode[1]) != -1) {
            stop("else if are not supported in nlmixr")
          }
        }
        if (any(regexpr(.regRx, rx.pred, perl = TRUE) != -1)) {
          stop("Mixed estimation types and ODEs #2.")
        }
      } else {
        stop("Mixed estimation types and ODEs #3.")
      }
    }
    ## Now separate out parameters for SAEM.
    .tmp <- saem.pars
    nlme.mu.fun2 <- saem.pars
    ## Now separate out for nlme.mu.fun
    rx.ode <- rx.txt[-(1:w)]
    rx.pred <- eval(parse(text = paste(c("function() {", rx.txt[1:w], "}"), collapse = "\n")))
    ## Now separate out parameters for SAEM.
    .no.mu.etas <- c()
    if (!is.null(saem.pars)) {
      w <- max(which(regexpr(reg, saem.pars, perl = TRUE) != -1))
      .tmp <- try(parse(text = paste(c(saem.pars[1:w], "})"), collapse = "\n")), silent = TRUE)
      if (inherits(.tmp, "try-error")) {
        .saemPars2 <- saem.pars[-(1:w)]
        .w <- which(regexpr(rex::rex("}"), .saemPars2) != -1)
        if (length(.w) > 0) {
          .i <- 1
          .w0 <- w + .w[.i]
          .saemPars2 <- saem.pars[-(1:.w0)]
          if (regexpr(rex::rex("else"), .saemPars2[1]) != -1) {
            .i <- .i + 1
            .w0 <- w + .w[.i]
            w <- .w0
          }
        }
      }
      saem.pars <- c(saem.pars[1:w], "")
      .no.mu.etas <- intersect(
        paste(ini$name[!is.na(ini$neta1)]),
        allVars(eval(parse(text = paste0("quote({", paste(saem.pars[-1], collapse = "\n"), "})"))))
      )
      if (length(.no.mu.etas) > 0) {
        saem.pars <- NULL
        nlme.mu.fun2 <- NULL
        nlme.mu.fun <- NULL
      } else {
        ## sapply(saem.pars, message)
        nlme.mu.fun2 <- saem.pars
        w <- max(which(regexpr(reg, nlme.mu.fun, perl = TRUE) != -1))
        .tmp <- try(parse(text = paste(c(nlme.mu.fun[1:w], "})"), collapse = "\n")), silent = TRUE)
        if (inherits(.tmp, "try-error")) {
          .saemPars2 <- nlme.mu.fun[-(1:w)]
          .w <- which(regexpr(rex::rex("}"), .saemPars2) != -1)
          if (length(.w) > 0) {
            .i <- 1
            .w0 <- w + .w[.i]
            .saemPars2 <- nlme.mu.fun[-(1:.w0)]
            if (regexpr(rex::rex("else"), .saemPars2[1]) != -1) {
              .i <- .i + 1
              .w0 <- w + .w[.i]
              w <- .w0
            }
          }
        }
        nlme.mu.fun <- c(nlme.mu.fun[1:w], "")
      }
    }

    rxode <- paste(rx.ode, collapse = "\n")
    if (.linCmt) {
      rxode <- RxODE::rxNorm(RxODE::rxGetLin(rxode))
    }
    rest <- rx.pred
    all.vars <- all.vars[!(all.vars %in% RxODE::rxState(rxode))]
    rest.vars <- rest.vars[!(rest.vars %in% RxODE::rxState(rxode))]
    all.covs <- setdiff(rest.vars, paste0(bounds$name))
    all.covs <- all.covs[!(all.covs %in% RxODE::rxLhs(rxode))]
    all.covs <- all.covs[!(all.covs %in% RxODE::rxState(rxode))]
    all.covs <- setdiff(all.covs, c(
      "t", "time", "podo", "M_E", "M_LOG2E", "M_LOG10E", "M_LN2", "M_LN10", "M_PI", "M_PI_2", "M_PI_4", "M_1_PI",
      "M_2_PI", "M_2_SQRTPI", "M_SQRT2", "M_SQRT1_2", "M_SQRT_3", "M_SQRT_32", "M_LOG10_2", "M_2PI", "M_SQRT_PI",
      "M_1_SQRT_2PI", "M_SQRT_2dPI", "M_LN_SQRT_PI", "M_LN_SQRT_2PI", "M_LN_SQRT_PId2", "pi",
      nlmixrfindLhs(body(rest))
    ))
    .mv <- RxODE::rxModelVars(rxode)
    .state <- c(.mv$state, .mv$stateExtra)
    .tmp <- .data.frame(cmt = seq_along(.state), cond = .state)
    .predDf$dvid <- seq_along(.predDf$var)
    .predDf <- merge(.predDf, .tmp, all.x = TRUE, by = "cond")
    if (any(is.na(.predDf$cmt))) {
      .predDfNa <- .predDf[is.na(.predDf$cmt), names(.predDf) != "cmt"]
      .predDf <- .predDf[!is.na(.predDf$cmt), ]
      .tmp <- .data.frame(cmt = seq_along(.state), var = .state)
      .predDf <- rbind(
        .predDf,
        merge(.predDfNa, .tmp, all.x = TRUE, by = "var")[names(.predDf)]
      )
    }
    .w <- which(is.na(.predDf$cmt))
    .predDf$cond <- paste(.predDf$cond)
    .cmtEndpoints <- integer(0)
    if (length(.w) > 0) {
      .cmtEndpoints <- paste(.predDf$cmt[-.w])
    }
    .predDf$cmt[.w] <- length(.state) + seq_along(.w)
    .w <- .predDf$cond == "" | .predDf$cond != .predDf$var
    .predDf$cond <- paste(.predDf$cond)
    .predDf$cond[.predDf$cond == ""] <- paste(.predDf$var[.predDf$cond == ""])
    .predDf$cmt[.predDf$cond != .predDf$var] <- NA
    .extra <- paste(.predDf$cond[.w])
    .extra <- .extra[regexpr("[()+-]", .extra) == -1]
    if (length(.extra) > 0) {
      rxode <- paste0(rxode, ";\n", paste(c(paste0("cmt(", .extra, ");\n")), collapse = ""))
    }
    if (any(is.na(.predDf$cmt))) {
      .predDfNa <- .predDf[is.na(.predDf$cmt), names(.predDf) != "cmt"]
      .predDf <- .predDf[!is.na(.predDf$cmt), ]
      .mv <- RxODE::rxModelVars(rxode)
      .state <- c(.mv$state, .mv$stateExtra)
      .tmp <- .data.frame(cmt = seq_along(.state), cond = .state)
      .predDf <- rbind(
        .predDf,
        merge(.predDfNa, .tmp, all.x = TRUE, by = "cond")[names(.predDf)]
      )
    }
    if (length(.predDf$cond) > 1) {
      rxode <- paste0(rxode, ";\ndvid(", paste(.predDf[order(.predDf$dvid), "cmt"], collapse = ","), ");\n")
    }
    if (length(.predDf$cmt) > 1L) {
      if (length(.w) > 0L) {
        .nums <- suppressWarnings(as.numeric(paste(.predDf[.w, "cond"])))
        .w2 <- which(!is.na(.nums))
        if (length(.w2) > 0L) {
          .predDf[.w[.w2], "cmt"] <- .nums[.w2]
          .w <- which(is.na(.predDf$cmt))
        }
      }
    }
  } else {
    .predDf$cmt <- -1
    all.covs <- setdiff(rest.vars, paste0(bounds$name))
    nlme.mu.fun2 <- saem.pars
    rxode <- NULL
  }
  fun2 <- c("function(){", .bodyDewrap(as.character(attr(fun, "srcref"), useSource = TRUE)), "}")
  fun2 <- eval(parse(text = paste0(fun2, collapse = "\n"), keep.source = TRUE))
  fun2 <- as.character(attr(fun2, "srcref"), useSource = TRUE)
  fun3 <- c(.bodyDewrap(as.character(attr(.fun000, "srcref"), useSource = TRUE)))
  fun3 <- paste0(fun3, collapse = "\n")

  misplaced.dists <- intersect(rest.funs, c(names(dists), unsupported.dists))
  if (length(misplaced.dists) == 1) {
    if (misplaced.dists == "dt") {
      if (!any(regexpr("[^/]\\bdt[(]", .deparse(rest), perl = TRUE) != -1)) {
        misplaced.dists <- character()
      }
    }
    if (length(misplaced.dists) == 1) {
      if (misplaced.dists == "f") {
        if (!any(regexpr(rex::rex(one_of("Ff"), any_spaces, "(", except_some_of(")\n"), ")", any_spaces, or("<-", "=")), .deparse(rest), perl = TRUE) != -1)) {
          misplaced.dists <- character()
        }
      }
    }
  } else if (length(misplaced.dists) == 2) {
    .tmp <- order(misplaced.dists)
    .tmp <- misplaced.dists[.tmp]
    if (all(misplaced.dists == c("dt", "f"))) {
      if (!any(regexpr("[^/]\\bdt[(]", .deparse(rest), perl = TRUE) != -1)) {
        misplaced.dists <- character()
      }
      if (!any(regexpr(rex::rex(one_of("Ff"), any_spaces, "(", except_some_of(")\n"), ")", any_spaces, or("<-", "=")), .deparse(rest), perl = TRUE) != -1)) {
        misplaced.dists <- character()
      }
    }
  }
  if (length(misplaced.dists) > 0) {
    stop(sprintf("Distributions need to be on residual model lines (like f ~ add(add.err)).\nMisplaced Distribution(s): %s", paste(misplaced.dists, collapse = ", ")))
  }
  tmp <- gsub(rex::rex("linCmt(", any_spaces, ")"), "nlmixr_lincmt_pred", .deparse(pred))
  pred <- eval(parse(text = paste(tmp, collapse = "\n")))
  tmp <- tmp[regexpr(rex::rex("nlmixr_pred <- "), tmp) != -1]
  if (!is.null(saem.pars)) {
    saem.pars <- new.fn(saem.pars)
    saem.theta.trans <- rep(NA, length(theta.names))
  } else {
    saem.pars <- NULL
    saem.theta.trans <- NULL
  }
  if (!is.null(nlme.mu.fun)) {
    nlme.mu.fun <- new.fn(nlme.mu.fun)
  }
  if (!is.null(nlme.mu.fun2)) {
    nlme.mu.fun2 <- new.fn(nlme.mu.fun2)
  }

  cov.theta.pars <- gsub(rex::rex(or(all.covs), "."), "", names(unlist(cov.ref)))
  for (i in seq_along(theta.names)) {
    if (!any(theta.names[i] == cov.theta.pars)) {
      w <- which(theta.names[i] == theta.ord)
      if (length(w) == 1) {
        if (!any(theta.names[i] == cov.theta.pars)) {
          saem.theta.trans[i] <- w
        }
      }
    }
  }
  cur <- 1
  if (!all(is.na(saem.theta.trans))) {
    while (cur <= max(saem.theta.trans, na.rm = TRUE)) {
      while (!any(saem.theta.trans[!is.na(saem.theta.trans)] == cur)) {
        w <- which(saem.theta.trans > cur)
        saem.theta.trans[w] <- saem.theta.trans[w] - 1
      }
      cur <- cur + 1
    }
  }
  env <- new.env(parent = emptyenv())
  env$sum.prod <- FALSE
  ## Split out inPars
  saem.all.covs <- all.covs[all.covs %in% names(cov.ref)]
  saem.inPars <- all.covs[!(all.covs %in% names(cov.ref))]
  .subs <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], quote(`==`)) &&
        all(as.character(x[[3]]) == .what)) {
        if (x[[2]] != "CMT") {
          stop("Multiple endpoints can only be defined in terms of CMT")
        }
        x[[3]] <- eval(parse(text = sprintf("quote(%s)", .with)))
      }
      as.call(lapply(x, .subs))
    } else if (is.pairlist(x)) {
      as.pairlist(lapply(x, .subs))
    } else {
      return(x)
    }
  }
  if (length(.predDf$cond) > 1) {
    for (i in seq_along(.predDf$cond)) {
      .what <- (.predDf$cond[i])
      .with <- paste(.predDf$cmt[i])
      body(pred) <- .subs(body(pred))
      body(err) <- .subs(body(err))
    }
  } else {
    ## Strip if clauses out.
    .strip <- function(x) {
      if (is.call(x)) {
        if (identical(x[[1]], quote(`if`))) {
          return(x[[3]][[2]])
        }
        as.call(lapply(x, .strip))
      } else if (is.pairlist(x)) {
        as.pairlist(lapply(x, .strip))
      } else {
        return(x)
      }
    }
    body(pred) <- .strip(body(pred))
    body(err) <- .strip(body(err))
  }
  .predDf <- .predDf[order(.predDf$cmt), ]
  .w <- which(is.na(.predDf$cond))
  if (length(.w) > 0) .predDf$cond[.w] <- ""
  .predDf$var <- gsub(rex::rex("linCmt(", any_spaces, ")"), "nlmixr_lincmt_pred", paste(.predDf$var))
  .predSaem <- eval(parse(text = sprintf("function(){\n%s;\n}", paste(paste(.predDf$var), collapse = ";\n"))))
  .w <- which(!is.na(bounds$err))
  .tmpErr <- paste(bounds$name[.w])
  .errReg <- rex::rex("nlmixr_", or(.tmpErr), "_par", any_spaces, or("<-", "="), any_spaces, or(.tmpErr))
  .subs <- function(x) {
    if (is.atomic(x)) {
      x
    } else if (is.name(x)) {
      .w <- which(as.character(x) == paste0("nlmixr_", .tmpErr, "_par"))
      if (length(.w) == 1) {
        return(eval(parse(text = sprintf("quote(%s)", .tmpErr[.w]))))
      } else {
        return(x)
      }
    } else if (is.call(x)) {
      as.call(lapply(x, .subs))
    } else if (is.pairlist(x)) {
      as.pairlist(lapply(x, .subs))
    } else {
      return(x)
    }
  }
  ## Now fix the errors again...
  fun2 <- fun2[regexpr(.errReg, fun2) == -1]
  fun2 <- .deparse(.subs(body(eval(parse(text = paste(fun2, collapse = "\n"))))))
  fun2[1] <- paste("function(){")
  body(err) <- .subs(body(err))
  .tmp <- .deparse(body(rest))
  .tmp <- .tmp[regexpr(.errReg, .tmp) == -1]
  .tmp[1] <- "function(){"
  rest <- eval(parse(text = paste(.tmp, collapse = "\n")))
  if (!is.null(saem.pars)) {
    .tmp <- .deparse(body(saem.pars))
    .tmp <- .tmp[regexpr(.errReg, .tmp) == -1]
    .tmp[1] <- "function(){"
    saem.pars <- eval(parse(text = paste(.tmp, collapse = "\n")))
  }
  if (!is.null(nlme.mu.fun)) {
    .tmp <- .deparse(body(nlme.mu.fun))
    .tmp <- .tmp[regexpr(.errReg, .tmp) == -1]
    .tmp[1] <- "function(){"
    nlme.mu.fun <- eval(parse(text = paste(.tmp, collapse = "\n")))
  }
  if (!is.null(nlme.mu.fun2)) {
    .tmp <- .deparse(body(nlme.mu.fun2))
    .tmp <- .tmp[regexpr(.errReg, .tmp) == -1]
    .tmp[1] <- "function(){"
    nlme.mu.fun2 <- eval(parse(text = paste(.tmp, collapse = "\n")))
  }
  .saemErr <- ""
  if (length(intersect(c("t", "time"), allVars(body(rest)))) != 0) {
    .saemErr <- "Initial parameters defined based on ini({}) block variables cannot be defined in terms of time;\nTo fix:\n1. Define ini({}) variables/relationships first.\n2. Use the new variables in an time-based if/else clause later in the `model({})`"
  }
  if (length(.predDf$cond) == 1) {
    .w <- which(bounds$condition == "")
    bounds$condition[.w] <- paste(.predDf$cond)
  }
  ret <- list(
    ini = bounds, model = bigmodel,
    nmodel = list(
      fun = fun2, fun.txt = fun3, pred = pred, error = err, rest = rest, rxode = rxode,
      all.vars = all.vars, rest.vars = rest.vars, all.names = all.names,
      all.funs = all.funs, all.lhs = all.lhs,
      all.covs = all.covs, saem.all.covs = saem.all.covs,
      saem.inPars = saem.inPars,
      errs.specified = errs.specified, add.prop.errs = add.prop.errs,
      grp.fn = grp.fn, mu.ref = .mu.ref, cov.ref = cov.ref, cov.theta = cov.theta,
      saem.pars = saem.pars, nlme.mu.fun = nlme.mu.fun, nlme.mu.fun2 = nlme.mu.fun2,
      log.theta = log.theta,
      log.eta = log.eta, theta.ord = theta.ord, saem.theta.trans = saem.theta.trans,
      predDf = .predDf, predSaem = .predSaem, env = env, predSys = .pred,
      noMuEtas = .no.mu.etas,
      saemErr = .saemErr, cmtEndpoints = .cmtEndpoints,
      oneTheta = .oneTheta, extra = .extra
    )
  )
  if (.linCmt) {
    ret$nmodel$lin.solved <- TRUE
  } else {
    ret$nmodel$lin.solved <- NULL
  }
  return(ret)
}

##' Create the nlme specs list for nlmixr nlme solving
##' @inheritParams nlmixrUI.nlmefun
##' @param mu.type is the mu-referencing type of model hat nlme will be using.
##' @return specs list for nlme
##' @author Matthew L. Fidler
nlmixrUI.nlme.specs <- function(object, mu.type = c("thetas", "covariates", "none")) {
  mu.type <- match.arg(mu.type)
  if (mu.type == "thetas") {
    return(list(
      fixed = object$fixed.form,
      random = object$random.mu,
      start = object$theta
    ))
  } else if (mu.type == "covariates") {
    theta <- names(object$theta)
    cov.ref <- object$cov.ref
    cov.theta <- unique(as.vector(unlist(cov.ref)))
    cov.base <- theta[!(theta %in% cov.theta)]
    cov.base <- cov.base[!(cov.base %in% unlist(lapply(names(cov.ref), function(x) {
      names(cov.ref[[x]])
    })))]
    cov.lst <- list()
    new.theta <- cov.base
    for (n in names(cov.ref)) {
      cov.base <- cov.base[!(cov.base %in% (names(cov.ref[[n]])))]
      cur <- cov.ref[[n]]
      for (i in seq_along(cur)) {
        m <- cur[i]
        cov.lst[[m]] <- c(cov.lst[[m]], n)
        new.theta <- c(new.theta, as.vector(m), names(m))
      }
    }
    e1 <- paste(paste(cov.base, collapse = "+"), "~ 1")
    fixed.form <- paste(c(e1, sapply(names(cov.lst), function(x) {
      paste(x, "~", paste(cov.lst[[x]], collapse = "+"))
    })), collapse = ", ")
    fixed.form <- eval(parse(text = sprintf("list(%s)", fixed.form)))
    if (length(cov.base) == 0) {
      fixed.form <- fixed.form[-1]
    }
    theta <- theta[new.theta]
    return(list(
      fixed = fixed.form,
      random = object$random.mu,
      start = object$theta
    ))
  } else {
    return(list(
      fixed = object$fixed.form,
      random = object$random,
      start = object$theta
    ))
  }
}
##' Create the nlme parameter transform function from the UI object.
##'
##' @param object UI object
##' @param mu Is the model mu referenced?
##' \itemize{
##'
##' \item With the "thetas" only the population parameters are
##' mu-referenced; All covariates are included in the model parameter
##' function.  The between subject variability pieces are specified in
##' the \code{random} specs parameter.
##'
##' \item With the "covariates" option, the population parameters are
##' mu referenced and covariates are removed from the model function.
##' The covariates will be specified used in the fixed effects
##' parameterization of nlme, like \code{list(lKA+lCL~1, lV~WT)}
##'
##' \item With the "none" option, the model function is given to nlme
##' without any modification.
##'
##' }
##' @return Parameter function for nlme
##' @author Matthew L. Fidler
##' @keywords internal
nlmixrUI.nlmefun <- function(object, mu.type = c("thetas", "covariates", "none")) {
  ## create nlme function
  mu.type <- match.arg(mu.type)
  if (mu.type == "thetas") {
    if (length(object$mu.ref) == 0L) {
      return(NULL)
    }
    .all <- object$ini$name[which(object$ini$neta1 == object$ini$neta2)] %in% names(object$mu.ref)
    if (!all(.all)) {
      return(NULL)
    }
    fn <- eval(parse(text = sprintf("function(%s) NULL", paste(unique(c(names(object$ini$theta), object$all.covs)), collapse = ", "))))
    body(fn) <- body(object$nlme.mu.fun)
  } else if (mu.type == "covariates") {
    if (length(object$mu.ref) == 0L) {
      return(NULL)
    }
    vars <- unique(c(unlist(object$mu.ref), unlist(object$cov.ref)))
    fn <- eval(parse(text = sprintf("function(%s) NULL", paste(vars, collapse = ", "))))
    body(fn) <- body(object$nlme.mu.fun2)
    vars2 <- allVars(body(fn))
    if (length(vars) != length(vars2)) {
      return(NULL)
    }
  } else {
    fn <- eval(parse(text = sprintf("function(%s) NULL", paste(object$rest.vars, collapse = ", "))))
    body(fn) <- body(object$rest)
  }
  return(fn)
}
##' Return dynmodel variable translation function
##'
##' @param object nlmixr ui object
##' @return nlmixr dynmodel translation
##' @author Matthew Fidler
nlmixrUI.dynmodelfun <- function(object) {
  .fn <- nlmixrUI.nlmefun(object, "none")
  .fn <- deparse(body(.fn))
  .fn[1] <- paste0("{\n.env <-environment();\nsapply(names(..par),function(x){assign(x,setNames(..par[x],NULL),envir=.env)})\n")
  .fn[length(.fn)] <- paste("return(unlist(as.list(environment())))}")
  .fn <- eval(parse(text = paste0("function(..par)", paste(.fn, collapse = "\n"))))
  return(.fn)
}

##' Return dynmodel variable translation function
##'
##' @param object nlmixr ui object
##' @return nlmixr dynmodel translation
##' @author Matthew Fidler
nlmixrUI.dynmodelfun2 <- function(object) {
  .fn <- nlmixrUI.nlmefun(object, "none")
  .bfn <- body(.fn)
  .fn <- deparse(.bfn)
  .extra <- paste(deparse(nlmixrfindLhs(.bfn)), collapse = " ")
  .fn[1] <- paste0("{\n.env <-environment();\nsapply(names(..par),function(x){assign(x,setNames(..par[[x]],NULL),envir=.env)})\n")
  .fn[length(.fn)] <- paste0(
    ".names <- unique(c(names(..par),",
    .extra, "));\nreturn(as.data.frame(setNames(lapply(.names,function(x){get(x,setNames(..par[x],NULL), envir=.env)}),.names)));\n}"
  )
  .fn <- eval(parse(text = paste0("function(..par)", paste(.fn, collapse = "\n"))))
  return(.fn)
}


##' Get the variance for the nlme fit process based on UI
##'
##' @param object UI object
##' @return nlme/lme variance object
##' @author Matthew L. Fidler
##' @keywords internal
nlmixrUI.nlme.var <- function(object) {
  ## Get the variance for the nlme object
  add.prop.errs <- object$add.prop.errs
  w.no.add <- which(!add.prop.errs$add)
  w.no.prop <- which(!add.prop.errs$prop)
  const <- grp <- ""
  power <- ", fixed=c(1)"
  powera <- ", fixed=list(power=1)"
  if (length(add.prop.errs$y) > 1) {
    grp <- " | nlmixr.grp"
  }
  if (length(w.no.add) > 0) {
    const <- sprintf(", fixed=list(%s)", paste(paste0(add.prop.errs$y[w.no.add], "=0"), collapse = ", "))
  }
  if (length(w.no.prop) > 0) {
    power <- sprintf(", fixed=list(%s)", paste(paste0(add.prop.errs$y, "=", ifelse(add.prop.errs$prop, 1, 0)), collapse = ", "))
    powera <- sprintf(", fixed=list(power=list(%s))", paste(paste0(add.prop.errs$y, "=", ifelse(add.prop.errs$prop, 1, 0)), collapse = ", "))
  }
  tmp <- sprintf("varConstPower(form=~fitted(.)%s%s)", grp, powera)
  if (all(!add.prop.errs$prop)) {
    tmp <- sprintf("varIdent(form = ~ 1%s)", grp)
    if (tmp == "varIdent(form = ~ 1)") {
      warning("Initial condition for additive error ignored with nlme")
      return(NULL)
    }
  } else if (all(!add.prop.errs$add)) {
    tmp <- sprintf("varPower(form = ~ fitted(.)%s%s)", grp, power)
  }
  return(eval(parse(text = tmp)))
}
##' Return RxODE model with predictions appended
##'
##' @param object UI object
##' @return String or NULL if RxODE is not specified by UI.
##' @author Matthew L. Fidler
nlmixrUI.rxode.pred <- function(object) {
  if (is.null(object$rxode)) {
    return(NULL)
  } else {
    tmp <- .deparse1(body(object$pred))
    return(paste(c(object$rxode, tmp), collapse = "\n"))
  }
}
##' Get the Parameter  function with THETA/ETAs defined
##'
##' @param obj UI object
##' @return parameters function defined in THETA[#] and ETA[#]s.
##' @author Matthew L. Fidler
nlmixrUI.theta.pars <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .dft <- .df[!is.na(.df$ntheta), ]
  .unfixed <- with(.dft, sprintf("%s=THETA[%d]", name, seq_along(.dft$name)))
  .eta <- .df[!is.na(.df$neta1), ]
  .eta <- .eta[.eta$neta1 == .eta$neta2, ]
  .eta <- with(.eta, sprintf("%s=ETA[%d]", name, .eta$neta1))
  .f <- .deparse1(body(obj$rest))
  .f <- eval(parse(text = paste(c("function(){", .unfixed, .eta, .f, "}"), collapse = "\n")))
  return(.f)
}
##' Get SAEM distribution
##'
##' @param obj UI object
##' @return Character of distribution
##' @author Matthew L. Fidler
nlmixrUI.saem.distribution <- function(obj) {
  .df <- obj$ini$err
  .df <- paste(.df[which(!is.na(.df))])
  if (any(.df %in% c("dpois", "pois"))) {
    return("poisson")
  }
  if (any(.df %in% c("dbern", "bern", "dbinom", "binom"))) {
    if (.df %in% c("dbinom", "binom")) {
      .df <- obj$ini
      .w <- which(.df$err %in% c("dbinom", "binom"))
      if (length(.w) != 1L) stop("Distribution unsupported by SAEM")
      if (!is.na(.df$name[.w])) stop("Distribution unsupported by SAEM")
    }
    return("binomial")
  }
  if (any(.df %in% c("dnorm", "norm", "prop", "add"))) {
    return("normal")
  }
  if (any(.df %in% c("dlnorm", "lnorm", "logn", "dlogn"))) {
    return("lnorm")
  }
  stop("Distribution unsupported by SAEM")
}
##' Get parameters that are fixed
##'
##' @param obj UI object
##' @return logical vector of fixed THETA parameters
##' @author Matthew L. Fidler
nlmixrUI.focei.fixed <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .dft <- .df[!is.na(.df$ntheta), ]
  .fix <- .dft$fix
  .dft <- .df[is.na(.df$ntheta), ]
  .fix <- c(.fix, .dft$fix)
  return(.fix)
}
##' Get parameters that are fixed for SAEM
##'
##' @param obj UI object
##' @return List of parameters that are fixed.
##' @author Matthew L. Fidler
nlmixrUI.saem.fixed <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .dft <- .df[!is.na(.df$ntheta), ]
  .fixError <- .dft[!is.na(.dft$err), ]
  if (any(.fixError$fix)) {
    stop("Residuals cannot be fixed in SAEM.")
  }
  .dft <- .dft[is.na(.dft$err), ]
  .dft <- setNames(.dft$fix, paste(.dft$name))
  .dft <- .dft[obj$saem.theta.name]
  return(setNames(which(.dft), NULL))
}

##' Get the FOCEi initializations
##'
##' @param obj UI object
##' @return list with FOCEi style initializations
##' @author Matthew L. Fidler
nlmixrUI.focei.inits <- function(obj) {
  df <- .as.data.frame(obj$ini)
  dft <- df[!is.na(df$ntheta), ]
  eta <- df[!is.na(df$neta1), ]
  len <- length(eta$name)
  cur.lhs <- character()
  cur.rhs <- numeric()
  ome <- character()
  for (i in seq_along(eta$name)) {
    last.block <- FALSE
    if (i == len) {
      last.block <- TRUE
    } else if (eta$neta1[i + 1] == eta$neta2[i + 1]) {
      last.block <- TRUE
    }
    if (eta$neta1[i] == eta$neta2[i]) {
      cur.lhs <- c(cur.lhs, sprintf("ETA[%d]", eta$neta1[i]))
      cur.rhs <- c(cur.rhs, eta$est[i])
      if (last.block) {
        ome[length(ome) + 1] <- sprintf(
          "%s ~ %s", paste(cur.lhs, collapse = " + "),
          paste(.deparse(cur.rhs), collapse = " ")
        )
        cur.lhs <- character()
        cur.rhs <- numeric()
      }
    } else {
      cur.rhs <- c(cur.rhs, eta$est[i])
    }
  }
  ome <- eval(parse(text = sprintf("list(%s)", paste(ome, collapse = ","))))
  return(list(
    THTA = dft$est,
    OMGA = ome
  ))
}
##' Get the eta->eta.trans for SAEM
##'
##' @param obj ui object
##' @return list of eta to eta.trans
##' @author Matthew L. Fidler
nlmixrUI.saem.eta.trans <- function(obj) {
  eta.names <- obj$eta.names
  theta.names <- obj$theta.names
  theta.trans <- obj$saem.theta.trans
  mu.ref <- obj$mu.ref
  trans <- rep(NA, length(eta.names))
  for (i in seq_along(eta.names)) {
    ref <- mu.ref[[eta.names[i]]]
    if (!is.null(ref)) {
      w <- which(ref == theta.names)
      if (length(w) == 1) {
        trans[i] <- theta.trans[w]
      }
    }
  }
  if (any(is.na(trans))) {
    stop("Could not figure out the mu-referencing for this model.")
  }
  return(trans)
}
##' Get the SAEM model Omega
##'
##' @param obj UI model
##' @return SAEM model$omega spec
##' @author Matthew L. Fidler
nlmixrUI.saem.model.omega <- function(obj) {
  dm <- sum(!is.na(obj$saem.theta.trans))
  et <- obj$saem.eta.trans
  mat <- matrix(rep(0, dm * dm), dm)
  etd <- which(!is.na(obj$neta1))
  for (i in etd) {
    mat[et[obj$neta1[i]], et[obj$neta2[i]]] <- mat[et[obj$neta2[i]], et[obj$neta1[i]]] <- 1
  }
  return(mat)
}
##' Get the SAEM model$res.mod code
##'
##' @param obj UI model
##' @return SAEM model$res.mod spec
##' @author Matthew L. Fidler
nlmixrUI.saem.res.mod <- function(obj) {
  if (any(obj$saem.distribution == c("poisson", "binomial"))) {
    return(1)
  }
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .hasAdd <- any(.tmp$err == "add") | any(.tmp$err == "norm") | any(.tmp$err == "dnorm") |
      any(.tmp$err == "dlnorm") | any(.tmp$err == "lnorm") | any(.tmp$err == "logn") |
      any(.tmp$err == "dlogn")
    .hasProp <- any(.tmp$err == "prop")
    if (.hasAdd & .hasProp) {
      return(3)
    }
    if (.hasAdd) {
      return(1)
    }
    return(2)
  }))
}
##' Get error names for SAEM
##'
##' @param obj SAEM user interface function.
##' @return Names of error estimates for SAEM
##' @author Matthew L. Fidler
nlmixrUI.saem.res.name <- function(obj) {
  w <- which(sapply(obj$err, function(x) any(x == c("add", "norm", "dnorm", "dlnorm", "lnorm", "logn", "dlogn"))))
  ret <- c()
  if (length(w) == 1) {
    ret[length(ret) + 1] <- paste(obj$name[w])
  }
  w <- which(obj$err == "prop")
  if (length(w) == 1) {
    ret[length(ret) + 1] <- paste(obj$name[w])
  }
  return(ret)
}

##' Get initial estimate for ares SAEM.
##'
##' @param obj UI model
##' @return SAEM model$ares spec
##' @author Matthew L. Fidler
nlmixrUI.saem.ares <- function(obj) {
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(sapply(.tmp$err, function(x) {
      any(x == c(
        "add", "norm", "dnorm", "dpois",
        "pois", "dbinom", "binom", "dbern", "bern",
        "lnorm", "dlnorm", "logn", "dlogn"
      ))
    }))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      return(10)
    }
  }))
}

##' Get initial estimate for bres SAEM.
##'
##' @param obj UI model
##' @return SAEM model$ares spec
##' @author Matthew L. Fidler
nlmixrUI.saem.bres <- function(obj) {
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(sapply(.tmp$err, function(x) any(x == "prop")))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      return(1)
    }
  }))
}

##' Get model$log.eta for SAEM
##'
##' @param obj UI model
##' @return SAEM model$log.eta
##' @author Matthew L. Fidler
nlmixrUI.saem.log.eta <- function(obj) {
  lt <- obj$log.theta
  dm <- sum(!is.na(obj$saem.theta.trans))
  ret <- rep(FALSE, dm)
  theta.trans <- obj$saem.theta.trans
  theta.names <- obj$theta.names
  for (n in lt) {
    w <- which(n == theta.names)
    if (length(w) == 1) {
      ret[theta.trans[w]] <- TRUE
    }
  }
  return(ret)
}

##' Generate saem.fit user function.
##'
##' @param obj UI object
##' @return saem user function
##' @author Matthew L. Fidler
nlmixrUI.saem.fit <- function(obj) {
  if (any(ls(envir = obj$env) == "saem.fit")) {
    return(obj$env$saem.fit)
  } else if (!is.null(obj$rxode)) {
    ## RxODE function
    if (is.null(obj$saem.pars)) {
      stop("SAEM requires mu-referenced parameters")
    }
    message("Compiling RxODE equations...", appendLF = FALSE)
    if (obj$env$sum.prod) {
      ode <- RxODE::RxODE(RxODE::rxSumProdModel(obj$rxode))
    } else {
      ode <- RxODE::RxODE(obj$rxode)
    }
    RxODE::rxLoad(ode)
    obj$env$saem.ode <- ode
    RxODE::rxLoad(ode)
    message("done.")
    inPars <- obj$saem.inPars
    if (length(inPars) == 0) {
      inPars <- NULL
    } else {
      ## Check for inPars in Covariates in RxODE model
      .extra <- RxODE::rxModelVars(obj$rxode)
      .extra <- intersect(obj$saem.all.covs, .extra$params)
      if (length(.extra) > 0) {
        inPars <- c(inPars, .extra)
      }
    }
    saem.fit <- gen_saem_user_fn(model = ode, obj$saem.pars, pred = obj$predSaem, inPars = inPars)
    obj$env$saem.fit <- saem.fit
    return(obj$env$saem.fit)
  }
}
##' Generate SAEM model list
##'
##' @param obj  nlmixr UI object
##' @return SAEM model list
##' @author Matthew L. Fidler
nlmixrUI.saem.model <- function(obj) {
  mod <- list(saem_mod = obj$saem.fit)
  if (length(obj$saem.all.covs > 0)) {
    mod$covars <- obj$saem.all.covs
  }
  mod$res.mod <- obj$saem.res.mod
  mod$log.eta <- obj$saem.log.eta
  ## if (FALSE){
  ## FIXME option/warning
  mod$ares <- obj$saem.ares
  mod$bres <- obj$saem.bres
  ## }
  mod$omega <- obj$saem.model.omega
  return(mod)
}
##' Get THETA names for nlmixr's SAEM
##'
##' @param uif nlmixr UI object
##' @return SAEM theta names
##' @author Matthew L. Fidler
nlmixrUI.saem.theta.name <- function(uif) {
  .trans <- uif$saem.theta.trans
  .df <- .as.data.frame(uif$ini)
  .df <- .df[!is.na(.df$ntheta), ]
  .transName <- paste(.df$name[which(!is.na(.trans))])
  .trans <- .trans[!is.na(.trans)]
  theta.name <- .transName[order(.trans)]
  all.covs <- uif$saem.all.covs
  lc <- length(all.covs)
  if (lc > 0) {
    m <- matrix(rep(NA, length(theta.name) * (lc + 1)), nrow = lc + 1)
    dimnames(m) <- list(c("_name", all.covs), theta.name)
    m["_name", ] <- theta.name
    for (cn in names(uif$cov.ref)) {
      v <- uif$cov.ref[[cn]]
      for (var in names(v)) {
        rn <- v[var]
        m[cn, rn] <- var
      }
    }
    ret <- unlist(m)
    ret <- ret[!is.na(ret)]
    return(ret)
  }
  return(theta.name)
}
##' Generate SAEM initial estimates for THETA.
##'
##' @param obj nlmixr UI object
##' @return SAEM theta initial estimates
##' @author Matthew L. Fidler
nlmixrUI.saem.init.theta <- function(obj) {
  theta.name <- obj$saem.theta.name
  cov.names <- unique(names(unlist(structure(obj$cov.ref, .Names = NULL))))
  theta.name <- theta.name[!(theta.name %in% cov.names)]
  nm <- paste(obj$ini$name)
  lt <- obj$log.theta
  i <- 0
  this.env <- environment()
  theta.ini <- sapply(theta.name, function(x) {
    w <- which(x == nm)
    assign("i", i + 1, this.env)
    if (any(lt == x)) {
      return(exp(obj$ini$est[w]))
    } else {
      return(obj$ini$est[w])
    }
  })
  all.covs <- obj$saem.all.covs
  lc <- length(all.covs)
  if (lc > 0) {
    m <- matrix(rep(NA, lc * length(theta.name)), ncol = lc)
    dimnames(m) <- list(theta.name, all.covs)
    for (cn in names(obj$cov.ref)) {
      v <- obj$cov.ref[[cn]]
      for (var in names(v)) {
        rn <- v[var]
        w <- which(var == nm)
        m[rn, cn] <- obj$ini$est[w]
      }
    }
    return(as.vector(c(theta.ini, as.vector(m))))
  }
  return(as.vector(theta.ini))
}
##' SAEM's init$omega
##'
##' @param obj nlmixr UI object
##' @param names When \code{TRUE} return the omega names.  By default
##'     this is \code{FALSE}.
##' @return Return initial matrix
##' @author Matthew L. Fidler
nlmixrUI.saem.init.omega <- function(obj, names = FALSE) {
  dm <- sum(!is.na(obj$saem.theta.trans))
  et <- obj$saem.eta.trans
  ret <- rep(NA, dm)
  etd <- which(obj$neta1 == obj$neta2)
  for (i in etd) {
    if (names) {
      ret[et[obj$neta1[i]]] <- paste(obj$name[i])
    } else {
      ret[et[obj$neta1[i]]] <- obj$est[i]
    }
  }
  if (names) {
    ret <- ret[!is.na(ret)]
    return(ret)
  } else {
    tmp <- unique(ret[!is.na(ret)])
    if (length(tmp) == 1) {
      ret[is.na(ret)] <- tmp
    } else {
      ret[is.na(ret)] <- 1
    }
  }
  return(ret)
}
##' Get saem initilization list
##'
##' @param obj nlmixr UI object
##' @return Return SAEM inits list.
##' @author Matthew L. Fidler
nlmixrUI.saem.init <- function(obj) {
  ret <- list()
  ret$theta <- obj$saem.init.theta
  ## if (FALSE){
  ret$omega <- obj$saem.init.omega
  ## }
  return(ret)
}

nlmixrUI.focei.mu.ref <- function(obj) {
  .muRef <- obj$mu.ref
  .tn <- obj$focei.names
  sapply(seq_along(.muRef), function(x) {
    .cur <- .muRef[[x]]
    .w <- which(.tn == .cur)
    if (length(.w) == 1) {
      return(.w - 1)
    } else {
      return(-1)
    }
  })
}

nlmixrUI.model.desc <- function(obj) {
  .mv <- RxODE::rxModelVars(obj$rxode.pred)
  if (obj$predSys) {
    return("RxODE-based Pred model")
  } else if (.mv$extraCmt == 0) {
    return("RxODE-based ODE model")
  } else {
    .pred <- obj$rxode.pred
    if (regexpr(rex::rex("rx_B~0;"), .pred)) {
      .ncmt <- 1
    } else if (regexpr(rex::rex("rx_C~0;"), .pred)) {
      .ncmt <- 2
    } else {
      .ncmt <- 3
    }
    return(sprintf(
      "RxODE-based %s-compartment model%s", .ncmt,
      ifelse(.mv$extraCmt == 2, " with first-order absorption", "")
    ))
  }
}

nlmixrUI.lincmt.dvdx <- function(obj) {
  if (is.null(obj$rxode.pred)) {
    .df <- .as.data.frame(obj$ini)
    .dft <- .df[!is.na(.df$ntheta), ]
    .unfixed <- with(.dft, sprintf("%s=THETA[%d]", name, seq_along(.dft$name)))
    .eta <- .df[!is.na(.df$neta1), ]
    .eta <- .eta[.eta$neta1 == .eta$neta2, ]
    .eta <- with(.eta, sprintf("%s=ETA[%d]", name, .eta$neta1))
    .txt <- .deparse1(body(obj$rest))
    .txt[length(.txt) + 1] <- obj$lin.solved$extra.lines
    .txt <- paste(c(.unfixed, .eta, .txt), collapse = "\n")
    .txt <- substring(.txt, 0, nchar(.txt) - 1)
    return(RxODE::rxSymPyLincmtDvdx(.txt, obj$lin.solved$ncmt, obj$lin.solved$parameterization))
  }
  return(NULL)
}

nlmixrUI.poped.notfixed_bpop <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .tmp <- .df[!is.na(.df$ntheta) & is.na(.df$err), ]
  return(setNames(1 - .tmp$fix * 1, paste(.tmp$name)))
}

nlmixrUI.poped.d <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .tmp <- .df[which(is.na(.df$ntheta) & .df$neta1 == .df$neta2), ]
  return(setNames(.tmp$est, paste(.tmp$name)))
}

nlmixrUI.poped.sigma <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .tmp <- .df[!which(is.na(.df$err) & .df$neta1 == .df$neta2), ]
  return(setNames(.tmp$est * .tmp$est, paste(.tmp$name)))
}

nlmixUI.logThetasList <- function(obj) {
  .ini <- .as.data.frame(obj$ini)
  .logThetas <- as.integer(which(setNames(sapply(obj$focei.names, function(x) any(x == obj$log.theta)), NULL)))
  .thetas <- .ini[!is.na(.ini$ntheta), ]
  .one <- obj$oneTheta
  .logThetasF <- .thetas[.thetas$name %in% .one, "ntheta"]
  .logThetasF <- intersect(.logThetas, .logThetasF)
  list(.logThetas, .logThetasF)
}


nlmixrUI.poped.ff_fun <- function(obj) {
  if (!is.null(obj$lin.solved)) {
    stop("Solved system not supported yet.")
  } else {
    .df <- .as.data.frame(obj$ini)
    .dft <- .df[!is.na(.df$ntheta) & is.na(.df$err), ]
    .unfixed <- with(.dft, sprintf("%s=bpop[%d]", name, seq_along(.dft$name)))
    .eta <- .df[!is.na(.df$neta1), ]
    .eta <- .eta[.eta$neta1 == .eta$neta2, ]
    .eta <- with(.eta, sprintf("%s=b[%d]", name, .eta$neta1))
    .lhs <- nlmixrfindLhs(body(obj$rest))
    .f <- .deparse(body(obj$rest))[-1]
    .lhs <- sprintf("return(c(%s))", paste(sprintf("\"%s\"=%s", .lhs, .lhs), collapse = ", "))
    .f <- eval(parse(text = paste(c("function(x,a,bpop,b,bocc){", .unfixed, .eta, .f[-length(.f)], .lhs, "}"), collapse = "\n")))
    return(.f)
  }
}

##' @export
`$.nlmixrUI` <- function(obj, arg, exact = TRUE) {
  x <- obj
  .cls <- class(x)
  class(x) <- "list"
  if (arg == "ini") {
    return(x$ini)
  } else if (arg == "nmodel") {
    return(x$nmodel)
  } else if (arg == "model") {
    return(x$model)
  } else if (arg == "nlme.fun.mu") {
    return(nlmixrUI.nlmefun(obj, "thetas"))
  } else if (arg == "dynmodel.fun") {
    return(nlmixrUI.dynmodelfun(obj))
  } else if (arg == "dynmodel.fun.df") {
    return(nlmixrUI.dynmodelfun2(obj))
  } else if (arg == "nlme.fun") {
    return(nlmixrUI.nlmefun(obj, "none"))
  } else if (arg == "nlme.fun.mu.cov") {
    return(nlmixrUI.nlmefun(obj, "covariates"))
  } else if (arg == "nlme.specs") {
    return(nlmixrUI.nlme.specs(obj, "none"))
  } else if (arg == "nlme.specs.mu") {
    return(nlmixrUI.nlme.specs(obj, "thetas"))
  } else if (arg == "nlme.specs.mu.cov") {
    return(nlmixrUI.nlme.specs(obj, "covariates"))
  } else if (arg == "nlme.var") {
    return(nlmixrUI.nlme.var(obj))
  } else if (arg == "rxode.pred") {
    return(nlmixrUI.rxode.pred(obj))
  } else if (arg == "theta.pars") {
    return(nlmixrUI.theta.pars(obj))
  } else if (arg == "focei.inits") {
    return(nlmixrUI.focei.inits(obj))
  } else if (arg == "focei.fixed") {
    return(nlmixrUI.focei.fixed(obj))
  } else if (arg == "focei.mu.ref") {
    return(nlmixrUI.focei.mu.ref(obj))
  } else if (arg == "saem.fixed") {
    return(nlmixrUI.saem.fixed(obj))
  } else if (arg == "saem.theta.name") {
    return(nlmixrUI.saem.theta.name(obj))
  } else if (arg == "saem.eta.trans") {
    return(nlmixrUI.saem.eta.trans(obj))
  } else if (arg == "saem.model.omega") {
    return(nlmixrUI.saem.model.omega(obj))
  } else if (arg == "saem.res.mod") {
    return(nlmixrUI.saem.res.mod(obj))
  } else if (arg == "saem.ares") {
    return(nlmixrUI.saem.ares(obj))
  } else if (arg == "saem.bres") {
    return(nlmixrUI.saem.bres(obj))
  } else if (arg == "saem.log.eta") {
    return(nlmixrUI.saem.log.eta(obj))
  } else if (arg == "saem.fit") {
    return(nlmixrUI.saem.fit(obj))
  } else if (arg == "saem.model") {
    return(nlmixrUI.saem.model(obj))
  } else if (arg == "saem.init.theta") {
    return(nlmixrUI.saem.init.theta(obj))
  } else if (arg == "saem.init.omega") {
    return(nlmixrUI.saem.init.omega(obj))
  } else if (arg == "saem.init") {
    return(nlmixrUI.saem.init(obj))
  } else if (arg == "saem.omega.name") {
    return(nlmixrUI.saem.init.omega(obj, TRUE))
  } else if (arg == "saem.res.name") {
    return(nlmixrUI.saem.res.name(obj))
  } else if (arg == "model.desc") {
    return(nlmixrUI.model.desc(obj))
  } else if (arg == "meta") {
    return(x$meta)
  } else if (arg == "saem.distribution") {
    return(nlmixrUI.saem.distribution(obj))
  } else if (arg == "lincmt.dvdx") {
    return(nlmixrUI.lincmt.dvdx(obj))
  } else if (arg == "notfixed_bpop" || arg == "poped.notfixed_bpop") {
    return(nlmixrUI.poped.notfixed_bpop(obj))
  } else if (arg == "poped.ff_fun") {
    return(nlmixrUI.poped.ff_fun(obj))
  } else if (arg == "poped.d") {
    return(nlmixrUI.poped.d(obj))
  } else if (arg == "poped.sigma") {
    return(nlmixrUI.poped.sigma(obj))
  } else if (arg == "logThetasList") {
    return(nlmixUI.logThetasList(obj))
  } else if (arg == ".clean.dll") {
    if (exists(".clean.dll", envir = x$meta)) {
      clean <- x$meta$.clean.dll
      if (is(clean, "logical")) {
        return(clean)
      }
    }
    return(TRUE)
  } else if (arg == "random.mu") {
    return(nlmixrBoundsOmega(x$ini, x$nmodel$mu.ref))
  } else if (arg == "bpop") {
    arg <- "theta"
  } else if (arg == "multipleEndpoint") {
    return(nlmixrUI.multipleEndpoint(x))
  } else if (arg == "muRefTable") {
    class(x) <- .cls
    return(.nmMuTable(x))
  }
  m <- x$ini
  ret <- `$.nlmixrBounds`(m, arg, exact = exact)
  if (is.null(ret)) {
    m <- x$nmodel
    ret <- m[[arg, exact = exact]]
    if (is.null(ret)) {
      if (exists(arg, envir = x$meta)) {
        ret <- get(arg, envir = x$meta)
      }
    }
  }
  ret
}



##' @export
str.nlmixrUI <- function(object, ...) {
  obj <- object
  class(obj) <- "list"
  str(obj$ini)
  str(obj$nmodel)
  cat(" $ ini       : Model initilizations/bounds object\n")
  cat(" $ model     : Original Model\n")
  cat(" $ model.desc: Model description\n")
  cat(" $ nmodel    : Parsed Model List\n")
  cat(" $ nlme.fun  : The nlme model function.\n")
  cat(" $ nlme.specs: The nlme model specs.\n")
  cat(" $ nlme.var  : The nlme model varaince.\n")
  cat(" $ rxode.pred: The RxODE block with pred attached (final pred is nlmixr_pred)\n")
  cat(" $ theta.pars: Parameters in terms of THETA[#] and ETA[#]\n")
  cat(" $ focei.inits: Initialization for FOCEi style blocks\n")
  cat(" $ focei.fixed: Logical vector of FOCEi fixed parameters\n")
  cat(" $ focei.mu.ref: Integer Vector of focei.mu.ref\n")
  cat(" $ saem.eta.trans: UI ETA -> SAEM ETA\n")
  cat(" $ saem.model.omega: model$omega for SAEM\n")
  cat(" $ saem.res.mod: model$res.mod for SAEM\n")
  cat(" $ saem.ares: model$ares for SAEM\n")
  cat(" $ saem.bres: model$bres for SAEM\n")
  cat(" $ saem.log.eta: model$log.eta for SAEM\n")
  cat(" $ saem.fit  : The SAEM fit user function\n")
  cat(" $ saem.model: The SAEM model list\n")
  cat(" $ saem.init.theta: The SAEM init$theta\n")
  cat(" $ saem.init.omega: The SAEM init$omega\n")
  cat(" $ saem.init : The SAEM inits list\n")
  cat(" $ saem.theta.name : The SAEM theta names\n")
  cat(" $ saem.omega.name : The SAEM theta names\n")
  cat(" $ saem.res.name : The SAEM omega names\n")
  cat(" $ saem.distribution: SAEM distribution\n")
  cat(" $ .clean.dll : boolean representing if dlls are cleaned after running.\n")
  cat(" $ logThetasList: List of logThetas:\n")
  cat("     first element are scaling log thetas;\n")
  cat("     second element are back-transformed thetas;\n")
  cat(" $ multipleEndpoint: table/huxtable of multiple endpoint translations in nlmixr\n")
  cat(" $ muRefTable: table/huxtable of mu-referenced items in a model\n")
}


## Local Variables:
## ess-indent-offset: 2
## indent-tabs-mode: nil
## End:
