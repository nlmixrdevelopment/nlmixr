## utils.R: population PK/PD modeling library
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

# Utilities for nlmixr ####################################################

# ####################################################################### #
#
## Utilities for building nlmixr
#
# ####################################################################### #

refresh <- function() {
  ## nocov start
  source(devtools::package_file("build/refresh.R"))
  ## nocov end
}

nsis <- function() { ## build installer...
  ## nocov start
  source(devtools::package_file("build/nsis.R"))
  ## nocov end
}
# ########################################################################

# .collectWarnings --------------------------------------------------------
##' Collect warnings and just warn once.
##'
##' @param expr R expression
##' @param lst When \code{TRUE} return a list with
##'     list(object,warnings) instead of issuing the warnings.
##'     Otherwise, when \code{FALSE} issue the warnings and return the
##'     object.
##' @return The value of the expression or a list with the value of
##'     the expression and a list of warning messages
##' @author Matthew L. Fidler
##' @noRd
.collectWarnings <- function(expr, lst = FALSE) {
  ws <- c()
  this.env <- environment()
  ret <-
    suppressWarnings(withCallingHandlers(
      expr,
      warning = function(w) {
        assign("ws", unique(c(w$message, ws)), this.env)
      }
    ))
  if (lst) {
    return(list(ret, ws))
  } else {
    for (w in ws) {
      warning(w)
    }
    return(ret)
  }
}
# #########################################################################

# nlmixrPrint() -----------------------------------------------------------
##' Print x using the message facility
##'
##' This allows the suppressMessages to work on print functions.  This
##' captures the output function sends it through the message routine.
##'
##' catpureOutput was used since it is much faster than the internal
##' capture.output see https://www.r-bloggers.com/performance-captureoutput-is-much-faster-than-capture-output/
##' @param x object to print
##' @param ... Other things output
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
nlmixrPrint <- function(x, ...) {
  this.env <- environment()
  message(invisible(paste(
    .captureOutput(assign("x", print(x, ...), this.env)),
    collapse = "\n"
  )),
  appendLF = TRUE
  )
  invisible(x)
}
# #########################################################################

.dontRun <- function(...) {
  ## This is for r checks, though they need to be loaded...
  vpc::vpc(...)
  dparser::dparse(...)
}

# cholSE() ----------------------------------------------------------------
##' Generalized Cholesky Matrix Decomposition
##'
##'  Performs a (modified) Cholesky factorization of the form
##'
##'   t(P) \%*\% A \%*\% P  + E = t(R) \%*\% R
##'
##'  As detailed in Schnabel/Eskow (1990)
##'
##' @param matrix Matrix to be Factorized.
##' @param tol Tolerance; Algorithm suggests (.Machine$double.eps) ^ (1 / 3), default
##' @return Generalized Cholesky decomposed matrix.
##' @author Matthew L. Fidler (translation), Johannes Pfeifer, Robert
##'     B. Schnabel and Elizabeth Eskow
##'
##' @references
##'
##' matlab source: http://www.dynare.org/dynare-matlab-m2html/matlab/chol_SE.html; Slightly different return values
##'
##' Robert B. Schnabel and Elizabeth
##' Eskow. 1990. "A New Modified Cholesky Factorization," SIAM Journal
##' of Scientific Statistical Computing, 11, 6: 1136-58.
##'
##' Elizabeth Eskow and Robert B. Schnabel
##' 1991. "Algorithm 695 - Software for a New Modified Cholesky Factorization,"
##' ACM Transactions on Mathematical Software, Vol 17, No 3: 306-312
##'
##' @note
##'
##' This version does not pivot or return the E matrix
##'
##' @export
cholSE <- function(matrix, tol = (.Machine$double.eps)^(1 / 3)) {
  .Call(`_nlmixr_cholSE_`, matrix, tol)
}
# #########################################################################

.setRoot <- function() {
  setwd("c:/")
}

#' Generate a data.frame using the R4.0 convention
#'
#' @param ... Passed to \code{base::data.frame()} or
#'   \code{base::as.data.frame()}
#' @param stringsAsFactors Captured so that it can be ignored and always set to
#'   \code{FALSE}
#' @return A data.frame with strings not converted to factors
#' @noRd
.data.frame <- function(..., stringsAsFactors = FALSE) {
  base::data.frame(..., stringsAsFactors = FALSE)
}
.as.data.frame <- function(..., stringsAsFactors = FALSE) {
  base::as.data.frame(..., stringsAsFactors = FALSE)
}


.isTestthat <- function() {
  return(regexpr("/tests/testthat/", getwd(), fixed = TRUE) != -1)
}

nlmixrTest <- function(expr, silent = .isTestthat(), test = "cran") {
  .Call(`_nlmixr_setSilentErr`, 1L, PACKAGE = "nlmixr")
  RxODE::rxSetSilentErr(1L)
  do.it <- TRUE
  .test <- .test0 <- Sys.getenv("NOT_CRAN")
  if (Sys.getenv("nmCran") != "") {
    .test <- .test0 <- Sys.getenv("nmCran")
  }
  on.exit({
    .Call(`_nlmixr_setSilentErr`, 0L, PACKAGE = "nlmixr")
    RxODE::rxSetSilentErr(0L)
  })
  if (any(.test == c("false", "", "cran"))) {
    if (any(test == c("false", "", "cran"))) {
      do.it <- TRUE
    }
    else {
      do.it <- FALSE
    }
  }
  else {
    if (test == .test || .test == "true") {
      do.it <- TRUE
    }
    else {
      do.it <- FALSE
    }
  }
  if (do.it) {
    .lastCran <- Sys.getenv("NOT_CRAN")
    Sys.setenv(NOT_CRAN = "true")
    on.exit(
      {
        Sys.setenv(NOT_CRAN = .lastCran)
      },
      add = TRUE
    )
    if (is(substitute(expr), "{")) {
      if (silent) {
        return(suppressMessages(eval(substitute(expr),
          envir = parent.frame(1)
        )))
      }
      else {
        return(eval(substitute(expr), envir = parent.frame(1)))
      }
    }
  }
}
