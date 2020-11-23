## instant.stan.extension.R: population PK/PD modeling library
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

#' instant.stan.extension.
#'
#' instant.stan.extension
#'
#' @param ode_str ODE equations in a string
#' @param covar a character vector of covariates
#' @return NULL
instant.stan.extension <- function(ode_str = NULL, covar = NULL) {
  if (is.null(ode_str)) {
    stop("please provide ODE string")
  }

  if (is.null(covar)) {
    .tmpl <- system.file("include/generic_ode_interface_template.txt", package = "nlmixr")
    cat(ode_str, file = "model.txt")
    nvar <- 0
  } else {
    .tmpl <- system.file("include/generic_ode_interface_template_cov.txt", package = "nlmixr")
    covar <- strsplit(covar, "[,| \t]+")[[1]]
    nvar <- length(covar)
    if (prod(nchar(covar)) * nvar == 0) {
      stop("unrecoganized covar string")
    }
    paste(covar, "= 9999.999 + -9999.999;")
    cat(paste(covar, "= 9999.999 + -9999.999;"), file = "model.txt", sep = "\n")
    cat(ode_str, file = "model.txt", append = TRUE)
  }
  .extn <- system.file("include/stan/math/PMXStan", package = "StanHeaders")
  x <- .C("parse_ode", .tmpl, "model.txt", sprintf("%s/generic_ode_interface.hpp", .extn), as.character(nvar))

  pars <- scan("ODE_PARS.txt", what = "", quiet = TRUE)
  cat("A new ODE extension for Stan has been created.\n")
  cat(sprintf("System parameters are: %s\n", paste(pars, collapse = " ")))
  invisible()
}
