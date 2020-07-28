## gauss.quad.R: population PK/PD modeling library
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


#' Sets nodes and weights of Gauss-Hermite quadrature
#'
#' @param n number of nodes
#' @return a list of nodes and weights of Gauss-Hermite quadrature
#' @examples
#' gauss.quad(5)
#' @export
gauss.quad <- function(n) {
  s <- fastGHQuad::gaussHermiteData(n)
  names(s) <- c("nodes", "weights")
  s
}
