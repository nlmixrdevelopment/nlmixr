##' Simulated Data Set for comparing objective functions
##'
##' This is a simulated dataset from Wang2007 where various NONMEM
##' estimation methods (Lapclace FO, FOCE with and without
##' interaction) are described.
##'
##' @format A data frame with 20 rows and 3 rows
##'
##' \describe{
##'   \item{ID}{Siumlated Subject ID}
##'   \item{Time}{Simulated Time}
##'   \item{Time}{Simulated Objective Function}
##' }
##'
##' @source Table 1 from Wang, Y \emph{Derivation of Various NONMEM estimation methods}. J Pharmacokinet Pharmacodyn (2007) 34:575-593.
"Wang2007"

##' @templateVar est.method FO
##' @templateVar est.err Additive
##' @template w2007
"Wang2007_add_fo"

##' @templateVar est.method FOCE
##' @templateVar est.err Additive
##' @template w2007
"Wang2007_add_foce"

##' @templateVar est.method FOCE
##' @templateVar est.err Proportional
##' @template w2007
"Wang2007_prop_foce"

##' @templateVar est.method FO
##' @templateVar est.err Proportional
##' @template w2007
"Wang2007_prop_fo"

##' @templateVar est.method FOCE with Interaction
##' @templateVar est.err Proportional
##' @template w2007
"Wang2007_prop_focei"

##' @templateVar est.method FO
##' @templateVar est.err Exponential
##' @template w2007
"Wang2007_exp_fo"

##' @templateVar est.method FOCE
##' @templateVar est.err Exponential
##' @template w2007
"Wang2007_exp_foce"


##' @templateVar est.method FOCE with Interaction
##' @templateVar est.err Exponential
##' @template w2007
"Wang2007_exp_focei"
