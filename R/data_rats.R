##' Pregnant Rat Diet Experiment
##'
##' 16 pregnant rats have a control diet, and 16 have a chemically
##' treated diet.  The litter size for each rat is recorded after 4
##' and 21 days.  This dataset is used in the SAS Probit-model with
##' binomial data, and saved in the nlmixr package as rats.
##'
##' @format A data frame with 32 rows and 6 columns
##'
##' \describe{
##'   \item{trt}{Treatment; c= control diet; t=treated diet}
##'   \item{m}{Litter size after 4 days}
##'   \item{x}{Litter size after 21 days}
##'   \item{x1}{Indicator for trt=c}
##'   \item{x2}{Indicator for trt=t}
##'   \item{ID}{Rat ID}
##' }
##'
##' @source \url{https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_nlmixed_sect040.htm}
##'
##' @references Weil, C.S., 1970. Selection of the valid number of
##'     sampling units and a consideration of their combination in
##'     toxicological studies involving reproduction, teratogenesis or
##'     carcinogenesis. Fd. Cosmet. Toxicol. 8, 177-182.
##'
##' @references Williams, D.A., 1975. The analysis of binary responses
##'     from toxicological experiments involving reproduction and
##'     teratogenicity. Biometrics 31, 949-952.
##'
##' @references McCulloch, C. E. (1994),
##'     "Maximum Likelihood Variance Components Estimation for Binary Data,"
##'     Journal of the American Statistical Association, 89, 330 -
##'     335.
##'
##' @references Ochi, Y. and Prentice, R. L. (1984),
##'     "Likelihood Inference in a Correlated Probit Regression Model,"
##'     Biometrika, 71, 531-543.
##'
##' @family nlmixr datasets
"rats"
