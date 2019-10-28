##' Multiple dose Theophilline PK data
##'
##' This datastart is  the day 1 concentrations of the
##' Theophiline data that is included in the nlme/NONMEM.
##'
##' @format A data frame with 144 rows by 7 columns
##'
##' \describe{
##'   \item{ID}{Subject ID}
##'   \item{TIME}{Time (hrs)}
##'   \item{DV}{Dependant Variable, Theophiline Concentration}
##'   \item{AMT}{Dose Amount/kg}
##'   \item{EVID}{RxODE/nlmixr event ID (not NONMEM's)}
##'   \item{CMT}{Compartment Number}
##'   \item{WT}{Weight (kg)}
##' }
##'
##' @source NONMEM/nlme
"theo_sd"
