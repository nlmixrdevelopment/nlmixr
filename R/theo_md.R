##' Multiple dose Theophiline PK data
##'
##' This datastart starts with the day 1 concentrations of the
##' Theophiline data that is included in the nlme/NONMEM. After day 7
##' concentrations were simulated with once a day regimen for 7 days
##' (QD).
##'
##' @format A data frame with 348 rows by 7 columns
##'
##' \describe{
##'   \item{ID}{Subject ID}
##'   \item{TIME}{Time (hrs)}
##'   \item{DV}{Dependant Variable, Theophiline Concentration}
##'   \item{AMT}{Dose Amount}
##'   \item{EVID}{RxODE/nlmixr event ID (not NONMEM's)}
##'   \item{CMT}{Compartment number}
##'   \item{WT}{Weight (kg)}
##' }
##'
##' @source NONMEM/nlme
"theo_md"
