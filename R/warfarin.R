##' Warfarin PK/PD data
##'
##' @format A data frame with 519 rows and 9 columns
##'
##' \describe{
##'   \item{id}{Patient identifier (n=32)}
##'   \item{time}{Time [h]}
##'   \item{amt}{Total drug administered [mg]}
##'   \item{dv}{Warfarin concentrations [mg/L] or PCA measurement}
##'   \item{dvid}{Dependent identifier Information (cp: Dose or PK, pca: PCA, factor)}
##'   \item{evid}{Event identifier}
##'   \item{wt}{Weight [kg]}
##'   \item{age}{Age [yr]}
##'   \item{sex}{Gender (male or female, factor)}
##' }
##'
##' @source Funaki T, Holford N, Fujita S (2018). Population PKPD
##'     analysis using nlmixr and NONMEM. PAGJA 2018
##'
##' @references
##'
##' O'Reilly RA, Aggeler PM, Leong LS. Studies of the coumarin
##'   anticoagulant drugs: The pharmacodynamics of warfarin in man.
##'   Journal of Clinical Investigation 1963; 42(10): 1542-1551
##'
##' O'Reilly RA, Aggeler PM. Studies on coumarin anticoagulant drugs
##'   Initiation of warfarin therapy without a loading dose. Circulation
##'   1968; 38: 169-177.
##'
"warfarin"
