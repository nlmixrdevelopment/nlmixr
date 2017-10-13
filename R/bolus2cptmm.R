##' 2 Compartment Model with Michelis Menton Clearance
##'
##' This is a simulated dataset from the ACOP 2016 poster.  All
##' Datasets were simulated with the following methods.
##'
##' Richly sampled profiles were simulated for 4 different dose levels
##' (10, 30, 60 and 120 mg) of 30 subjects each as single dose (over
##' 72h), multiple dose (4 daily doses), single and multiple dose
##' combined, and steady state dosing, for a range of test models: 1-
##' and 2-compartment disposition, with and without 1st order
##' absorption, with either linear or Michaelis-Menten (MM)
##' clearance(MM without steady state dosing). This provided a total
##' of 42 test cases. All inter-individual variabilities (IIVs) were
##' set at 30%, residual error at 20% and overlapping PK parameters
##' were the same for all models. A similar set of models was
##' previously used to compare NONMEM and Monolix4. Estimates of
##' population parameters, standard errors for fixed-effect
##' parameters, and run times were compared both for closed-form
##' solutions and using ODEs. Additionally, a sparse data estimation
##' situation was investigated where 500 datasets of 600 subjects each
##' (150 per dose) were generated consisting of 4 random time point
##' samples in 24 hours per subject, using a first-order absorption,
##' 1-compartment disposition, linear elimination model.
##'
##' @format A data frame with 7,920 rows and 15 columns
##'
##' \describe{
##'   \item{ID}{Siumlated Subject ID}
##'   \item{TIME}{Simulated Time}
##'   \item{DV}{Simulated Dependent Variable}
##'   \item{LNDV}{Simulated log(Dependent Variable)}
##'   \item{MDV}{Missing DV data item}
##'   \item{AMT}{Dosing AMT}
##'   \item{EVID}{NONMEM Event ID}
##'   \item{DOSE}{Dose}
##'   \item{V}{Individual Central Compartment Volume}
##'   \item{VM}{Individual Vmax}
##'   \item{KM}{Individual Km}
##'   \item{Q}{Individual Q}
##'   \item{V2}{Individual Peripheral Compartment Volume}
##'   \item{SD}{Single Dose Flag}
##'   \item{CMT}{Compartment Indicator}
##' }
##'
##' @source Schoemaker R, Xiong Y, Wilkins J, Laveille C, Wang W.
##'     nlmixr: an open-source package for pharmacometric modelling in
##'     R. ACOP 2016
"Bolus_2CPTMM"
