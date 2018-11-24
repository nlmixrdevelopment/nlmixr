source("tests/testthat/models/helper-prep_fit.R")
context("NLME02: one-compartment bolus, multiple-dose")
runno <- "N002"

datr <- Bolus_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

specs1 <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.6, lV = 4.5)
  )

dat <- datr[datr$SD == 0,]

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs1,
    ncmt = 1,
    oral = FALSE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-26811.54, 53633.07, 53665.32),
    param=c(1.3593, 4.1966),
    stdev_param=c(1.3234, 1.5182),
    sigma=c(0.20503)
  )
