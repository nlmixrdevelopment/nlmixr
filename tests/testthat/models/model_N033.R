source("tests/testthat/models/helper-prep_fit.R")
context("NLME33: two-compartment bolus, multiple-dose")
runno <- "N033"

datr <- Bolus_2CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]
specs6 <-
  list(
    fixed = lCL + lV + lCLD + lVT ~ 1,
    random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
    start = c(
      lCL = 1.37,
      lV = 4.19,
      lCLD = 1.37,
      lVT = 3.87
    )
  )

dat <- datr[datr$SD == 0,]

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs6,
    ncmt = 2,
    oral = FALSE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-27501.23, 55020.46, 55078.52),
    param=c(1.3326, 4.1936, 1.29, 3.8091),
    stdev_param=c(0.34247, 0.30254, 0.00067117, 0.27539),
    sigma=0.20531
  )
