source("helper-prep_fit.R")
context("NLME35: two-compartment oral, multiple-dose")
runno <- "N035"

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

dat <- datr

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
    lik=c(-39312.84, 78643.67, 78705.3),
    param=c(1.3429, 4.1966, 1.3044, 3.896),
    stdev_param=c(0.33882, 0.31267, 0.00056146, 0.30415),
    sigma=0.20288
  )
