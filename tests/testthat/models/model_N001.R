source("helper-prep_fit.R")
context("NLME01: one-compartment bolus, single-dose")
runno <- "N001"

datr <- Bolus_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2, ]

specs1 <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.6, lV = 4.5)
  )

dat <- datr[datr$SD == 1, ]

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs1,
    ncmt = 1,
    verbose = verbose_minimization,
    oral = FALSE,
    weight = varPower(fixed = c(1)),
    control = default_control
  )

# Generate this with dput(generate_expected_values(fit[[runno]]))
expected_values[[runno]] <-
  list(
    lik = c(-12119.4, 24248.8, 24277.45),
    param = c(1.3641, 4.2025),
    stdev_param = c(1.3436, 1.5226),
    sigma = 0.19888
  )
