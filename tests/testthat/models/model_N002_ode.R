source("tests/testthat/models/helper-prep_fit.R")
context("NLME02: one-compartment bolus, multiple-dose")
runno <- "N002_ode"

datr <- Bolus_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2, ]

specs1 <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.6, lV = 4.5)
  )

dat <- datr[datr$SD == 0, ]

ode1 <- "
    d/dt(centr)  = -(CL/V)*centr;
    "

mypar1 <- function(lCL, lV)
{
  CL <- exp(lCL)
  V <- exp(lV)
}

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1,
    par_model = specs1,
    par_trans = mypar1,
    response = "centr",
    response.scaler = "V",
    verbose = verbose_minimization,
    weight = varPower(fixed = c(1)),
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values <-
  list(
    lik=c(-26811.52, 53633.04, 53665.29),
    param=c(1.3594, 4.1967),
    stdev_param=c(1.3233, 1.5183),
    sigma=c(0.20503)
  )
