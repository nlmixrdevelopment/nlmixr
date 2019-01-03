source("helper-prep_fit.R")
context("NLME01: one-compartment bolus, single-dose")
runno <- "N001_ode"

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
expected_values[[runno]] <-
  list(
    lik = c(-12119.41, 24248.81, 24277.46),
    param = c(1.3641, 4.2025),
    stdev_param = c(1.3436, 1.5226),
    sigma = 0.19888
  )
