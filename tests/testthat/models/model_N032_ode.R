source("tests/testthat/models/helper-prep_fit.R")
context("NLME32: two-compartment bolus, single-dose")
runno <- "N032_ode"

datr <- Bolus_2CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

ode2 <- "
d/dt(centr)  = K21*periph-K12*centr-K10*centr;
d/dt(periph) =-K21*periph+K12*centr;
"

mypar6 <- function(lCL, lV, lCLD, lVT)
{
  CL <- exp(lCL)
  V  <- exp(lV)
  CLD <- exp(lCLD)
  VT <- exp(lVT)
  K10 <- CL / V
  K12 <- CLD / V
  K21 <- CLD / VT
}

specs6 <-
  list(
    fixed = lCL + lV + lCLD + lVT ~ 1,
    random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
    start = c(
      lCL = 1.6,
      lV = 4.5,
      lCLD = 1.5,
      lVT = 3.9
    )
  )
dat <- datr[datr$SD == 1,]

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode2,
    par_model = specs6,
    par_trans = mypar6,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values <-
  list(
    lik=c(-13285.23, 26580.47, 26609.12),
    param=c(1.4221, 4.3410),
    stdev_param=c(0.84730, 0),
    sigma=0.43887
  )
