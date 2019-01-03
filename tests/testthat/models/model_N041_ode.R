source("helper-prep_fit.R")
context("NLME41: two-compartment bolus Michaelis-Menten, multiple-dose")
runno <- "N041_ode"

datr <- Bolus_2CPTMM
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

ode2MM <- "
d/dt(centr)  = K21*periph-K12*centr-(VM*centr/V)/(KM+centr/V);
d/dt(periph) =-K21*periph+K12*centr;
"

mypar7 <- function(lVM, lKM, lV, lCLD, lVT)
{
  VM <- exp(lVM)
  KM <- exp(lKM)
  V <- exp(lV)
  CLD  <- exp(lCLD)
  VT <- exp(lVT)
  K12 <- CLD / V
  K21 <- CLD / VT
}
specs7 <-
  list(
    fixed = lVM + lKM + lV + lCLD + lVT ~ 1,
    random = pdDiag(lVM + lKM + lV + lCLD + lVT ~ 1),
    start = c(
      lVM = 7.6,
      lKM = 6.83,
      lV = 4.29,
      lCLD = 1.26,
      lVT = 3.58
    )
  )

dat <- datr[datr$SD == 0,]

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode2MM,
    par_model = specs7,
    par_trans = mypar7,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-29949.65, 59921.3, 59992.26),
    param=c(7.6111, 6.9606, 4.252, 1.4127, 3.76),
    stdev_param=c(0.0032084, 3.6923, 1.378, 2.6674e-05, 1.4656),
    sigma=c(0.21318)
  )
