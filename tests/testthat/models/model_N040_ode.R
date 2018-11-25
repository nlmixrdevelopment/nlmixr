source("helper-prep_fit.R")
context("NLME40: two-compartment bolus Michaelis-Menten, single-dose")
runno <- "N040_ode"

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
      lVM = 7,
      lKM = 6,
      lV = 4,
      lCLD = 1.5,
      lVT = 4
    )
  )

dat <- datr[datr$SD == 1,]

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
    lik=c(-11996.52, 24015.04, 24078.1),
    param=c(6.5832, 4.9069, 4.2465, 1.451, 3.992),
    stdev_param=c(1.7026, 6.9182e-16, 1.4727, 1.5418, 1.4888),
    sigma=c(0.19485)
  )
