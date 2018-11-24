source("helper-prep_fit.R")
context("NLME10: one-compartment bolus, Michaelis-Menten, multiple-dose")
runno <- "N010.ode"

ode1MM <- "
d/dt(centr)  = -(VM*centr/V)/(KM+centr/V);
"

mypar3 <- function(lVM, lKM, lV)
{
  VM <- exp(lVM)
  KM <- exp(lKM)
  V <- exp(lV)
}
specs3 <-
  list(
    fixed = lVM + lKM + lV ~ 1,
    random = pdDiag(lVM + lKM + lV ~ 1),
    start = c(lVM = 7, lKM = 6, lV = 4)
  )

datr <-Bolus_1CPTMM
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)

datr <- datr[datr$EVID != 2,]

dat <- datr[datr$SD == 0,]

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1MM,
    par_model = specs3,
    par_trans = mypar3,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-31282.25, 62578.5, 62623.65),
    param=c(7.003, 5.7414, 4.1216),
    stdev_param=c(1.2782, 1.567, 1.4919),
    sigma=c(0.20506)
  )
