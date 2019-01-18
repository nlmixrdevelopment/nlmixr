source("helper-prep_fit.R")
context("NLME68: two-compartment oral Michaelis-Menten, single-dose")
runno <- "N068_ode"

datr <-
  read.csv("../Oral_2CPTMM.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

ode2MMKA <- "
d/dt(abs)    =-KA*abs;
d/dt(centr)  = KA*abs+K21*periph-K12*centr-exp(log(VM)+log(centr)-log(V)-log(KM+centr/V));
d/dt(periph) =-K21*periph+K12*centr;
"

mypar9 <- function(lVM, lKM, lV, lCLD, lVT, lKA)
{
  VM <- exp(lVM)
  KM <- exp(lKM)
  V <- exp(lV)
  CLD  <- exp(lCLD)
  VT <- exp(lVT)
  KA <- exp(lKA)
  K12 <- CLD / V
  K21 <- CLD / VT
}
specs9 <-
  list(
    fixed = lVM + lKM + lV + lCLD + lVT + lKA ~ 1,
    random = pdDiag(lVM + lKM + lV + lCLD + lVT + lKA ~ 1),
    start = c(
      lVM = 7.3,
      lKM = 6.1,
      lV = 4.3,
      lCLD = 1.4,
      lVT = 3.75,
      lKA = -0.01
    )
  )

dat <- datr[datr$SD == 1,]

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode2MMKA,
    par_model = specs9,
    par_trans = mypar9,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-11772.86, 23571.72, 23646.24),
    param=c(7.6964, 6.5899, 4.2466, 1.3497, 3.6854, -0.023332),
    stdev_param=c(1.7751e-06, 2.0083, 1.6421, 0.0011, 1.2913, 1.4265),
    sigma=c(0.20446)
  )
