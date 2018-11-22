source("tests/testthat/models/helper-prep_fit.R")
context("NLME69: two-compartment oral Michaelis-Menten, multiple-dose")
runno <- "N069_ode"

datr <-
  read.csv("Oral_2CPTMM.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

ode2MMKA <- "
d/dt(abs)    =-KA*abs;
d/dt(centr)  = KA*abs+K21*periph-K12*centr-(VM*centr/V)/(KM+centr/V);
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
      lVM = 7,
      lKM = 6,
      lV = 4,
      lCLD = 1.5,
      lVT = 4,
      lKA = 0.1
    )
  )

dat <- datr[datr$SD == 0,]

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
expected_values <-
  list(
    lik=c(-13285.23, 26580.47, 26609.12),
    param=c(1.4221, 4.3410),
    stdev_param=c(0.84730, 0),
    sigma=0.43887
  )
