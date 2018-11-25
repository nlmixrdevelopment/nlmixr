source("helper-prep_fit.R")
context("NLME61: two-compartment oral, multiple-dose")
runno <- "N061_ode"

datr <- read.csv("../Oral_2CPT.csv",
                 header = TRUE,
                 stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]


ode2KA <- "
d/dt(abs)    = -KA*abs;
d/dt(centr)  =  KA*abs+K21*periph-K12*centr-K10*centr;
d/dt(periph) =        -K21*periph+K12*centr;
"

mypar8 <- function(lCL, lV, lCLD, lVT, lKA)
{
  CL <- exp(lCL)
  V  <- exp(lV)
  CLD <- exp(lCLD)
  VT <- exp(lVT)
  KA <- exp(lKA)
  K10 <- CL / V
  K12 <- CLD / V
  K21 <- CLD / VT
}

specs8 <-
  list(
    fixed = lCL + lV + lCLD + lVT + lKA ~ 1,
    random = pdDiag(lCL + lV + lCLD + lVT + lKA ~ 1),
    start = c(
      lCL = 1.6,
      lV = 4.5,
      lCLD = 1.5,
      lVT = 3.9,
      lKA = 0.1
    )
  )

dat <- datr[datr$SD == 0,]

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode2KA,
    par_model = specs8,
    par_trans = mypar8,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-13285.23, 26580.47, 26609.12),
    param=c(1.4221, 4.3410),
    stdev_param=c(0.84730, 0),
    sigma=0.43887
  )
