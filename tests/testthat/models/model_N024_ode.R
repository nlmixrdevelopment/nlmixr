source("helper-prep_fit.R")
context("NLME24: one-compartment oral, multiple-dose")
runno <- "N024_ode"

datr <-Oral_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

ode1KA <- "
d/dt(abs)    = -KA*abs;
d/dt(centr)  =  KA*abs-(CL/V)*centr;
"

mypar4 <- function(lCL, lV, lKA)
{
  CL <- exp(lCL)
  V <- exp(lV)
  KA <- exp(lKA)
}

specs4i <-
  list(
    fixed = lCL + lV + lKA ~ 1,
    random = pdDiag(lCL + lV + lKA ~ 1),
    start = c(lCL = 1, lV = 4, lKA = 0)
  )

dat <- datr[datr$SD == 0,]

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1KA,
    par_model = specs4i,
    par_trans = mypar4,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-26158.59, 52331.17, 52376.32),
    param=c(1.3892, 4.2036, -0.013188),
    stdev_param=c(1.3388, 1.4266, 1.6809),
    sigma=c(0.19738)
  )
