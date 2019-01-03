source("helper-prep_fit.R")
context("NLME30: one-compartment oral, Michaelis-Menten, multiple-dose")
runno <- "N030_ode"

ode1MMKA <- "
    d/dt(abs)    = -KA*abs;
d/dt(centr)  =  KA*abs-(VM*centr/V)/(KM+centr/V);
"

mypar5 <- function(lVM, lKM, lV, lKA)
{
  VM <- exp(lVM)
  KM <- exp(lKM)
  V <- exp(lV)
  KA <- exp(lKA)
  
}
specs5i <-
  list(
    fixed = lVM + lKM + lV + lKA ~ 1,
    random = pdDiag(value = diag(c(2, 2, 2, 2)), form = lVM + lKM + lV + lKA ~
                      1),
    start = c(
      lVM = 7,
      lKM = 6.2,
      lV = 4.5,
      lKA = -0.2
    )
  )

datr <-
  read.csv("../Oral_1CPTMM.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

dat <- datr[datr$SD == 0,]

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1MMKA,
    par_model = specs5i,
    par_trans = mypar5,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-30897.83, 61813.66, 61871.72),
    param=c(7.0674, 6.0633, 4.1802, -0.077726),
    stdev_param=c(0.8028, 3.5625, 1.4762, 2.4018e-06),
    sigma=c(0.20274)
  )
