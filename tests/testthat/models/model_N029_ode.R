source("tests/testthat/models/helper-prep_fit.R")
context("NLME29: one-compartment oral, Michaelis-Menten, multiple-dose")
runno <- "N029_ode"

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
specs5 <-
  list(
    fixed = lVM + lKM + lV + lKA ~ 1,
    random = pdDiag(lVM + lKM + lV + lKA ~ 1),
    start = c(
      lVM = 7,
      lKM = 6,
      lV = 4,
      lKA = 0
    )
  )
datr <-
  read.csv("Oral_1CPTMM.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

dat <- datr[datr$SD == 1,]

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1MMKA,
    par_model = specs5,
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
    lik=c(-13285.23, 26580.47, 26609.12),
    param=c(1.4221, 4.3410),
    stdev_param=c(0.84730, 0),
    sigma=0.43887
  )
