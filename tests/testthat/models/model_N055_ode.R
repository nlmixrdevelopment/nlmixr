source("tests/testthat/models/helper-prep_fit.R")
context("NLME55: two-compartment infusion Michaelis-Menten, multiple-dose")
runno <- "N055_ode"

datr <-
  read.csv("Infusion_2CPTMM.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

datr <- datr[datr$EVID != 2,]
datIV <- datr[datr$AMT > 0,]
datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
datIV$AMT <- -1*datIV$AMT
datr <- rbind(datr, datIV)
datr <- datr[order(datr$ID, datr$TIME),]

dat <- datr[datr$SD == 0,]

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

specs7i <-
  list(
    fixed = lVM + lKM + lV + lCLD + lVT ~ 1,
    random = pdDiag(value = diag(c(3, 3, 3, 3, 3)), form = lVM + lKM + lV +
                      lCLD + lVT ~ 1),
    start = c(
      lVM = 7,
      lKM = 5,
      lV = 4,
      lCLD = 1.2,
      lVT = 4
    )
  )

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode2MM,
    par_model = specs7i,
    par_trans = mypar7,
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
