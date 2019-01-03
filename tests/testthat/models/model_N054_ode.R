source("helper-prep_fit.R")
context("NLME54: two-compartment infusion Michaelis-Menten, single-dose")
runno <- "N054_ode"

datr <-
  read.csv("../Infusion_2CPTMM.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

datr <- datr[datr$EVID != 2,]
datIV <- datr[datr$AMT > 0,]
datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
datIV$AMT <- -1*datIV$AMT
datr <- rbind(datr, datIV)
datr <- datr[order(datr$ID, datr$TIME),]

dat <- datr[datr$SD == 1,]

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
    lik=c(-12629.59, 25281.17, 25344.22),
    param=c(7.9617, 6.941, 4.2537, 1.2796, 3.6305),
    stdev_param=c(0.00013892, 2.501, 1.4523, 0.0011695, 1.4857),
    sigma=c(0.2041)
  )
