source("helper-prep_fit.R")
context("NLME22: one-compartment infusion, multiple-dose, Michaelis-Menten")
runno <- "N022_ode"

datr <-
  read.csv("../Infusion_1CPTMM.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

datr <- subset(datr, EVID != 2)
datIV <- subset(datr, AMT>0)
datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
datIV$AMT  <- -1*datIV$AMT
#datIV <- datr[AMT > 0][, TIME := TIME + AMT / RATE][, AMT := -1 * AMT]
datr <- rbind(datr, datIV)
datr <- datr[order(datr$ID, datr$TIME),]

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

dat <- datr

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
    lik=c(-45848.64, 91711.28, 91759.21),
    param=c(6.9242, 5.5095, 4.249),
    stdev_param=c(1.2636, 1.3739, 1.4725),
    sigma=c(0.20421)
  )
