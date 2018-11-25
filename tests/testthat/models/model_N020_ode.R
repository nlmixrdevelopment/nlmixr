source("helper-prep_fit.R")
context("NLME20: one-compartment infusion, single-dose, Michaelis-Menten")
runno <- "N020_ode"

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

dat <- datr[datr$SD == 1,]

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
    lik=c(-13025, 26063.99, 26104.11),
    param=c(6.8493, 5.3436, 4.2594),
    stdev_param=c(1.52, 0.0012525, 1.4839),
    sigma=c(0.20363)
  )
