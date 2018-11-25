source("helper-prep_fit.R")
context("NLME15: one-compartment infusion, multiple-dose")
runno <- "N015_ode"

datr <-
  read.csv("../Infusion_1CPT.csv",
           header = TRUE,
           stringsAsFactors = F)

datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

datr <- datr[datr$EVID != 2, ]

datIV <- datr[datr$AMT > 0, ]
datIV$TIME <- datIV$TIME + (datIV$AMT / datIV$RATE)
datIV$AMT  <- -1 * datIV$AMT

datr <- rbind(datr, datIV)
datr <- datr[order(datr$ID, datr$TIME), ]

specs1m <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.3, lV = 4)
  )

ode1 <- "
    d/dt(centr)  = -(CL/V)*centr;
    "

mypar1 <- function(lCL, lV)
{
  CL <- exp(lCL)
  V <- exp(lV)
}

dat <- datr

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1,
    par_model = specs1m,
    par_trans = mypar1,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-37861.07, 75732.14, 75766.37),
    param=c(1.387, 4.2689),
    stdev_param=c(1.3973, 1.5217),
    sigma=c(0.20047)
  )
