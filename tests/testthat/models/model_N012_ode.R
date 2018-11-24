source("helper-prep_fit.R")
context("NLME12: one-compartment infusion, single-dose")
runno <- "N012.ode"

datr <-
  read.csv("Infusion_1CPT.csv",
           header = TRUE,
           stringsAsFactors = F)

datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

datr <- datr[datr$EVID != 2, ]

datIV <- datr[datr$AMT > 0, ]
datIV$TIME <- datIV$TIME + (datIV$AMT / datIV$RATE)
datIV$AMT  <- -1 * datIV$AMT

datr <- rbind(datr, datIV)
datr <- datr[order(datr$ID, datr$TIME), ]

ode1 <- "
d/dt(centr)  = -(CL/V)*centr;
"

mypar1 <- function(lCL, lV)
{
  CL <- exp(lCL)
  V <- exp(lV)
}

specs1m <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.3, lV = 4)
  )

dat <- datr[datr$SD == 1,]

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
    lik=c(-13285.23, 26580.47, 26609.12),
    param=c(1.4221, 4.3410),
    stdev_param=c(0.84730, 0),
    sigma=0.43887
  )
