source("tests/testthat/models/helper-prep_fit.R")
context("NLME12: one-compartment infusion, single-dose")
runno <- "N012"

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

specs1 <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.5, lV = 4)
  )

dat <- datr[datr$SD == 1, ]

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs1,
    ncmt = 1,
    oral = FALSE,
    infusion = TRUE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-11867.37, 23744.74, 23773.39),
    param=c(1.3898, 4.2657),
    stdev_param=c(0.27648, 0.30737),
    sigma=0.20145
  )
