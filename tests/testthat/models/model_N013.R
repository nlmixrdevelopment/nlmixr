source("tests/testthat/models/helper-prep_fit.R")
context("NLME13: one-compartment infusion, multiple-dose")
runno <- "N013"

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

dat <- datr[datr$SD == 0, ]

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
expected_values <-
  list(
    lik=c(-26406.19, 52822.38, 52854.63),
    param=c(1.3855, 4.2658),
    stdev_param=c(0.28146, 0.29631),
    sigma=0.20063
  )
