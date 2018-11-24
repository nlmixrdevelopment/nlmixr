source("helper-prep_fit.R")
context("NLME46: two-compartment infusion, single-dose")
runno <- "N046"

datr <-
  read.csv("Infusion_2CPT.csv",
           header = TRUE,
           stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

datr <- datr[datr$EVID != 2,]

datIV <- datr[datr$AMT > 0,]
datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
datIV$AMT  <- -1*datIV$AMT

datr <- rbind(datr, datIV)
datr <- datr[order(datr$ID, datr$TIME),]

specs6 <-
  list(
    fixed = lCL + lV + lCLD + lVT ~ 1,
    random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
    start = c(
      lCL = 1.6,
      lV = 4.5,
      lCLD = 1.5,
      lVT = 3.9
    )
  )

dat <- datr[datr$SD == 1,]

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs6,
    ncmt = 2,
    oral = FALSE,
    infusion = TRUE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-11943.84, 23905.67, 23957.26),
    param=c(1.386, 4.2173, 1.3682, 3.8964),
    stdev_param=c(0.29887, 0.30316, 0.24351, 0.27795),
    sigma=0.19856
  )
