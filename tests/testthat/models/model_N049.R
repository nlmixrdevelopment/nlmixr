source("helper-prep_fit.R")
context("NLME49: two-compartment infusion, multiple-dose")
runno <- "N049"

datr <-
  read.csv("../Infusion_2CPT.csv",
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
      lCL = 1.36,
      lV = 4.2,
      lCLD = 1.47,
      lVT = 3.9
    )
  )
    
    dat <- datr
    
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
    lik=c(-38881.34, 77780.69, 77842.32),
    param=c(1.3627, 4.2142, 1.3564, 3.9297),
    stdev_param=c(0.2983, 0.3034, 0.27183, 0.31848),
    sigma=0.20382
  )
