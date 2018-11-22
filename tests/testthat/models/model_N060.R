source("tests/testthat/models/helper-prep_fit.R")
context("NLME60: two-compartment oral, single-dose")
runno <- "N060"

datr <- read.csv("Oral_2CPT.csv",
                 header = TRUE,
                 stringsAsFactors = F)
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]
specs8 <-
  list(
    fixed = lCL + lV + lCLD + lVT + lKA ~ 1,
    random = pdDiag(lCL + lV + lCLD + lVT + lKA ~ 1),
    start = c(
      lCL = 1.6,
      lV = 4.5,
      lCLD = 1.5,
      lVT = 3.9,
      lKA = 0.1
    )
  )

dat <- datr[datr$SD == 1,]

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs8,
    ncmt = 2,
    oral = TRUE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values <-
  list(
    lik=c(-11862.13, 23746.25, 23809.31),
    param=c(1.3744, 4.1784, 1.6059, 3.8443, -0.13271),
    stdev_param=c(0.30773, 0.27772, 0.51161, 0.00087904, 0.33321),
    sigma=0.20344
  )
