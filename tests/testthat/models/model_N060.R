source("helper-prep_fit.R")
context("NLME60: two-compartment oral, single-dose")
runno <- "N060"

datr <- read.csv("../Oral_2CPT.csv",
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
expected_values[[runno]] <-
  list(
    lik=c(-11842.77, 23707.54, 23770.59),
    param=c(1.3771, 4.1467, 1.7061, 3.886, -0.17619),
    stdev_param=c(1.5404, 1.3587, 2.1849, 1.0597, 1.6278),
    sigma=c(0.19981)
  )
