source("tests/testthat/models/helper-prep_fit.R")
context("NLME61: two-compartment oral, multiple-dose")
runno <- "N061"

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

dat <- datr[datr$SD == 0,]

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
    lik=c(-27010.13, 54042.26, 54113.22),
    param=c(1.3404, 4.1502, 1.6643, 3.8834, -0.18143),
    stdev_param=c(0.31712, 0.25673, 0.42758, 0.23294, 0.32122),
    sigma=0.19942
  )
