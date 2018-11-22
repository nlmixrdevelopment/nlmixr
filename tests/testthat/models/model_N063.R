source("tests/testthat/models/helper-prep_fit.R")
context("NLME63: two-compartment oral, multiple-dose")
runno <- "N063"

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

dat <- datr

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
    lik=c(-38499.4, 77020.79, 77096.12),
    param=c(1.3494, 4.2014, 1.4879, 3.8904, -0.10277),
    stdev_param=c(0.31561, 0.26139, 0.38163, 0.27402, 0.32006),
    sigma=0.19978
  )
