source("tests/testthat/models/helper-prep_fit.R")
context("NLME26: one-compartment oral, multiple-dose")
runno <- "N026"

datr <- Oral_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

specs4i <-
    list(
      fixed = lCL + lV + lKA ~ 1,
      random = pdDiag(value = diag(c(3, 3, 3)), form = lCL + lV + lKA ~ 1),
      start = c(lCL = 1.6, lV = 4.5, lKA = 0.2)
     )

dat <- datr
        
fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs4i,
    ncmt = 1,
    oral = TRUE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-37484.02, 74982.05, 75029.97),
    param=c(1.3882, 4.2007, -0.0062572),
    stdev_param=c(0.25942, 0.27712, 0.32461),
    sigma=0.19852
  )
