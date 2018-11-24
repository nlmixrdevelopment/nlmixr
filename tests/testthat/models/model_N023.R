source("helper-prep_fit.R")
context("NLME23: one-compartment oral, single-dose")
runno <- "N023"

datr <- Oral_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

        specs4 <-
            list(
                fixed = lCL + lV + lKA ~ 1,
                random = pdDiag(lCL + lV + lKA ~ 1),
                start = c(lCL = 1, lV = 4, lKA = 0)
            )

        dat <- datr[datr$SD == 1,]
        
fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs4,
    ncmt = 1,
    oral = TRUE,
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-11768.94, 23551.88, 23591.99),
    param=c(1.3857, 4.1871,-0.030783),
    stdev_param=c(0.25193, 0.26979,0.31213),
    sigma=0.2002
  )
