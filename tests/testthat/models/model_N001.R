# This header variable assignment simplifies testing and comparison across
# platforms
if (!exists("fit")) {
  fit <- list()
  z <- list()
}
if (!exists("verbose_minimization")) verbose_minimization <- FALSE
default_control <-
  nlmeControl(
    returnObject=TRUE,
    msMaxIter=1L,
    maxIter=1L,
    pnlsMaxIter=1L
  )

runno <- "N001"

datr <- Bolus_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2, ]

specs1 <-
  list(
    fixed = lCL + lV ~ 1,
    random = pdDiag(lCL + lV ~ 1),
    start = c(lCL = 1.6, lV = 4.5)
  )

dat <- datr[datr$SD == 1, ]

fit[[runno]] <-
  nlme_lin_cmpt(
    dat,
    par_model = specs1,
    ncmt = 1,
    verbose = verbose_minimization,
    oral = FALSE,
    weight = varPower(fixed = c(1)),
    control = default_control
  )

expected_values <-
  list(
    lik=c(-13285.15, 26580.30, 26608.95),
    param=c(1.4221, 4.3410),
    stdev_param=c(0.84660, 0),
    sigma=0.43887
  )
