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

runno <- "N001_ode"

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

ode1 <- "
d/dt(centr)  = -(CL/V)*centr;
"

mypar1 <- function(lCL, lV)
{
  CL <- exp(lCL)
  V <- exp(lV)
}

fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1,
    par_model = specs1,
    par_trans = mypar1,
    response = "centr",
    response.scaler = "V",
    verbose = verbose_minimization,
    weight = varPower(fixed = c(1)),
    control = default_control
  )

expected_values <-
  list(
    # Generate with dput(round(c(logLik(fit[[runno]]), AIC(fit[[runno]]), BIC(fit[[runno]])), 2))
    lik=c(-13285.23, 26580.47, 26609.12),
    # Generate with dput(unname(signif(fixef(fit[[runno]]), 5)))
    param=c(1.4221, 4.3410),
    # Generate with dput(unname(signif(VarCorr(fit[[runno]])[1:length(fixef(fit[[runno]])), "StdDev"], 5)))
    stdev_param=c(0.84730, 0),
    # Generate with signif(fit[[runno]]$sigma, 5)
    sigma=0.43887
  )
