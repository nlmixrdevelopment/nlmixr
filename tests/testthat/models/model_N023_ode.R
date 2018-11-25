source("helper-prep_fit.R")
context("NLME23: one-compartment oral, single-dose")
runno <- "N023_ode"

datr <- Oral_1CPT
datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
datr <- datr[datr$EVID != 2,]

ode1KA <- "
d/dt(abs)    = -KA*abs;
d/dt(centr)  =  KA*abs-(CL/V)*centr;
"

        mypar4 <- function(lCL, lV, lKA)
        {
            CL <- exp(lCL)
            V <- exp(lV)
            KA <- exp(lKA)
        }

        specs4 <-
            list(
                fixed = lCL + lV + lKA ~ 1,
                random = pdDiag(lCL + lV + lKA ~ 1),
                start = c(lCL = 1, lV = 4, lKA = 0)
            )

        dat <- datr[datr$SD == 1,]
        
fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1KA,
    par_model = specs4,
    par_trans = mypar4,
    response = "centr",
    response.scaler = "V",
    weight = varPower(fixed = c(1)),
    verbose = verbose_minimization,
    control = default_control
  )

# Generate this with generate_expected_values(fit[[runno]])
expected_values[[runno]] <-
  list(
    lik=c(-11840.4, 23694.8, 23734.9),
    param=c(1.3863, 4.1879, -0.028869),
    stdev_param=c(1.2515, 1.3394, 1.5383),
    sigma=c(0.19919)
  )
