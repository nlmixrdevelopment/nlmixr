source("helper-prep_fit.R")
context("NLME26: one-compartment oral, multiple-dose")
runno <- "N026_ode"

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

        specs4i <-
            list(
                fixed = lCL + lV + lKA ~ 1,
                random = pdDiag(lCL + lV + lKA ~ 1),
                start = c(lCL = 1.6, lV = 4.5, lKA = 0.2)
            )

        dat <- datr
        
fit[[runno]] <-
  nlme_ode(
    dat,
    model = ode1KA,
    par_model = specs4i,
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
    lik=c(-37561.87, 75137.75, 75185.67),
    param=c(1.3884, 4.2011, -0.0049212),
    stdev_param=c(1.3111, 1.3988, 1.6119),
    sigma=c(0.19818)
  )
