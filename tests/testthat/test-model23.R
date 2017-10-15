library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME23: one-compartment oral, single-dose")
    test_that("Closed-form", {

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

        runno <- "N023"

        dat <- datr[datr$SD == 1,]

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs4,
                ncmt = 1,
                verbose = TRUE,
                oral = TRUE,
                weight = varPower(fixed = c(1))
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -11768.9)
        expect_equal(signif(AIC(fit), 6), 23551.9)
        expect_equal(signif(BIC(fit), 6), 23592)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.19)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), -0.0308)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.252)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.27)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.312)

        expect_equal(signif(fit$sigma, 3), 0.200)
    })
    test_that("ODE", {

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

        runno <- "N023"

        dat <- datr[datr$SD == 1,]

        fitODE <-
            nlme_ode(
                dat,
                model = ode1KA,
                par_model = specs4,
                par_trans = mypar4,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(
                    pnlsTol = .3,
                    tolerance = 1e-3,
                    msVerbose = TRUE
                )
            )

        z <- VarCorr(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik),6), -11775.3)
        expect_equal(signif(AIC(fitODE), 6), 23564.6)
        expect_equal(signif(BIC(fitODE), 6), 23604.7)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.19)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), -0.0302)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.252)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.27)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.309)

        expect_equal(signif(fitODE$sigma, 3), 0.200)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
