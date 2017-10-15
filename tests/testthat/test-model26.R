library(testthat)
library(nlmixr)

rxPermissive({
    context("NLME26: one-compartment oral, multiple-dose")
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

        specs4i <-
            list(
                fixed = lCL + lV + lKA ~ 1,
                random = pdDiag(value = diag(c(3, 3, 3)), form = lCL + lV + lKA ~ 1),
                start = c(lCL = 1.6, lV = 4.5, lKA = 0.2)
            )

        runno <- "N026"

        dat <- datr

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs4i,
                ncmt = 1,
                verbose = TRUE,
                oral = TRUE,
                weight = varPower(fixed = c(1))
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -37484)
        expect_equal(signif(AIC(fit), 6), 74982)
        expect_equal(signif(BIC(fit), 6), 75030)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.2)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), -0.00626)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.259)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.277)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.325)

        expect_equal(signif(fit$sigma, 3), 0.199)

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

        specs4i <-
            list(
                fixed = lCL + lV + lKA ~ 1,
                random = pdDiag(lCL + lV + lKA ~ 1),
                start = c(lCL = 1.6, lV = 4.5, lKA = 0.2)
            )

        runno <- "N026"

        dat <- datr

        fitODE <-
            nlme_ode(
                dat,
                model = ode1KA,
                par_model = specs4i,
                par_trans = mypar4,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .2, msVerbose = TRUE)
            )

        z <- VarCorr(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik),6), -37490.9)
        expect_equal(signif(AIC(fitODE), 6), 74995.8)
        expect_equal(signif(BIC(fitODE), 6), 75043.7)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.2)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), -0.0061)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.260)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.277)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.321)

        expect_equal(signif(fitODE$sigma, 3), 0.198)
    })
}, on.validate="NLMIXR_VALIDATION")
