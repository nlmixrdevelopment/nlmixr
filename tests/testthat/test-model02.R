library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME02: one-compartment bolus, multiple-dose")
    test_that("Closed-form", {
        datr <- Bolus_1CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        specs1 <-
            list(
                fixed = lCL + lV ~ 1,
                random = pdDiag(lCL + lV ~ 1),
                start = c(lCL = 1.6, lV = 4.5)
            )

        runno <- "N002"

        dat <- datr[datr$SD == 0,]

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs1,
                ncmt = 1,
                verbose = TRUE,
                oral = FALSE,
                weight = varPower(fixed = c(1))
            )

        z <- summary(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -26811.5)
        expect_equal(signif(AIC(fit), 6), 53633.1)
        expect_equal(signif(BIC(fit), 6), 53665.3)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.36)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.20)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.271)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.311)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.205)
    })

    test_that("ODE", {
        datr <- Bolus_1CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2, ]

        specs1 <-
            list(
                fixed = lCL + lV ~ 1,
                random = pdDiag(lCL + lV ~ 1),
                start = c(lCL = 1.6, lV = 4.5)
            )

        runno <- "N002"

        dat <- datr[datr$SD == 0, ]

        ode1 <- "
    d/dt(centr)  = -(CL/V)*centr;
    "

        mypar1 <- function(lCL, lV)
        {
            CL <- exp(lCL)
            V <- exp(lV)
        }

        fitODE <-
            nlme_ode(
                dat,
                model = ode1,
                par_model = specs1,
                par_trans = mypar1,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .01, msVerbose = TRUE)
            )

        z <- summary(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik), 6),-26811.5)
        expect_equal(signif(AIC(fitODE), 6), 53633)
        expect_equal(signif(BIC(fitODE), 6), 53665.2)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]), 3), 1.36)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]), 3), 4.20)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.271)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.311)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.205)


    })
}, on.validate="NLMIXR_VALIDATION_FULL", silent=TRUE)
