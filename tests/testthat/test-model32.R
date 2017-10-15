library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME32: two-compartment bolus, single-dose")
    test_that("Closed-form", {

        datr <- Bolus_2CPT

        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]


        ode2 <- "
    d/dt(centr)  = K21*periph-K12*centr-K10*centr;
    d/dt(periph) =-K21*periph+K12*centr;
    "

        mypar6 <- function(lCL, lV, lCLD, lVT)
        {
            CL <- exp(lCL)
            V  <- exp(lV)
            CLD <- exp(lCLD)
            VT <- exp(lVT)
            K10 <- CL / V
            K12 <- CLD / V
            K21 <- CLD / VT
        }

        specs6 <-
            list(
                fixed = lCL + lV + lCLD + lVT ~ 1,
                random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
                start = c(
                    lCL = 1.37,
                    lV = 4.19,
                    lCLD = 1.37,
                    lVT = 3.87
                )
            )

        runno <- "N032"

        dat <- datr[datr$SD == 1,]

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs6,
                ncmt = 2,
                verbose = TRUE,
                oral = FALSE,
                weight = varPower(fixed = c(1)),
                control=nlmeControl(pnlsTol=0.01)
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -12177.9)
        expect_equal(signif(AIC(fit), 6), 24373.8)
        expect_equal(signif(BIC(fit), 6), 24425.3)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.37)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.18)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 1.43)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), 3.88)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.316)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.314)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.299)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.252)

        expect_equal(signif(fit$sigma, 3), 0.198)
    })
    test_that("ODE", {

        datr <- Bolus_2CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        ode2 <- "
    d/dt(centr)  = K21*periph-K12*centr-K10*centr;
    d/dt(periph) =-K21*periph+K12*centr;
    "

        mypar6 <- function(lCL, lV, lCLD, lVT)
        {
            CL <- exp(lCL)
            V  <- exp(lV)
            CLD <- exp(lCLD)
            VT <- exp(lVT)
            K10 <- CL / V
            K12 <- CLD / V
            K21 <- CLD / VT
        }

        specs6 <-
            list(
                fixed = lCL + lV + lCLD + lVT ~ 1,
                random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
                start = c(
                    lCL = 1.6,
                    lV = 4.5,
                    lCLD = 1.5,
                    lVT = 3.9
                )
            )

        runno <- "N032"

        dat <- datr[datr$SD == 1,]

        fitODE <-
            nlme_ode(
                dat,
                model = ode2,
                par_model = specs6,
                par_trans = mypar6,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .3, msVerbose = TRUE)
            )

        z <- VarCorr(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik),6), -12177.3)
        expect_equal(signif(AIC(fitODE), 6), 24372.6)
        expect_equal(signif(BIC(fitODE), 6), 24424.2)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.37)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.2)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), 1.33)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]),3), 3.83)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.317)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.314)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.284)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.244)

        expect_equal(signif(fitODE$sigma, 3), 0.198)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
