library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME33: two-compartment bolus, multiple-dose")
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

        runno <- "N033"

        dat <- datr[datr$SD == 0,]

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs6,
                ncmt = 2,
                verbose = TRUE,
                oral = FALSE,
                weight = varPower(fixed = c(1))
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -27501.2)
        expect_equal(signif(AIC(fit), 6), 55020.5)
        expect_equal(signif(BIC(fit), 6), 55078.5)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.33)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.19)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 1.29)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), 3.81)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.342)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.303)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.000671)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.275)

        expect_equal(signif(fit$sigma, 3), 0.205)
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
                    lCL = 1.3,
                    lV = 4.19,
                    lCLD = 1.5,
                    lVT = 3.9
                )
            )

        runno <- "N033"

        dat <- datr[datr$SD == 0,]

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

        expect_equal(signif(as.numeric(fitODE$logLik),6), -27502.7)
        expect_equal(signif(AIC(fitODE), 6), 55023.4)
        expect_equal(signif(BIC(fitODE), 6), 55081.5)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.33)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.19)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), 1.3)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]),3), 3.82)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.343)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.303)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.000625)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.277)

        expect_equal(signif(fitODE$sigma, 3), 0.205)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
