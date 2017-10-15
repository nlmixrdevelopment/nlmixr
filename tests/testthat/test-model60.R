library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME60: two-compartment oral, single-dose")
    test_that("Closed-form", {

        datr <- read.csv("Oral_2CPT.csv",
                         header = TRUE,
                         stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]


        ode2KA <- "
    d/dt(abs)    = -KA*abs;
    d/dt(centr)  =  KA*abs+K21*periph-K12*centr-K10*centr;
    d/dt(periph) =        -K21*periph+K12*centr;
    "

        mypar8 <- function(lCL, lV, lCLD, lVT, lKA)
        {
            CL <- exp(lCL)
            V  <- exp(lV)
            CLD <- exp(lCLD)
            VT <- exp(lVT)
            KA <- exp(lKA)
            K10 <- CL / V
            K12 <- CLD / V
            K21 <- CLD / VT
        }

        specs8 <-
            list(
                fixed = lCL + lV + lCLD + lVT + lKA ~ 1,
                random = pdDiag(lCL + lV + lCLD + lVT + lKA ~ 1),
                start = c(
                    lCL = 1.6,
                    lV = 4.5,
                    lCLD = 1.5,
                    lVT = 3.9,
                    lKA = 0.1
                )
            )

        runno <- "N060"

        dat <- datr[datr$SD == 1,]

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs8,
                ncmt = 2,
                verbose = TRUE,
                oral = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .1, msVerbose = TRUE)
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -11842.3)
        expect_equal(signif(AIC(fit), 6), 23706.5)
        expect_equal(signif(BIC(fit), 6), 23769.6)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.37)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.19)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 1.56)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), 3.85)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[5]),3), -0.123)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.310)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.267)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.446)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.240)
        expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.317)

        expect_equal(signif(fit$sigma, 3), 0.2)
    })

    test_that("ODE", {

        datr <- read.csv("Oral_2CPT.csv",
                         header = TRUE,
                         stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]


        ode2KA <- "
    d/dt(abs)    = -KA*abs;
    d/dt(centr)  =  KA*abs+K21*periph-K12*centr-K10*centr;
    d/dt(periph) =        -K21*periph+K12*centr;
    "

        mypar8 <- function(lCL, lV, lCLD, lVT, lKA)
        {
            CL <- exp(lCL)
            V  <- exp(lV)
            CLD <- exp(lCLD)
            VT <- exp(lVT)
            KA <- exp(lKA)
            K10 <- CL / V
            K12 <- CLD / V
            K21 <- CLD / VT
        }

        specs8 <-
            list(
                fixed = lCL + lV + lCLD + lVT + lKA ~ 1,
                random = pdDiag(lCL + lV + lCLD + lVT + lKA ~ 1),
                start = c(
                    lCL = 1.6,
                    lV = 4.5,
                    lCLD = 1.5,
                    lVT = 3.9,
                    lKA = 0.1
                )
            )

        runno <- "N060"

        dat <- datr[datr$SD == 1,]

        fitODE <-
            nlme_ode(
                dat,
                model = ode2KA,
                par_model = specs8,
                par_trans = mypar8,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .3, msVerbose = TRUE)
            )

        z <- VarCorr(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik),6), -11844.4)
        expect_equal(signif(AIC(fitODE), 6), 23710.7)
        expect_equal(signif(BIC(fitODE), 6), 23773.8)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.38)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.18)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), 1.63)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]),3), 3.84)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[5]),3), -0.139)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.311)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.264)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.436)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.256)
        expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.319)

        expect_equal(signif(fitODE$sigma, 3), 0.199)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
