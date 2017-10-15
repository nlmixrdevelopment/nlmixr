library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME61: two-compartment oral, multiple-dose")
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

        runno <- "N061"

        dat <- datr[datr$SD == 0,]

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

        expect_equal(signif(as.numeric(fit$logLik),6), -27008.6)
        expect_equal(signif(AIC(fit), 6), 54039.2)
        expect_equal(signif(BIC(fit), 6), 54110.1)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.34)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.16)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 1.64)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), 3.87)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[5]),3), -0.164)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.317)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.258)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.410)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.231)
        expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.321)

        expect_equal(signif(fit$sigma, 3), 0.199)
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

        runno <- "N061"

        dat <- datr[datr$SD == 0,]

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

        expect_equal(signif(as.numeric(fitODE$logLik),6), -27007.7)
        expect_equal(signif(AIC(fitODE), 6), 54037.3)
        expect_equal(signif(BIC(fitODE), 6), 54108.3)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.34)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.17)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), 1.61)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]),3), 3.86)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[5]),3), -0.15)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.317)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.258)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.409)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.229)
        expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.321)

        expect_equal(signif(fitODE$sigma, 3), 0.20)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
