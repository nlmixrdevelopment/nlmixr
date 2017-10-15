library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME46: two-compartment infusion, single-dose")
    test_that("Closed-form", {

        datr <-
            read.csv("Infusion_2CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

        datr <- datr[datr$EVID != 2,]

        datIV <- datr[datr$AMT > 0,]
        datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
        datIV$AMT  <- -1*datIV$AMT

        datr <- rbind(datr, datIV)
        datr <- datr[order(datr$ID, datr$TIME),]

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

        runno <- "N046"

        dat <- datr[datr$SD == 1,]

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs6,
                ncmt = 2,
                verbose = TRUE,
                oral = FALSE,
                infusion = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(
                    pnlsTol = .1,
                    msVerbose = TRUE,
                    maxIter = 200
                )
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -11943)
        expect_equal(signif(AIC(fit), 6), 23904.1)
        expect_equal(signif(BIC(fit), 6), 23955.7)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.22)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 1.37)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), 3.89)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.299)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.303)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.235)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.278)

        expect_equal(signif(fit$sigma, 3), 0.199)
    })
    test_that("ODE", {

        datr <-
            read.csv("Infusion_2CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

        datr <- datr[datr$EVID != 2,]

        datIV <- datr[datr$AMT > 0,]
        datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
        datIV$AMT  <- -1*datIV$AMT

        datr <- rbind(datr, datIV)
        datr <- datr[order(datr$ID, datr$TIME),]

        ode2 <- "
    d/dt(centr)  = K21*periph-K12*centr-K10*centr;
    d/dt(periph) =-K21*periph+K12*centr;
    "

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


        runno <- "N046"

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
                control = nlmeControl(pnlsTol = .1, msVerbose = TRUE)
            )

        z <- VarCorr(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik),6), -11943)
        expect_equal(signif(AIC(fitODE), 6), 23904.1)
        expect_equal(signif(BIC(fitODE), 6), 23955.7)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.22)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), 1.37)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]),3), 3.89)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.299)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.303)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.235)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.278)

        expect_equal(signif(fitODE$sigma, 3), 0.199)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
