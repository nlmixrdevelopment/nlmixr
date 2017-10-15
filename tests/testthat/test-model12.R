library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME12: one-compartment infusion, single-dose")
    test_that("Closed-form", {
        datr <-
            read.csv("Infusion_1CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)

        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

        datr <- datr[datr$EVID != 2, ]

        datIV <- datr[datr$AMT > 0, ]
        datIV$TIME <- datIV$TIME + (datIV$AMT / datIV$RATE)
        datIV$AMT  <- -1 * datIV$AMT

        datr <- rbind(datr, datIV)
        datr <- datr[order(datr$ID, datr$TIME), ]

        specs1 <-
            list(
                fixed = lCL + lV ~ 1,
                random = pdDiag(lCL + lV ~ 1),
                start = c(lCL = 1.5, lV = 4)
            )


        runno <- "N012"

        dat <- datr[datr$SD == 1, ]

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs1,
                ncmt = 1,
                verbose = TRUE,
                oral = FALSE,
                infusion = TRUE,
                weight = varPower(fixed = c(1))
            )

        z <- summary(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -11867.4)
        expect_equal(signif(AIC(fit), 6), 23744.7)
        expect_equal(signif(BIC(fit), 6), 23773.4)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.27)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.276)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.307)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.201)
    })
    test_that("ODE", {

        datr <-
            read.csv("Infusion_1CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)

        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

        datr <- datr[datr$EVID != 2, ]

        datIV <- datr[datr$AMT > 0, ]
        datIV$TIME <- datIV$TIME + (datIV$AMT / datIV$RATE)
        datIV$AMT  <- -1 * datIV$AMT

        datr <- rbind(datr, datIV)
        datr <- datr[order(datr$ID, datr$TIME), ]

        ode1 <- "
    d/dt(centr)  = -(CL/V)*centr;
    "

        mypar1 <- function(lCL, lV)
        {
            CL <- exp(lCL)
            V <- exp(lV)
        }

        specs1m <-
            list(
                fixed = lCL + lV ~ 1,
                random = pdDiag(lCL + lV ~ 1),
                start = c(lCL = 1.3, lV = 4)
            )

        runno <- "N012"

        dat <- datr[datr$SD == 1,]

        fitODE <-
            nlme_ode(
                dat,
                model = ode1,
                par_model = specs1m,
                par_trans = mypar1,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .1, msVerbose = TRUE)
            )

        z <- summary(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik),6), -11867.3)
        expect_equal(signif(AIC(fitODE), 6), 23744.6)
        expect_equal(signif(BIC(fitODE), 6), 23773.3)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.27)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.276)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.307)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.201)

    })
}, on.validate="NLMIXR_VALIDATION_FULL", silent=TRUE)
