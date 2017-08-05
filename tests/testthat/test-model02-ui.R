library(testthat)
library(nlmixr)

context("UI-NLME02: one-compartment bolus, multiple-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {

    test_that("Closed-form", {

        datr <- Bolus_1CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        uif <- function(){
            ini({
                lCl <- 1.6
                lV <- 4.5
                prop.err <- 0.1
                eta.V ~ 0.1
                eta.Cl ~ 0.1
            })
            model({
                v <- exp(lV + eta.V)
                cl <- exp(lCl + eta.Cl)
                linCmt() ~ prop(prop.err)
            })
        }

        runno <- "N002"

        dat <- datr[datr$SD == 0,]

        fit <- uif %>% nlme(dat)

        z <- summary(as.nlme(fit))

        expect_equal(signif(as.numeric(as.nlme(fit)$logLik),6), -26811.5)
        expect_equal(signif(AIC(as.nlme(fit)), 6), 53633.1)
        expect_equal(signif(BIC(as.nlme(fit)), 6), 53665.3)

        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[1]),3), 1.36)
        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[2]),3), 4.20)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.271)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.311)
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

        uif <- function(){
            ini({
                lCl <- 1.6
                lV <- 4.5
                prop.err <- 0.1
                eta.V ~ 0.1
                eta.cl ~ 0.1
            })
            model({
                v <- exp(lV+ eta.V)
                cl <- exp(lCl + eta.cl)
                d / dt(centr) = -(cl / v) * centr;
                cp = centr / v;
                cp ~ prop(prop.err)
            })
        }

        fitODE <- uif %>% nlme(dat, control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))

        z <- summary(as.nlme(fitODE))

        expect_equal(signif(as.numeric(as.nlme(fitODE)$logLik), 6),-26811.5)
        expect_equal(signif(AIC(as.nlme(fitODE)), 6), 53633)
        expect_equal(signif(BIC(as.nlme(fitODE)), 6), 53665.2)

        expect_equal(signif(as.numeric(as.nlme(fitODE)$coefficients$fixed[1]), 3), 1.36)
        expect_equal(signif(as.numeric(as.nlme(fitODE)$coefficients$fixed[2]), 3), 4.20)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.271)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.311)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.205)

    })
}
