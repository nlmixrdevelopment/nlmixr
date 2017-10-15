library(testthat)
library(nlmixr)

rxPermissive({
    context("NLME10: one-compartment bolus, Michaelis-Menten, multiple-dose")
    test_that("ODE", {
        ode1MM <- "
d/dt(centr)  = -(VM*centr/V)/(KM+centr/V);
  "

        mypar3 <- function(lVM, lKM, lV)
        {
            VM <- exp(lVM)
            KM <- exp(lKM)
            V <- exp(lV)
        }
        specs3 <-
            list(
                fixed = lVM + lKM + lV ~ 1,
                random = pdDiag(lVM + lKM + lV ~ 1),
                start = c(lVM = 7, lKM = 6, lV = 4)
            )

        datr <-Bolus_1CPTMM
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)

        datr <- datr[datr$EVID != 2,]

        runno <- "N010"

        dat <- datr[datr$SD == 0,]

        fit <-
            nlme_ode(
                dat,
                model = ode1MM,
                par_model = specs3,
                par_trans = mypar3,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .01, msVerbose = TRUE)
            )

        z <- summary(fit)

        expect_equal(signif(as.numeric(fit$logLik), 6), -31282.7)
        expect_equal(signif(AIC(fit), 6), 62579.4)
        expect_equal(signif(BIC(fit), 6), 62624.5)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]), 3), 7.00)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]), 3), 5.74)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]), 3), 4.12)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.261)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.305)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.306)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[4], 3)), 0.205)
})
}, on.validate="NLMIXR_VALIDATION_FULL", silent=TRUE)
