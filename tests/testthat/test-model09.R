library(testthat)
library(nlmixr)

rxPermissive({
    context("NLME09: one-compartment bolus, Michaelis-Menten, single-dose")
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

        datr <- Bolus_1CPTMM
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)

        datr <- datr[datr$EVID != 2,]

        runno <- "N009"

        dat <- datr[datr$SD == 1,]

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

        expect_equal(signif(as.numeric(fit$logLik), 6), -12566.5)
        expect_equal(signif(AIC(fit), 6), 25147)
        expect_equal(signif(BIC(fit), 6), 25187.1)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]), 3), 6.80)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]), 3), 5.33)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]), 3), 4.18)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.301)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.259)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.297)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[4], 3)), 0.200)

})
}, on.validate="NLMIXR_VALIDATION_FULL", silent=TRUE)
