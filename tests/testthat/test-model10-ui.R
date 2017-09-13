library(testthat)
library(nlmixr)

context("UI-NLME10: one-compartment bolus, Michaelis-Menten, multiple-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {

    test_that("ODE", {

        uif <- function(){
            ini({
                lV <- 4
                lVm <- 7
                lKm <- 6
                prop.err <- 3
                eta.Vm ~ 0.2
                eta.Km ~ 0.2
                eta.V ~ 0.2
            })
            model({
                VM <- exp(lVm + eta.Vm)
                KM <- exp(lKm + eta.Km)
                V <- exp(lV + eta.V)
                d/dt(centr) <- -(VM*centr/V)/(KM+centr/V);
                cp = centr / V
                cp ~ prop(prop.err)
            })
        }

        datr <- Bolus_1CPTMM
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)

        datr <- datr[datr$EVID != 2,]

        runno <- "N010"

        dat <- datr[datr$SD == 0,]

        fit <- nlmixr(uif, datr, est="nlme", control = nlmeControl(pnlsTol = .1, msVerbose = TRUE))

        z <- summary(as.nlme(fit))

        expect_equal(signif(as.numeric(as.nlme(fit)$logLik), 6), -43484.8)
        expect_equal(signif(AIC(as.nlme(fit)), 6), 86983.5)
        expect_equal(signif(BIC(as.nlme(fit)), 6), 87031.5)

        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[1]), 3), 4.18)
        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[2]), 3), 6.9)
        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[3]), 3), 5.53)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.293)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.295)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.3)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[4], 3)), 0.202)

    })
}


