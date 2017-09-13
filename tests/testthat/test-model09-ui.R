library(testthat)
library(nlmixr)

context("UI-NLME09: one-compartment bolus, Michaelis-Menten, single-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {

    test_that("ODE", {

        uif <- function(){
            ini({
                lVM = 7
                lKM = 6
                lV = 4
                et.VM ~ 0.1
                et.V ~ 0.1
                prop.err = 0.1
            })
            model({
                VM <- exp(lVM + et.VM)
                KM <- exp(lKM)
                V <- exp(lV + et.V)
                d/dt(centr)  = -(VM*centr/V)/(KM+centr/V);
                cp = centr / V;
                cp ~ prop(prop.err)
            })
        }

        datr <- Bolus_1CPTMM
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)

        datr <- datr[datr$EVID != 2,]

        runno <- "N009"

        dat <- datr[datr$SD == 1,]

        fit <- uif %>% nlme(dat, control=nlmeControl(pnlsTol = 1, msVerbose = TRUE));

        z <- summary(as.nlme(fit))

        expect_equal(signif(as.numeric(as.nlme(fit)$logLik), 6), -12573.9)
        expect_equal(signif(AIC(as.nlme(fit)), 6), 25159.8)
        expect_equal(signif(BIC(as.nlme(fit)), 6), 25194.2)

        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[1]), 3), 6.81)
        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[2]), 3), 5.38)
        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[3]), 3), 4.18)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.328)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.297)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.2)

    })
}


