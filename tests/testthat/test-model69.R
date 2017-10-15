library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME69: two-compartment oral Michaelis-Menten, multiple-dose")
    test_that("ODE", {

        datr <-
            read.csv("Oral_2CPTMM.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        ode2MMKA <- "
    d/dt(abs)    =-KA*abs;
    d/dt(centr)  = KA*abs+K21*periph-K12*centr-(VM*centr/V)/(KM+centr/V);
    d/dt(periph) =-K21*periph+K12*centr;
    "

        mypar9 <- function(lVM, lKM, lV, lCLD, lVT, lKA)
        {
            VM <- exp(lVM)
            KM <- exp(lKM)
            V <- exp(lV)
            CLD  <- exp(lCLD)
            VT <- exp(lVT)
            KA <- exp(lKA)
            K12 <- CLD / V
            K21 <- CLD / VT
        }

        specs9 <-
            list(
                fixed = lVM + lKM + lV + lCLD + lVT + lKA ~ 1,
                random = pdDiag(lVM + lKM + lV + lCLD + lVT + lKA ~ 1),
                start = c(
                    lVM = 7,
                    lKM = 6,
                    lV = 4,
                    lCLD = 1.5,
                    lVT = 4,
                    lKA = 0.1
                )
            )

        runno <- "N069"

        dat <- datr[datr$SD == 0,]

        fit <-
            nlme_ode(
                dat,
                model = ode2MMKA,
                par_model = specs9,
                par_trans = mypar9,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .1, msVerbose = TRUE)
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik), 6),-29758.9)
        expect_equal(signif(AIC(fit), 6), 59543.7)
        expect_equal(signif(BIC(fit), 6), 59627.6)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]), 3), 7.36)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]), 3), 6.51)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]), 3), 4.27)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]), 3), 1.21)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[5]), 3), 3.57)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[6]), 3), -0.0312)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.192)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.504)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.317)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.000534)
        expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.347)
        expect_equal(signif(as.numeric(z[6, "StdDev"]), 3), 0.000237)

        expect_equal(signif(fit$sigma, 3), 0.205)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
