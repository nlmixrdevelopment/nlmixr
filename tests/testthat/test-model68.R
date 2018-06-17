library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME68: two-compartment oral Michaelis-Menten, single-dose")
    test_that("ODE", {

        datr <-
            read.csv("Oral_2CPTMM.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

       ode2MMKA <- "
    d/dt(abs)    =-KA*abs;
    d/dt(centr)  = KA*abs+K21*periph-K12*centr-exp(log(VM)+log(centr)-log(V)-log(KM+centr/V));
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
                    lVM = 7.3,
                    lKM = 6.1,
                    lV = 4.3,
                    lCLD = 1.4,
                    lVT = 3.75,
                    lKA = -0.01
                )
            )

        runno <- "N068"

        dat <- datr[datr$SD == 1,]

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
                control = nlmeControl(pnlsTol = .1, msVerbose = TRUE, msMaxIter=1000)
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik), 6), -11769.1)
        expect_equal(signif(AIC(fit), 6), 23564.2)
        expect_equal(signif(BIC(fit), 6), 23638.7)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]), 3), 7.32)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]), 3), 6.09)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]), 3), 4.25)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]), 3), 1.38)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[5]), 3), 3.75)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[6]), 3), -0.021)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3),  0.19)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.33)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.336)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.000205)
        expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.253)
        expect_equal(signif(as.numeric(z[6, "StdDev"]), 3), 0.283)

        expect_equal(signif(fit$sigma, 3), 0.205)
    })
})
