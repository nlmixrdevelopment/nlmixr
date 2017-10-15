library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME30: one-compartment oral, Michaelis-Menten, multiple-dose")
    test_that("ODE", {

        ode1MMKA <- "
    d/dt(abs)    = -KA*abs;
    d/dt(centr)  =  KA*abs-(VM*centr/V)/(KM+centr/V);
    "

        mypar5 <- function(lVM, lKM, lV, lKA)
        {
            VM <- exp(lVM)
            KM <- exp(lKM)
            V <- exp(lV)
            KA <- exp(lKA)

        }
        specs5i <-
            list(
                fixed = lVM + lKM + lV + lKA ~ 1,
                random = pdDiag(value = diag(c(2, 2, 2, 2)), form = lVM + lKM + lV + lKA ~
                                                                 1),
                start = c(
                    lVM = 7,
                    lKM = 6.2,
                    lV = 4.5,
                    lKA = -0.2
                )
            )

        datr <-
            read.csv("Oral_1CPTMM.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        runno <- "N030"
        dat <- datr[datr$SD == 0,]

        fit <-
            nlme_ode(
                dat,
                model = ode1MMKA,
                par_model = specs5i,
                par_trans = mypar5,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .3, msVerbose = TRUE)
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -30897.9)
        expect_equal(signif(AIC(fit), 6), 61813.9)
        expect_equal(signif(BIC(fit), 6), 61871.9)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 7.09)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 6.11)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 4.18)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), -0.0836)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.158)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.728)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.299)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.000358)

        expect_equal(signif(fit$sigma, 3), 0.203)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
