library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME21: one-compartment infusion, multiple-dose, Michaelis-Menten")
    test_that("ODE", {

        datr <-
            read.csv("Infusion_1CPTMM.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)

        datr <- subset(datr, EVID != 2)
        datIV <- subset(datr, AMT>0)
        datIV$TIME <- datIV$TIME + (datIV$AMT/datIV$RATE)
        datIV$AMT  <- -1*datIV$AMT
                                        #datIV <- datr[AMT > 0][, TIME := TIME + AMT / RATE][, AMT := -1 * AMT]
        datr <- rbind(datr, datIV)
        datr <- datr[order(datr$ID, datr$TIME),]

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

    runno <- "N021"

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

    z <- VarCorr(fit)

    expect_equal(signif(as.numeric(fit$logLik),6), -33250.6)
    expect_equal(signif(AIC(fit), 6), 66515.2)
    expect_equal(signif(BIC(fit), 6), 66560.4)

    expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 7.08)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 5.83)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 4.16)

    expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.235)
    expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.466)
    expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.288)

    expect_equal(signif(fit$sigma, 3), 0.211)
})
}, on.validate="NLMIXR_VALIDATION_FULL")
