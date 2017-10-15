library(testthat)
library(nlmixr)
rxPermissive({
    context("NLME34: two-compartment bolus, steady-state")
    test_that("Closed-form", {
        datr <-Bolus_2CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2, ]


        ode2 <- "
    d/dt(centr)  = K21*periph-K12*centr-K10*centr;
    d/dt(periph) =-K21*periph+K12*centr;
    "

        mypar6 <- function(lCL, lV, lCLD, lVT)
        {
            CL <- exp(lCL)
            V  <- exp(lV)
            CLD <- exp(lCLD)
            VT <- exp(lVT)
            K10 <- CL / V
            K12 <- CLD / V
            K21 <- CLD / VT
        }

        specs6 <-
            list(
                fixed = lCL + lV + lCLD + lVT ~ 1,
                random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
                start = c(
                    lCL = 1.37,
                    lV = 4.19,
                    lCLD = 1.37,
                    lVT = 3.87
                )
            )

        runno <- "N034"

        datSS <- datr[datr$SS == 0, ]
        datSD <- datr[datr$SS == 1, ]

                                        #general solution to allow different times of SS dose and different II values per subject:
        datSSD <- datr[datr$SS == 1, c("ID", "TIME", "II")]

        datSSD$V0 <- datSSD$TIME
        datSSD$V1 <- datSSD$TIME - datSSD$II
        datSSD$V2 <- datSSD$TIME - 2 * datSSD$II
        datSSD$V3 <- datSSD$TIME - 3 * datSSD$II
        datSSD$V4 <- datSSD$TIME - 4 * datSSD$II
        datSSD$V5 <- datSSD$TIME - 5 * datSSD$II
        datSSD$V6 <- datSSD$TIME - 6 * datSSD$II
        datSSD$TIME <- NULL
        datSSD$II <- NULL

        index <- melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
        index$variable <- NULL
        index <- index[index$TIMED > 0, ]
        index <- index[order(index$ID, index$TIMED), ]

                                        #much easier solution if you know the time of SS dose and the II and if it is the same for all
                                        #index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))

        datSD2 <- merge(datSD, index, by = c("ID"), all = T)
        datSD2$TIME <- datSD2$TIMED
        datSD2$TIMED <- NULL

        datSS <- rbind(datSS, datSD2)
        datSS <- datSS[order(datSS$ID, datSS$TIME), ]
        dat <- datSS

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs6,
                ncmt = 2,
                verbose = TRUE,
                oral = FALSE,
                weight = varPower(fixed = c(1))
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik), 6),-13468.5)
        expect_equal(signif(AIC(fit), 6), 26955)
        expect_equal(signif(BIC(fit), 6), 27006.5)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]), 3), 1.34)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]), 3), 4.18)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]), 3), 1.31)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]), 3), 3.89)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.340)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.319)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.000775)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.292)

        expect_equal(signif(fit$sigma, 3), 0.2)
    })
    test_that("ODE", {
        datr <- Bolus_2CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2, ]

        ode2 <- "
    d/dt(centr)  = K21*periph-K12*centr-K10*centr;
    d/dt(periph) =-K21*periph+K12*centr;
    "

        mypar6 <- function(lCL, lV, lCLD, lVT)
        {
            CL <- exp(lCL)
            V  <- exp(lV)
            CLD <- exp(lCLD)
            VT <- exp(lVT)
            K10 <- CL / V
            K12 <- CLD / V
            K21 <- CLD / VT
        }

        specs6 <-
            list(
                fixed = lCL + lV + lCLD + lVT ~ 1,
                random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
                start = c(
                    lCL = 1.3,
                    lV = 4.19,
                    lCLD = 1.5,
                    lVT = 3.9
                )
            )

        runno <- "N034"

        datSS <- datr[datr$SS == 0, ]
        datSD <- datr[datr$SS == 1, ]

                                        #general solution to allow different times of SS dose and different II values per subject:
        datSSD <- datr[datr$SS == 1, c("ID", "TIME", "II")]

        datSSD$V0 <- datSSD$TIME
        datSSD$V1 <- datSSD$TIME - datSSD$II
        datSSD$V2 <- datSSD$TIME - 2 * datSSD$II
        datSSD$V3 <- datSSD$TIME - 3 * datSSD$II
        datSSD$V4 <- datSSD$TIME - 4 * datSSD$II
        datSSD$V5 <- datSSD$TIME - 5 * datSSD$II
        datSSD$V6 <- datSSD$TIME - 6 * datSSD$II
        datSSD$TIME <- NULL
        datSSD$II <- NULL

        index <- melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
        index$variable <- NULL
        index <- index[index$TIMED > 0, ]
        index <- index[order(index$ID, index$TIMED), ]

                                        #much easier solution if you know the time of SS dose and the II and if it is the same for all
                                        #index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))

        datSD2 <- merge(datSD, index, by = c("ID"), all = T)
        datSD2$TIME <- datSD2$TIMED
        datSD2$TIMED <- NULL

        datSS <- rbind(datSS, datSD2)
        datSS <- datSS[order(datSS$ID, datSS$TIME), ]
        dat <- datSS

        fitODE <-
            nlme_ode(
                dat,
                model = ode2,
                par_model = specs6,
                par_trans = mypar6,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .3, msVerbose = TRUE)
            )

        z <- VarCorr(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik), 6),-13466.2)
        expect_equal(signif(AIC(fitODE), 6), 26950.3)
        expect_equal(signif(BIC(fitODE), 6), 27001.9)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]), 3), 1.34)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]), 3), 4.18)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]), 3), 1.29)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]), 3), 3.91)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.338)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.318)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.000764)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.280)

        expect_equal(signif(fitODE$sigma, 3), 0.201)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
