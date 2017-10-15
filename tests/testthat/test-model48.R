library(testthat)
library(nlmixr)
library(reshape2)

rxPermissive({
    context("NLME48: two-compartment infusion, steady-state")
    test_that("Closed-form", {

        datr <-
            read.csv("Infusion_2CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        datSSobs <- datr[datr$SS == 0,]
        datSD <- datr[datr$SS == 1,]

                                        #general solution to allow different times of SS dose and different II values per subject:
        datSSD <- datSD[, c("ID","TIME","II")]

                                        #updates datSSD with 7 columns to account for the new dosing times
        datSSD$V0<-datSSD$TIME
        datSSD$V1<-datSSD$TIME-datSSD$II
        datSSD$V2<-datSSD$TIME-2*datSSD$II
        datSSD$V3<-datSSD$TIME-3*datSSD$II
        datSSD$V4<-datSSD$TIME-4*datSSD$II
        datSSD$V5<-datSSD$TIME-5*datSSD$II
        datSSD$V6<-datSSD$TIME-6*datSSD$II
        datSSD$TIME<-NULL
        datSSD$II<-NULL

        index <- melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
        index$variable <- NULL
        index <- index[index$TIMED > 0,]
        index<-index[order(index$ID,index$TIMED),]

        datSD2 <- merge(datSD, index, by = c("ID"), all = TRUE)
        datSD2$TIME <- datSD2$TIMED
        datSD2$TIMED <- NULL

        datSDoff <- datSD2
        datSDoff$TIME <- datSDoff$TIME + datSDoff$AMT / datSDoff$RATE
        datSDoff$AMT <- -1 * datSDoff$AMT

        datSD2 <- rbind(datSD2, datSDoff)

        datSS <- rbind(datSSobs, datSD2)
        datSS <- datSS[order(datSS$ID,datSS$TIME),]

        dat <- datSS

        specs6 <-
            list(
                fixed = lCL + lV + lCLD + lVT ~ 1,
                random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
                start = c(
                    lCL = 1.36,
                    lV = 4.2,
                    lCLD = 1.47,
                    lVT = 3.9
                )
            )

        runno <- "N048"

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs6,
                ncmt = 2,
                verbose = TRUE,
                oral = FALSE,
                infusion = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(
                    pnlsTol = .1,
                    msVerbose = TRUE,
                    maxIter = 200
                )
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -13384.3)
        expect_equal(signif(AIC(fit), 6), 26786.5)
        expect_equal(signif(BIC(fit), 6), 26838.1)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.36)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.21)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 1.34)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), 3.91)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.294)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.285)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.00177)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.363)

        expect_equal(signif(fit$sigma, 3), 0.207)
    })
    test_that("ODE", {

        datr <-
            read.csv("Infusion_2CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        datSSobs <- datr[datr$SS == 0,]
        datSD <- datr[datr$SS == 1,]

                                        #general solution to allow different times of SS dose and different II values per subject:
        datSSD <- datSD[, c("ID","TIME","II")]

                                        #updates datSSD with 7 columns to account for the new dosing times
        datSSD$V0<-datSSD$TIME
        datSSD$V1<-datSSD$TIME-datSSD$II
        datSSD$V2<-datSSD$TIME-2*datSSD$II
        datSSD$V3<-datSSD$TIME-3*datSSD$II
        datSSD$V4<-datSSD$TIME-4*datSSD$II
        datSSD$V5<-datSSD$TIME-5*datSSD$II
        datSSD$V6<-datSSD$TIME-6*datSSD$II
        datSSD$TIME<-NULL
        datSSD$II<-NULL

        index <- melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
        index$variable <- NULL
        index <- index[index$TIMED > 0,]
        index<-index[order(index$ID,index$TIMED),]

        datSD2 <- merge(datSD, index, by = c("ID"), all = TRUE)
        datSD2$TIME <- datSD2$TIMED
        datSD2$TIMED <- NULL

        datSDoff <- datSD2
        datSDoff$TIME <- datSDoff$TIME + datSDoff$AMT / datSDoff$RATE
        datSDoff$AMT <- -1 * datSDoff$AMT

        datSD2 <- rbind(datSD2, datSDoff)

        datSS <- rbind(datSSobs, datSD2)
        datSS <- datSS[order(datSS$ID,datSS$TIME),]

        dat <- datSS

        specs6 <-
            list(
                fixed = lCL + lV + lCLD + lVT ~ 1,
                random = pdDiag(lCL + lV + lCLD + lVT ~ 1),
                start = c(
                    lCL = 1.36,
                    lV = 4.21,
                    lCLD = 1.34,
                    lVT = 3.91
                )
            )

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

        runno <- "N048"

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
                control = nlmeControl(pnlsTol = .1, msVerbose = TRUE)
            )

        z <- VarCorr(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik),6), -13384.2)
        expect_equal(signif(AIC(fitODE), 6), 26786.4)
        expect_equal(signif(BIC(fitODE), 6), 26838)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.36)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.21)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), 1.34)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]),3), 3.91)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.294)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.285)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.0031)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.363)

        expect_equal(signif(fitODE$sigma, 3), 0.207)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
