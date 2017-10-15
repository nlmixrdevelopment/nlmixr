library(testthat)
library(nlmixr)
library(reshape2)
rxPermissive({
    context("NLME14: one-compartment infusion, steady state")
    test_that("Closed-form", {

        datr <-
            read.csv("Infusion_1CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        datSSobs <- datr[datr$SS == 0,]
        datSD <- datr[datr$SS == 1,]

                                        #general solution to allow different times of SS dose and different II values per subject:
        datSSD <- datSD[, c("ID","TIME","II")]

        specs1 <-
            list(
                fixed = lCL + lV ~ 1,
                random = pdDiag(lCL + lV ~ 1),
                start = c(lCL = 1.5, lV = 4)
            )

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

        runno <- "N014"

        dat <- datSS

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs1,
                ncmt = 1,
                verbose = TRUE,
                oral = FALSE,
                infusion = TRUE,
                weight = varPower(fixed = c(1))
            )

        z <- summary(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -12668.1)
        expect_equal(signif(AIC(fit), 6), 25346.2)
        expect_equal(signif(BIC(fit), 6), 25374.9)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.26)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.276)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.298)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.196)
    })
    test_that("ODE", {

        datr <-
            read.csv("Infusion_1CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        datSSobs <- datr[datr$SS == 0,]
        datSD <- datr[datr$SS == 1,]

                                        #general solution to allow different times of SS dose and different II values per subject:
        datSSD <- datSD[, c("ID","TIME","II")]

        specs1 <-
            list(
                fixed = lCL + lV ~ 1,
                random = pdDiag(lCL + lV ~ 1),
                start = c(lCL = 1.4, lV = 4.3)
            )

        ode1 <- "
    d/dt(centr)  = -(CL/V)*centr;
    "

        mypar1 <- function(lCL, lV)
        {
            CL <- exp(lCL)
            V <- exp(lV)
        }

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

        runno <- "N014"

        dat <- datSS

        fitODE <-
            nlme_ode(
                dat,
                model = ode1,
                par_model = specs1,
                par_trans = mypar1,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .01, msVerbose = TRUE)
            )

        z <- summary(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik),6), -12668.1)
        expect_equal(signif(AIC(fitODE), 6), 25346.2)
        expect_equal(signif(BIC(fitODE), 6), 25374.9)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.39)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.26)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.276)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.298)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.196)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")

