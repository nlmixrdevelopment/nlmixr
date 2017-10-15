library(testthat)
library(nlmixr)
library(reshape2)

rxPermissive({
    context("NLME03: one-compartment bolus, steady-state")

    test_that("Closed-form", {
        datr <- Bolus_1CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        specs1 <-
            list(
                fixed = lCL + lV ~ 1,
                random = pdDiag(lCL + lV ~ 1),
                start = c(lCL = 1.6, lV = 4.5)
            )

        runno <- "N003"

        datSS <- datr[datr$SS == 0,]
        datSD <- datr[datr$SS == 1,]

                                        # general solution to allow different times of SS dose and different II values per subject:
        datSSD <- datr[datr$SS == 1, c("ID", "TIME", "II")]

        nvar <- function(TIME, II) {
            for (i in seq(1, 7)) {
                r = TIME - (i - 1) * II
                assign(paste("ret", i, sep = ""), r)
            }
            return(list(
                r1 = ret1,
                r2 = ret2,
                r3 = ret3,
                r4 = ret4,
                r5 = ret5,
                r6 = ret6,
                r7 = ret7
            ))
        }

                                        #updates datSSD with 7 columns to account for the new dosing times
        datSSD[,(paste("V",seq(1,7)-1,sep=""))] <- nvar(datSSD$TIME, datSSD$II)
        datSSD <- datSSD[,c(1,4:10)]

        index <- melt(datSSD, id.vars="ID",value.name="TIMED")
        index <- index[,-2]
        index <- index[index$TIMED>0,]

                                        #much easier solution if you know the time of SS dose and the II and if it is the same for all
                                        #index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))

        datSD2 <- merge(datSD,index,by=c("ID"),all=TRUE)
        datSD2$TIME <- datSD2$TIMED
        datSD2 <- datSD2[,-15]

        datSS <- rbind(datSS,datSD2)
        datSS <- datSS[order(datSS$ID, datSS$TIME),]

        dat <- datSS

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs1,
                ncmt = 1,
                verbose = TRUE,
                oral = FALSE,
                weight = varPower(fixed = c(1))
            )

        z <- summary(fit)

        expect_equal(signif(as.numeric(fit$logLik),6), -12854.2)
        expect_equal(signif(AIC(fit), 6), 25718.4)
        expect_equal(signif(BIC(fit), 6), 25747)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.36)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.19)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.270)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.319)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.200)
    })

    test_that("ODE", {

        datr <- Bolus_1CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2, ]

        specs1 <-
            list(
                fixed = lCL + lV ~ 1,
                random = pdDiag(lCL + lV ~ 1),
                start = c(lCL = 1.6, lV = 4.5)
            )

        runno <- "N003"

        datSS <- datr[datr$SS == 0,]
        datSD <- datr[datr$SS == 1,]

                                        # general solution to allow different times of SS dose and different II values per subject:
        datSSD <- datr[datr$SS == 1, c("ID", "TIME", "II")]

        nvar <- function(TIME, II) {
            for (i in seq(1, 7)) {
                r = TIME - (i - 1) * II
                assign(paste("ret", i, sep = ""), r)
            }
            return(list(
                r1 = ret1,
                r2 = ret2,
                r3 = ret3,
                r4 = ret4,
                r5 = ret5,
                r6 = ret6,
                r7 = ret7
            ))
        }

                                        #updates datSSD with 7 columns to account for the new dosing times
        datSSD[,(paste("V",seq(1,7)-1,sep=""))] <- nvar(datSSD$TIME, datSSD$II)
        datSSD <- datSSD[,c(1,4:10)]

        index <- melt(datSSD, id.vars="ID",value.name="TIMED")
        index <- index[,-2]
        index <- index[index$TIMED>0,]

                                        #much easier solution if you know the time of SS dose and the II and if it is the same for all
                                        #index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))

        datSD2 <- merge(datSD,index,by=c("ID"),all=TRUE)
        datSD2$TIME <- datSD2$TIMED
        datSD2 <- datSD2[,-15]

        datSS <- rbind(datSS,datSD2)
        datSS <- datSS[order(datSS$ID, datSS$TIME),]

        dat <- datSS

        ode1 <- "
    d/dt(centr)  = -(CL/V)*centr;
    "

        mypar1 <- function(lCL, lV)
        {
            CL <- exp(lCL)
            V <- exp(lV)
        }

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

        expect_equal(signif(as.numeric(fitODE$logLik), 6),-12854.2)
        expect_equal(signif(AIC(fitODE), 6), 25718.4)
        expect_equal(signif(BIC(fitODE), 6), 25747)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]), 3), 1.36)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]), 3), 4.19)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.270)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.319)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.200)


    })
}, on.validate="NLMIXR_VALIDATION_FULL", silent=TRUE)
