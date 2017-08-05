library(testthat)
library(nlmixr)
library(reshape2)

context("UI-NLME03: one-compartment bolus, steady-state")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {

    test_that("Closed-form", {

        datr <- Bolus_1CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2,]

        uif <- function(){
            ini({
                lCl <- 1.6
                lV <- 4.5
                prop.err <- 0.1
                eta.V ~ 0.1
                eta.Cl ~ 0.1
            })
            model({
                v <- exp(lV + eta.V)
                cl <- exp(lCl + eta.Cl)
                linCmt() ~ prop(prop.err)
            })
        }

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

        fit <- uif %>% nlme(dat);

        z <- summary(as.nlme(fit))

        expect_equal(signif(as.numeric(as.nlme(fit)$logLik),6), -12854.2)
        expect_equal(signif(AIC(as.nlme(fit)), 6), 25718.4)
        expect_equal(signif(BIC(as.nlme(fit)), 6), 25747)

        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[1]),3), 1.36)
        expect_equal(signif(as.numeric(as.nlme(fit)$coefficients$fixed[2]),3), 4.19)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.270)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.319)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.200)

    })

    test_that("ODE", {

        datr <- Bolus_1CPT
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2, ]

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

        uif <- function(){
            ini({
                lCl <- 1.6
                lV <- 4.5
                prop.err <- 0.1
                eta.V ~ 0.1
                eta.cl ~ 0.1
            })
            model({
                v <- exp(lV+ eta.V)
                cl <- exp(lCl + eta.cl)
                d / dt(centr) = -(cl / v) * centr;
                cp = centr / v;
                cp ~ prop(prop.err)
            })
        }

        fitODE <- uif %>% nlme(dat, control = nlmeControl(pnlsTol = .01, msVerbose = TRUE))

        z <- summary(as.nlme(fitODE))

        expect_equal(signif(as.numeric(as.nlme(fitODE)$logLik), 6),-12854.2)
        expect_equal(signif(AIC(as.nlme(fitODE)), 6), 25718.4)
        expect_equal(signif(BIC(as.nlme(fitODE)), 6), 25747)

        expect_equal(signif(as.numeric(as.nlme(fitODE)$coefficients$fixed[1]), 3), 1.36)
        expect_equal(signif(as.numeric(as.nlme(fitODE)$coefficients$fixed[2]), 3), 4.19)

        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.270)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.319)
        expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.200)


    })
}
