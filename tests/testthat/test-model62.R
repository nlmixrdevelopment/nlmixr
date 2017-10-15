library(testthat)
library(nlmixr)

rxPermissive({
    context("NLME62: two-compartment oral, steady-state, multiple-dose")
    test_that("Closed-form", {
        datr <- read.csv("Oral_2CPT.csv",
                         header = TRUE,
                         stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2, ]

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

        specs8i <-
            list(
                fixed = lCL + lV + lCLD + lVT + lKA ~ 1,
                random = pdDiag(value = diag(c(6, 6, 6, 6, 6)), form = lCL + lV + lCLD +
                                                                    lVT + lKA ~ 1),
                start = c(
                    lCL = 1.4,
                    lV = 4.2,
                    lCLD = 1.3,
                    lVT = 3.9,
                    lKA = 0.1
                )
            )

        runno <- "N062"

        fit <-
            nlme_lin_cmpt(
                dat,
                par_model = specs8i,
                ncmt = 2,
                verbose = TRUE,
                oral = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .15, msVerbose = TRUE)
            )

        z <- VarCorr(fit)

        expect_equal(signif(as.numeric(fit$logLik), 6),-13243.6)
        expect_equal(signif(AIC(fit), 6), 26509.2)
        expect_equal(signif(BIC(fit), 6), 26572.3)

        expect_equal(signif(as.numeric(fit$coefficients$fixed[1]), 3), 1.35)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[2]), 3), 4.22)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[3]), 3), 1.48)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[4]), 3), 3.86)
        expect_equal(signif(as.numeric(fit$coefficients$fixed[5]), 3),-0.0586)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.322)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.301)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.00218)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.306)
        expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.00106)

        expect_equal(signif(fit$sigma, 3), 0.201)
    })
    test_that("ODE", {
        datr <- read.csv("Oral_2CPT.csv",
                         header = TRUE,
                         stringsAsFactors = F)
        datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
        datr <- datr[datr$EVID != 2, ]

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

        specs8i <-
            list(
                fixed = lCL + lV + lCLD + lVT + lKA ~ 1,
                random = pdDiag(value = diag(c(6, 6, 6, 6, 6)), form = lCL + lV + lCLD +
                                                                    lVT + lKA ~ 1),
                start = c(
                    lCL = 1.4,
                    lV = 4.2,
                    lCLD = 1.3,
                    lVT = 3.9,
                    lKA = 0.1
                )
            )

        runno <- "N062"

        ode2KA <- "
    d/dt(abs)    = -KA*abs;
    d/dt(centr)  =  KA*abs+K21*periph-K12*centr-K10*centr;
    d/dt(periph) =        -K21*periph+K12*centr;
    "

        mypar8 <- function(lCL, lV, lCLD, lVT, lKA)
        {
            CL <- exp(lCL)
            V  <- exp(lV)
            CLD <- exp(lCLD)
            VT <- exp(lVT)
            KA <- exp(lKA)
            K10 <- CL / V
            K12 <- CLD / V
            K21 <- CLD / VT
        }

        fitODE <-
            nlme_ode(
                dat,
                model = ode2KA,
                par_model = specs8i,
                par_trans = mypar8,
                response = "centr",
                response.scaler = "V",
                verbose = TRUE,
                weight = varPower(fixed = c(1)),
                control = nlmeControl(pnlsTol = .15, msVerbose = TRUE)
            )

        z <- VarCorr(fitODE)

        expect_equal(signif(as.numeric(fitODE$logLik), 6),-13243.4)
        expect_equal(signif(AIC(fitODE), 6), 26508.7)
        expect_equal(signif(BIC(fitODE), 6), 26571.8)

        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]), 3), 1.35)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]), 3), 4.23)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]), 3), 1.43)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]), 3), 3.85)
        expect_equal(signif(as.numeric(fitODE$coefficients$fixed[5]), 3),-0.0455)

        expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.322)
        expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.302)
        expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.00214)
        expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.303)
        expect_equal(round(as.numeric(z[5, "StdDev"]), 3), 0)

        expect_equal(signif(fitODE$sigma, 3), 0.201)
    })
}, on.validate="NLMIXR_VALIDATION_FULL")
