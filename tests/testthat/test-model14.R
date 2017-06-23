library(testthat)
library(nlmixr)
library(data.table)

context("NLME: one-compartment infusion, steady state")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {
  
  test_that("Closed-form", {
    
    datr <-
      read.csv("Infusion_1CPT.csv",
               header = TRUE,
               stringsAsFactors = F)
    datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)
    datr <- data.table(datr)
    datr <- datr[EVID != 2]
    datSSobs <- datr[SS == 0]
    datSD <- datr[SS == 1]
    setkey(datSD, ID, TIME)
    #general solution to allow different times of SS dose and different II values per subject:
    datSSD <- datSD[, .(ID, TIME, II)]
    
    specs1 <-
      list(
        fixed = lCL + lV ~ 1,
        random = pdDiag(lCL + lV ~ 1),
        start = c(lCL = 1.5, lV = 4)
      )
    
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
    datSSD[, (paste("V", seq(1, 7) - 1, sep = "")) := nvar(TIME, II)][, TIME :=
                                                                        NULL][, II := NULL]
    
    index <- melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
    index[, variable := NULL]
    index <- index[TIMED > 0]
    setkey(index, ID, TIMED)
    
    #much easier solution if you know the time of SS dose and the II and if it is the same for all
    #index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))
    
    datSD2 <- merge(datSD, index, by = c("ID"), all = TRUE)
    datSD2[, TIME := TIMED][, TIMED := NULL]
    
    datSDoff <- copy(datSD2)
    datSDoff[, TIME := TIME + AMT / RATE][, AMT := -1 * AMT]
    datSD2 <- rbind(datSD2, datSDoff)
    
    datSS <- rbind(datSSobs, datSD2)
    setkey(datSS, ID, TIME)
    
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
    datr <- data.table(datr)
    datr <- datr[EVID != 2]
    datSSobs <- datr[SS == 0]
    datSD <- datr[SS == 1]
    setkey(datSD, ID, TIME)
    #general solution to allow different times of SS dose and different II values per subject:
    datSSD <- datSD[, .(ID, TIME, II)]
    
    specs1 <-
      list(
        fixed = lCL + lV ~ 1,
        random = pdDiag(lCL + lV ~ 1),
        start = c(lCL = 1.6, lV = 4.5)
      )
    
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
    datSSD[, (paste("V", seq(1, 7) - 1, sep = "")) := nvar(TIME, II)][, TIME :=
                                                                        NULL][, II := NULL]
    
    index <- melt(datSSD, id.vars = c("ID"), value.name = "TIMED")
    index[, variable := NULL]
    index <- index[TIMED > 0]
    setkey(index, ID, TIMED)
    
    #much easier solution if you know the time of SS dose and the II and if it is the same for all
    #index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))
    
    datSD2 <- merge(datSD, index, by = c("ID"), all = TRUE)
    datSD2[, TIME := TIMED][, TIMED := NULL]
    
    datSDoff <- copy(datSD2)
    datSDoff[, TIME := TIME + AMT / RATE][, AMT := -1 * AMT]
    datSD2 <- rbind(datSD2, datSDoff)
    
    datSS <- rbind(datSSobs, datSD2)
    setkey(datSS, ID, TIME)
    
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
    expect_equal(signif(AIC(fitODE), 6), 25346.1)
    expect_equal(signif(BIC(fitODE), 6), 25374.8)  
    
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.39)
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.26)
    
    expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[1], 3)), 0.276)
    expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[2], 3)), 0.298)
    expect_equal(as.numeric(signif(exp(attr(z$apVar, "Pars"))[3], 3)), 0.196)
  }) 

}