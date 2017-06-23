library(testthat)
library(nlmixr)
library(data.table)

context("NLME: one-compartment oral, multiple-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {
  
  test_that("Closed-form", {
    
    datr <-
      read.csv("ORAL_1CPT.csv",
               header = TRUE,
               stringsAsFactors = F)
    datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
    datr <- data.table(datr)
    datr <- datr[EVID != 2]
    
    ode1KA <- "
    d/dt(abs)    = -KA*abs;
    d/dt(centr)  =  KA*abs-(CL/V)*centr;
    "
    
    mypar4 <- function(lCL, lV, lKA)
    {
      CL <- exp(lCL)
      V <- exp(lV)
      KA <- exp(lKA)
    }
    
    specs4 <-
      list(
        fixed = lCL + lV + lKA ~ 1,
        random = pdDiag(lCL + lV + lKA ~ 1),
        start = c(lCL = 1, lV = 4, lKA = 0)
      )
    
    runno<-"N025"
    
    datSS <- datr[SS == 0]
    datSD <- datr[SS == 1]
    
    #general solution to allow different times of SS dose and different II values per subject:
    datSSD <- datr[SS == 1, .(ID, TIME, II)]
    
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
    index <- melt(datSSD, id.vars = "ID", value.name = "TIMED")
    index[, variable := NULL]
    index <- index[TIMED > 0]
    
    #much easier solution if you know the time of SS dose and the II and if it is the same for all
    #index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))
    
    datSD2 <- merge(datSD, index, by = c("ID"), all = TRUE)
    datSD2[, TIME := TIMED][, TIMED := NULL]
    datSS <- rbind(datSS, datSD2)
    setkey(datSS, ID, TIME)
    dat <- datSS
    
    fit <-
      nlme_lin_cmpt(
        dat,
        par_model = specs4,
        ncmt = 1,
        verbose = TRUE,
        oral = TRUE,
        weight = varPower(fixed = c(1))
      )
    
    z <- VarCorr(fit)
    
    expect_equal(signif(as.numeric(fit$logLik),6), -12633.5)
    expect_equal(signif(AIC(fit), 6), 25281)
    expect_equal(signif(BIC(fit), 6), 25321.2)  
    
    expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.40)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.22)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 0.00131)
    
    expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.271)
    expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.284)
    expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.00042)
    
    expect_equal(signif(fit$sigma, 3), 0.209)
  })
  
  test_that("ODE", {
    
    datr <-
      read.csv("ORAL_1CPT.csv",
               header = TRUE,
               stringsAsFactors = F)
    datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
    datr <- data.table(datr)
    datr <- datr[EVID != 2]
    
    ode1KA <- "
    d/dt(abs)    = -KA*abs;
    d/dt(centr)  =  KA*abs-(CL/V)*centr;
    "
    
    mypar4 <- function(lCL, lV, lKA)
    {
      CL <- exp(lCL)
      V <- exp(lV)
      KA <- exp(lKA)
    }
    
    specs4 <-
      list(
        fixed = lCL + lV + lKA ~ 1,
        random = pdDiag(lCL + lV + lKA ~ 1),
        start = c(lCL = 1, lV = 4, lKA = 0)
      )
    
    runno<-"N025"
    
    datSS <- datr[SS == 0]
    datSD <- datr[SS == 1]
    
    #general solution to allow different times of SS dose and different II values per subject:
    datSSD <- datr[SS == 1, .(ID, TIME, II)]
    
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
    index <- melt(datSSD, id.vars = "ID", value.name = "TIMED")
    index[, variable := NULL]
    index <- index[TIMED > 0]
    
    #much easier solution if you know the time of SS dose and the II and if it is the same for all
    #index<-CJ(ID=datSSD$ID,TIMED=seq(192,0,-24))
    
    datSD2 <- merge(datSD, index, by = c("ID"), all = TRUE)
    datSD2[, TIME := TIMED][, TIMED := NULL]
    datSS <- rbind(datSS, datSD2)
    setkey(datSS, ID, TIME)
    dat <- datSS
    
    fitODE <-
      nlme_ode(
        dat,
        model = ode1KA,
        par_model = specs4,
        par_trans = mypar4,
        response = "centr",
        response.scaler = "V",
        verbose = TRUE,
        weight = varPower(fixed = c(1)),
        control = nlmeControl(pnlsTol = .1, msVerbose = TRUE)
      )
    
    z <- VarCorr(fitODE)
    
    expect_equal(signif(as.numeric(fitODE$logLik),6), -12633.4)
    expect_equal(signif(AIC(fitODE), 6), 25280.8)
    expect_equal(signif(BIC(fitODE), 6), 25320.9)  
    
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.4)
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.22)
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), 0.00439)
    
    expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.271)
    expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.284)
    expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.000419)
    
    expect_equal(signif(fitODE$sigma, 3), 0.209)
  })
  
}