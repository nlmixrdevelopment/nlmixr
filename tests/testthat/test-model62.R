library(testthat)
library(nlmixr)
library(data.table)

context("NLME: two-compartment oral, steady-state, multiple-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {
  
  test_that("Closed-form", {
    
    datr <- read.csv("ORAL_2CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)
    datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
    datr <- data.table(datr)
    datr <- datr[EVID != 2]
    
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
    
    expect_equal(signif(as.numeric(fit$logLik),6), -13224.6)
    expect_equal(signif(AIC(fit), 6), 26471.3)
    expect_equal(signif(BIC(fit), 6), 26534.3)  
    
    expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 1.34)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 4.28)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 1.33)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), 3.76)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[5]),3), -0.0101)
    
    expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.321)
    expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.297)
    expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.452)
    expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.0013)
    expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.333)
    
    expect_equal(signif(fit$sigma, 3), 0.199)
  })
  
  test_that("ODE", {
    
    datr <- read.csv("ORAL_2CPT.csv",
                     header = TRUE,
                     stringsAsFactors = F)
    datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
    datr <- data.table(datr)
    datr <- datr[EVID != 2]
    
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
    
    expect_equal(signif(as.numeric(fitODE$logLik),6), -13224.4)
    expect_equal(signif(AIC(fitODE), 6), 26470.9)
    expect_equal(signif(BIC(fitODE), 6), 26533.9)  
    
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[1]),3), 1.34)
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[2]),3), 4.28)
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[3]),3), 1.33)
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[4]),3), 3.76)
    expect_equal(signif(as.numeric(fitODE$coefficients$fixed[5]),3), -0.00741)
    
    expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.321)
    expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.297)
    expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.45)
    expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.00125)
    expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.333)
    
    expect_equal(signif(fit$sigma, 3), 0.199)
  })
  
}