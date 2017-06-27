library(testthat)
library(nlmixr)

context("NLME41: two-compartment bolus Michaelis-Menten, multiple-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {
  
  test_that("ODE", {

    datr <-
      read.csv("Bolus_2CPTMM.csv",
               header = TRUE,
               stringsAsFactors = F)
    datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
    datr <- datr[datr$EVID != 2,]
    
    ode2MM <- "
    d/dt(centr)  = K21*periph-K12*centr-(VM*centr/V)/(KM+centr/V);
    d/dt(periph) =-K21*periph+K12*centr;
    "
    
    mypar7 <- function(lVM, lKM, lV, lCLD, lVT)
    {
      VM <- exp(lVM)
      KM <- exp(lKM)
      V <- exp(lV)
      CLD  <- exp(lCLD)
      VT <- exp(lVT)
      K12 <- CLD / V
      K21 <- CLD / VT
    }
    specs7 <-
      list(
        fixed = lVM + lKM + lV + lCLD + lVT ~ 1,
        random = pdDiag(lVM + lKM + lV + lCLD + lVT ~ 1),
        start = c(
          lVM = 7.6,
          lKM = 6.83,
          lV = 4.29,
          lCLD = 1.26,
          lVT = 3.58
        )
      )
    
    runno <- "N041"
    
    dat <- datr[datr$SD == 0,]
    
    fit <-
      nlme_ode(
        dat,
        model = ode2MM,
        par_model = specs7,
        par_trans = mypar7,
        response = "centr",
        response.scaler = "V",
        verbose = TRUE,
        weight = varPower(fixed = c(1)),
        control = nlmeControl(pnlsTol = .1, msVerbose = TRUE)
      )
    
    z <- VarCorr(fit)
    
    expect_equal(signif(as.numeric(fit$logLik), 6),-29938.3)
    expect_equal(signif(AIC(fit), 6), 59898.5)
    expect_equal(signif(BIC(fit), 6), 59969.5)
    
    expect_equal(signif(as.numeric(fit$coefficients$fixed[1]), 3), 7.61)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[2]), 3), 6.91)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[3]), 3), 4.26)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[4]), 3), 1.41)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[5]), 3), 3.68)
    
    expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.00974)
    expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.727)
    expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.292)
    expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.00144)
    expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.417)
    
    expect_equal(signif(fit$sigma, 3), 0.211)
  })
  
}