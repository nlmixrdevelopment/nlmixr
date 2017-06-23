library(testthat)
library(nlmixr)
library(data.table)

context("NLME: two-compartment infusion Michaelis-Menten, multiple-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {
  
  test_that("ODE", {

    datr <-
      read.csv("INFUSION_2CPTMM.csv",
               header = TRUE,
               stringsAsFactors = F)
    datr$EVID <- ifelse(datr$EVID == 1, 10101, datr$EVID)
    datr <- data.table(datr)
    datr <- datr[EVID != 2]
    datIV <- datr[AMT > 0][, TIME := TIME + AMT / RATE][, AMT := -1 * AMT]
    datr <- rbind(datr, datIV)
    setkey(datr, ID, TIME)
    
    runno <- "N055"
    
    dat <- datr[SD == 0]
    
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
    
    specs7i <-
      list(
        fixed = lVM + lKM + lV + lCLD + lVT ~ 1,
        random = pdDiag(value = diag(c(3, 3, 3, 3, 3)), form = lVM + lKM + lV +
                          lCLD + lVT ~ 1),
        start = c(
          lVM = 7,
          lKM = 5,
          lV = 4,
          lCLD = 1.2,
          lVT = 4
        )
      )
    
    fit <-
      nlme_ode(
        dat,
        model = ode2MM,
        par_model = specs7i,
        par_trans = mypar7,
        response = "centr",
        response.scaler = "V",
        verbose = TRUE,
        weight = varPower(fixed = c(1)),
        control = nlmeControl(pnlsTol = .1, msVerbose = TRUE)
      )
    
    z <- VarCorr(fit)
    
    expect_equal(signif(as.numeric(fit$logLik), 6),-31637.2)
    expect_equal(signif(AIC(fit), 6), 63296.3)
    expect_equal(signif(BIC(fit), 6), 63367.3)
    
    expect_equal(signif(as.numeric(fit$coefficients$fixed[1]), 3), 7.54)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[2]), 3), 7.14)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[3]), 3), 4.25)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[4]), 3), 0.661)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[5]), 3), 3.84)
    
    expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.000193)
    expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.904)
    expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.279)
    expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.47)
    expect_equal(signif(as.numeric(z[5, "StdDev"]), 3), 0.369)
    
    expect_equal(signif(fit$sigma, 3), 0.202)
  })
  
}