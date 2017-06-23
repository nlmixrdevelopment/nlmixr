library(testthat)
library(nlmixr)
library(data.table)

context("NLME: one-compartment oral, Michaelis-Menten, multiple-dose")

if (identical(Sys.getenv("NLMIXR_VALIDATION_FULL"), "true")) {
  
  test_that("ODE", {
    
    ode1MMKA <- "
    d/dt(abs)    = -KA*abs;
    d/dt(centr)  =  KA*abs-(VM*centr/V)/(KM+centr/V);
    "
    
    mypar5 <- function(lVM, lKM, lV, lKA)
    {
      VM <- exp(lVM)
      KM <- exp(lKM)
      V <- exp(lV)
      KA <- exp(lKA)
      
    }
    specs5 <-
      list(
        fixed = lVM + lKM + lV + lKA ~ 1,
        random = pdDiag(lVM + lKM + lV + lKA ~ 1),
        start = c(
          lVM = 7,
          lKM = 6,
          lV = 4,
          lKA = 0
        )
      )
    
    datr <-
      read.csv("ORAL_1CPTMM.csv",
               header = TRUE,
               stringsAsFactors = F)
    datr$EVID <- ifelse(datr$EVID == 1, 101, datr$EVID)
    datr <- data.table(datr)
    datr <- datr[EVID != 2]
    
    runno <- "N029"
    dat <- datr[SD == 1]
    
    fit <-
      nlme_ode(
        dat,
        model = ode1MMKA,
        par_model = specs5,
        par_trans = mypar5,
        response = "centr",
        response.scaler = "V",
        verbose = TRUE,
        weight = varPower(fixed = c(1)),
        control = nlmeControl(pnlsTol = .3, msVerbose = TRUE)
      )
    
    z <- VarCorr(fit)
    
    expect_equal(signif(as.numeric(fit$logLik),6), -12042.5)
    expect_equal(signif(AIC(fit), 6), 24102.9)
    expect_equal(signif(BIC(fit), 6), 24154.5)  
    
    expect_equal(signif(as.numeric(fit$coefficients$fixed[1]),3), 6.87)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[2]),3), 5.44)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[3]),3), 4.23)
    expect_equal(signif(as.numeric(fit$coefficients$fixed[4]),3), -0.00591)
    
    expect_equal(signif(as.numeric(z[1, "StdDev"]), 3), 0.267)
    expect_equal(signif(as.numeric(z[2, "StdDev"]), 3), 0.323)
    expect_equal(signif(as.numeric(z[3, "StdDev"]), 3), 0.291)
    expect_equal(signif(as.numeric(z[4, "StdDev"]), 3), 0.301)
    
    expect_equal(signif(fit$sigma, 3), 0.194)
  })
  
}