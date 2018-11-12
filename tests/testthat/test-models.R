library(testthat)
library(nlmixr)

rxPermissive({
  context("Test all models")
  test_that("Models", {
    files <- list.files(path="tests/testthat/models/", full.names=TRUE, pattern="^model*")
    for (current_file in files) {
      source(current_file)
      z <- VarCorr(fit[[runno]])
    
      expect_equal(
        round(c(fit[[runno]]$logLik, AIC(fit[[runno]]), BIC(fit[[runno]])), 2),
        expected_values$lik,
        info=paste("Likelihood for", runno)
      )

      expect_equal(
        unname(fit[[runno]]$coefficients$fixed),
        expected_values$param,
        tol=1e-3,
        info=paste("Parameters for", runno)
      )

      expect_equal(
        unname(z[-nrow(z), "StdDev"]),
        expected_values$stdev_param,
        tol=1e-3,
        info=paste("Parameter stdev for", runno)
      )
    
      expect_equal(
        fit[[runno]]$sigma,
        expected_values$sigma,
        tol=1e-3,
        info=paste("Sigma for", runno)
      )
    }
  })
}, on.validate="NLMIXR_VALIDATION_FULL",silent=TRUE)
