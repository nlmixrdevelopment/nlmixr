nlmixrTest({

  context("full cov FOCEi (posthoc surrogate)")

  library(dplyr)

  warfarin %>%
    filter(dvid == "cp") ->
  dat

  #### doesn't work with FOCEI ### doesn't work with SAEM with CRWES=True but does iwth CWRES=FALSE
  One.SD.ODE <- function() {
    ini({
      # Where initial conditions/variables are specified
      lcl <- log(0.135) # log Cl (L/h)
      lv <- log(8) # log V (L)
      lmtt <- log(1.1) # log MTx

      prop.err <- 0.15 # proportional error (SD/mean)
      add.err <- 0.6 # additive error (mg/L)
      eta.cl + eta.v + eta.mtt ~ c(
        0.1,
        0.001, 0.1,
        0.001, 0.001, 0.1
      )
    })
    model({
      # Where the model is specified
      cl <- exp(lcl + eta.cl)
      v <- exp(lv + eta.v)
      mtt <- exp(lmtt + eta.mtt)
      ktr <- 6 / mtt

      ## ODE example
      d / dt(depot) <- -ktr * depot
      d / dt(central) <- ktr * trans5 - (cl / v) * central
      d / dt(trans1) <- ktr * depot - ktr * trans1
      d / dt(trans2) <- ktr * trans1 - ktr * trans2
      d / dt(trans3) <- ktr * trans2 - ktr * trans3
      d / dt(trans4) <- ktr * trans3 - ktr * trans4
      d / dt(trans5) <- ktr * trans4 - ktr * trans5

      ## Concentration is calculated
      cp <- central / v
      ## And is assumed to follow proportional and additive error
      cp ~ prop(prop.err) + add(add.err)
    })
  }


  f <- nlmixr(One.SD.ODE, dat, "posthoc")
  expect_true(inherits(f, "nlmixrPosthoc"))

}, test="focei")
