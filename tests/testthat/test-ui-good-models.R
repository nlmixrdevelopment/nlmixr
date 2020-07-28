nlmixrTest({
  context("Good UI models should not raise errors")

  one.compartment.saem <- function() {
    ini({
      tka <- .5 # Log Ka
      tcl <- -3.2 # Log Cl
      tv <- -1 # Log V
      eta.ka ~ 1
      eta.cl ~ 2
      eta.v ~ 1
      add.err <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot + exp(-k0 * t)
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err)
    })
  }

  expect_equal(nlmixr(one.compartment.saem)$all.covs, "k0")

  two.cmt.pd <- function() {
    ini({
      tKa <- log(0.64)
      tCl <- log(5.22)
      tV2 <- log(41.3)
      tV3 <- log(115)
      tQ <- log(11.96)
      BWef <- log(1.87)
      tSlope <- log(10) ## add for PD estimation
      tIntercept <- log(1) ## add for PD estimation
      eta.Ka ~ 1.18
      eta.Cl ~ 0.09
      eta.V2 ~ 0.2
      eta.V3 ~ 0.12
      eta.Q ~ 0.12
      eta.Slope ~ 0.1 ## add for PD estimation
      eta.Intercept ~ 0.1 ## add for PD estimation

      prop.err1 <- 0.1 ## Cp
      prop.err2 <- 0.3 ## Ef
    })
    model({
      Ka <- exp(tKa + eta.Ka)
      Cl <- exp(tCl + BWef * log.BW.70 + eta.Cl)
      V2 <- exp(tV2 + eta.V2)
      V3 <- exp(tV3 + eta.V3)
      Q <- exp(tQ + eta.Q)
      Slope <- exp(tSlope + eta.Slope) ## add for PD estimation
      Intercept <- exp(tIntercept + eta.Intercept) ## add for PD estimation

      d / dt(depot) <- -Ka * depot
      d / dt(center) <- Ka * depot - Cl / V2 * center + Q / V3 * periph - Q / V2 * center
      d / dt(periph) <- Q / V2 * center - Q / V3 * periph

      Cp <- center / V2
      Ef <- Cp * Slope + Intercept ## add for PD estimation

      Cp ~ prop(prop.err1) | center
      Ef ~ prop(prop.err2) ## add for PD estimation
    })
  }

  expect_equal("nlmixrUI", class(nlmixr(two.cmt.pd)))


  one.compartment.IV.model <- function() {
    ini({ # Where initial conditions/variables are specified
      # '<-' or '=' defines population parameters
      # Simple numeric expressions are supported
      Cl <- 1.6 # Cl (L/hr)
      Vc <- 4.5 # V (L)
      # Bounds may be specified by c(lower, est, upper), like NONMEM:
      # Residuals errors are assumed to be population parameters
      prop.err <- c(0, 0.3, 1)
      # Between subject variability estimates are specified by '~'
      # Semicolons are optional
      # eta.Vc ~ 0.1   #IIV V
      # eta.Cl ~ 0.1   #IIV Cl
    })
    model({ # Where the model is specified
      # The model uses the ini-defined variable names
      # Vc <- exp(lVc + eta.Vc)
      # Cl <- exp(lCl + eta.Cl)
      # RxODE-style differential equations are supported
      d / dt(centr) <- -(Cl / Vc) * centr
      ## Concentration is calculated
      cp <- centr / Vc
      # And is assumed to follow proportional error estimated by prop.err
      cp ~ prop(prop.err)
    })
  }

  expect_equal("nlmixrUI", class(nlmixr(one.compartment.IV.model)))


  model1 <- function() {
    ini({
      CL <- 2.2
      V <- 65
      add.err <- 0.01
      prop.err <- 0.01
    })
    model({
      kel <- CL / V
      X(0) <- 0
      d / dt(X) <- -kel * X
      cp <- X / V
      cp ~ add(add.err) + prop(prop.err)
    })
  }

  f <- nlmixr(model1)

  expect_true(regexpr(rex::rex("X(0)"), f$rxode.pred) != -1)
  expect_true(regexpr(rex::rex("X(0)"), paste(deparse(f$dynmodel.fun), collapse = "\n")) == -1)
},
test = "cran"
)
