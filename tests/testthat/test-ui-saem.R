nlmixrTest({
  context("UI SAEM tests")

  m1 <- function() {
    ini({ # Where initial conditions/variables are specified
      # '<-' or '=' defines population parameters
      lVM <- 7 # log Vmax (mg/hr)
      lKM <- 6 # log KM (mg/L)
      lVc <- 4 # log V (L)
      # Bounds may be specified by c(lower, est, upper), like NONMEM:
      # Residuals errors are assumed to be population parameters
      prop.err <- c(0, 0.2, 1)
      # Between subject variability estimates are specified by '~'
      # Semicolons are optional
      eta.KM ~ 0.3
      eta.Vc ~ 0.1
      eta.VM ~ 0.2
    })
    model({ # Where the model is specified
      # The model uses the ini-defined variable names
      Vc <- exp(lVc + eta.Vc)
      VM <- exp(lVM + eta.VM)
      KM <- exp(lKM + eta.KM)
      # RxODE-style differential equations are supported
      d / dt(centr) <- -(VM * centr / Vc) / (KM + centr / Vc)
      ## Concentration is calculated
      cp <- centr / Vc
      # And is assumed to follow proportional error estimated by prop.err
      cp ~ prop(prop.err)
    })
  }


  test_that("Initial estimate order is correct", {
    tmp <- nlmixr(m1)
    expect_equal(log(tmp$saem.init$theta), c(4, 7, 6))
    expect_equal(tmp$saem.init$omega, c(0.1, 0.2, 0.3))
  })

  m1 <- function() {
    ini({ # Where initial conditions/variables are specified
      # '<-' or '=' defines population parameters
      # Simple numeric expressions are supported
      lCl <- log(200) # log Cl/F (L/hr)
      lVc <- log(60) # log V/F (L)
      lKA <- log(1.5) # log KA (/hr)
      # Bounds may be specified by c(lower, est, upper), like NONMEM:
      # Residuals errors are assumed to be population parameters
      prop.err <- c(0, 0.2, 1)
      # Between subject variability estimates are specified by '~'
      # Semicolons are optional
      eta.Cl ~ 0.1
      eta.Vc ~ 0.1
      eta.KA ~ 0.1
      AllomCL <- 0.75
      AllomV <- 1
    })
    model({ # Where the model is specified
      # The model uses the ini-defined variable names
      Cl <- exp(lCl + eta.Cl + lWT70 * AllomCL)
      Vc <- exp(lVc + eta.Vc + lWT70 * AllomV)
      KA <- exp(lKA + eta.KA)
      # RxODE-style differential equations are supported
      d / dt(depot) <- -KA * depot
      d / dt(centr) <- KA * depot - (Cl / Vc) * centr
      ## Concentration is calculated
      cp <- centr / Vc
      # And is assumed to follow proportional error estimated by prop.err
      cp ~ prop(prop.err)
    })
  }

  test_that("Covariates make sense.", {
    m1 <- nlmixr(m1)
    expect_equal(c(200, 60, 1.5, 0.75, 1, NA), m1$saem.init.theta)
    expect_equal(c("lCl", "AllomCL", "lVc", "AllomV", "lKA"), m1$saem.theta.name)
  })
}, test="saem")
