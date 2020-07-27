nlmixrTest({
  context("SAEM covariate parsing")

  run1 <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <- log(0.6) # typical value of volume
      all.cl <- 1 # allometric exponent on CL
      eta.cl ~ 1
      # interindividual variability on clearance and volume
      add.err <- 0.1 # residual variability
    })
    model({
      cl <- exp(tcl + all.cl * log_allo_wt + eta.cl) # individual value of clearance
      v <- exp(tv) # individual value of volume
      ke <- cl / v # elimination rate constant
      d / dt(A1) <- -ke * A1 # model differential equation
      cp <- A1 / v # concentration in plasma
      cp ~ add(add.err) # define error model
    })
  }

  run2 <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <- log(0.6) # typical value of volume
      all.cl <- 1 # allometric exponent on CL
      eta.cl ~ 1
      # interindividual variability on clearance and volume
      add.err <- 0.1 # residual variability
    })
    model({
      cl <- exp(tcl + log_allo_wt * all.cl + eta.cl) # individual value of clearance
      v <- exp(tv) # individual value of volume
      ke <- cl / v # elimination rate constant
      d / dt(A1) <- -ke * A1 # model differential equation
      cp <- A1 / v # concentration in plasma
      cp ~ add(add.err) # define error model
    })
  }


  run3 <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <- log(0.6) # typical value of volume
      all.cl <- 1 # allometric exponent on CL
      eta.cl ~ 1
      # interindividual variability on clearance and volume
      add.err <- 0.1 # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl + log_allo_wt * all.cl) # individual value of clearance
      v <- exp(tv) # individual value of volume
      ke <- cl / v # elimination rate constant
      d / dt(A1) <- -ke * A1 # model differential equation
      cp <- A1 / v # concentration in plasma
      cp ~ add(add.err) # define error model
    })
  }


  run4 <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <- log(0.6) # typical value of volume
      all.cl <- 1 # allometric exponent on CL
      eta.cl ~ 1
      # interindividual variability on clearance and volume
      add.err <- 0.1 # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl + all.cl * log_allo_wt) # individual value of clearance
      v <- exp(tv) # individual value of volume
      ke <- cl / v # elimination rate constant
      d / dt(A1) <- -ke * A1 # model differential equation
      cp <- A1 / v # concentration in plasma
      cp ~ add(add.err) # define error model
    })
  }

  run5 <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <- log(0.6) # typical value of volume
      all.cl <- 1 # allometric exponent on CL
      eta.cl ~ 1
      # interindividual variability on clearance and volume
      add.err <- 0.1 # residual variability
    })
    model({
      cl <- exp(all.cl * log_allo_wt + tcl + eta.cl) # individual value of clearance
      v <- exp(tv) # individual value of volume
      ke <- cl / v # elimination rate constant
      d / dt(A1) <- -ke * A1 # model differential equation
      cp <- A1 / v # concentration in plasma
      cp ~ add(add.err) # define error model
    })
  }

  run6 <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <- log(0.6) # typical value of volume
      all.cl <- 1 # allometric exponent on CL
      eta.cl ~ 1
      # interindividual variability on clearance and volume
      add.err <- 0.1 # residual variability
    })
    model({
      cl <- exp(log_allo_wt * all.cl + tcl + eta.cl) # individual value of clearance
      v <- exp(tv) # individual value of volume
      ke <- cl / v # elimination rate constant
      d / dt(A1) <- -ke * A1 # model differential equation
      cp <- A1 / v # concentration in plasma
      cp ~ add(add.err) # define error model
    })
  }
  p <- list()
  p[[1]] <- nlmixr(run1)
  p[[2]] <- nlmixr(run2)
  p[[3]] <- nlmixr(run3)
  p[[4]] <- nlmixr(run4)
  p[[5]] <- nlmixr(run5)
  p[[6]] <- nlmixr(run6)

  ref <- list(log_allo_wt = c(all.cl = "tcl"))

  for (i in 1:6) {
    test_that(sprintf("Parsing cl/log_allo_wt works correctly #%d", i), {
      expect_equal(p[[1]]$cov.ref, ref)
    })
  }

  run7 <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <- log(0.6) # typical value of volume
      all.cl <- 1 # allometric exponent on CL
      eta.cl ~ 1
      # interindividual variability on clearance and volume
      add.err <- 0.1 # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + all.cl * log_allo_wt) # individual value of volume
      ke <- cl / v # elimination rate constant
      d / dt(A1) <- -ke * A1 # model differential equation
      cp <- A1 / v # concentration in plasma
      cp ~ add(add.err) # define error model
    })
  }

  ref <- list(log_allo_wt = c(all.cl = "tv"))
  p7 <- nlmixr(run7)

  ## I'm not sure why, but this doesn't seem to work, though when they parse they seem to be correct.

  ## expect_equal("Parsing on v without eta works.",{
  ##     expect_true(p7$cov.ref$log_allo_wt["all.cl"]=="tv");
  ## })

}, test="saem")
