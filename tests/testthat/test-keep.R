nlmixrTest({

  one.compartment <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
    keep = c("WT")
    drop = c("depot")
  }

  for (est in c("fo", "foi", "foce", "focei", "saem", "nlme", "posthoc")) {
    test_that(paste0("keep/drop in ", est), {
      if (est == "nlme") {
        fitF <- nlmixr(one.compartment, theo_sd, est="nlme", control=nlmeControl(pnlsTol=0.5))
      } else {
        fitF <- nlmixr(one.compartment, theo_sd, est=est)
      }
      expect_true(any(names(fitF) == "WT"))
      expect_true(!any(names(fitF) == "depot"))
      expect_true(!any(names(fitF) == "rxLambda"))
      expect_true(!any(names(fitF) == "rxYj"))
    })
  }

  one.compartment <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
    keep = c("WT")
    drop = c("depot")
  }

  fitF <- nlmixr(one.compartment, theo_sd, est="focei")
  expect_true(any(names(fitF) == "WT"))
  expect_true(!any(names(fitF) == "depot"))
  expect_true(!any(names(fitF) == "rxLambda"))
  expect_true(!any(names(fitF) == "rxYj"))




}, test="lvl2")
