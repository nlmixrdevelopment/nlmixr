nlmixrTest({
  context("mixed pk/ode")
  ## good
  one1 <- function() {
    ini({
      tka <- .5 # Log Ka
      tcl <- -3.2 # Log Cl
      tv <- -1 # Log V
      d0 <- 1
      eta.ka ~ 1
      eta.cl ~ 2
      eta.v ~ 1
      eta.d0 ~ 1
      add.err <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      depot(0) <- d0 + eta.d0
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err)
    })
  }

  one2 <- function() {
    ini({
      tka <- .5 # Log Ka
      tcl <- -3.2 # Log Cl
      tv <- -1 # Log V
      d0 <- 1
      eta.ka ~ 1
      eta.cl ~ 2
      eta.v ~ 1
      eta.d0 ~ 1
      add.err <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      depot(0) <- d0 + eta.d0 + v
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err)
    })
  }

  one3 <- function() {
    ini({
      tka <- .5 # Log Ka
      tcl <- -3.2 # Log Cl
      tv <- -1 # Log V
      d0 <- 1
      eta.ka ~ 1
      eta.cl ~ 2
      eta.v ~ 1
      eta.d0 ~ 1
      add.err <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      if (t < 10) {
        depot(0) <- d0 + eta.d0 + v
      } else {
        depot(0) <- d0 + v
      }
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err)
    })
  }


  one4 <- function() {
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
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      if (t < 10) {
        depot(0) <- v
      } else {
        depot(0) <- v + ka
      }
      cp ~ add(add.err)
    })
  }


  ## Good.
  one5 <- function() {
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
      d / dt(depot) <- -ka * depot
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err)
    })
  }


  one5 <- function() {
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
      d / dt(depot) <- -ka * depot
      if (depot > 0) {
        cl <- exp(tcl + eta.cl)
      } else {
        cl <- 3
      }
      v <- exp(tv + eta.v)
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.err)
    })
  }
}, test="cran")
