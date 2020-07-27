nlmixrTest({
    context("Bad UI models should raise errors")

    test_that("Duplicate parameters raise errors", {
      uif <- function() {
        ini({
          lCL <- 1.37
          lV <- 4.19
          lCLD <- 1.37
          lVT <- 3.87
          prop.err <- 1
          eta.Cl ~ 0.1
          eta.V ~ 0.1
          ## Duplicate CLs
          eta.Cl ~ 0.1
          eta.VT ~ 0.1
        })
        model({
          CL <- exp(lCL + eta.Cl)
          V <- exp(lV + eta.V)
          CLD <- exp(lCLD + eta.Cl)
          VT <- exp(lVT + eta.VT)
          ## FIXME possibly include both?
          ## K10 <- CL / V
          ## K12 <- CLD / V
          ## K21 <- CLD / VT
          linCmt() ~ prop(prop.err)
        })
      }

      expect_error(nlmixr(uif), rex::rex("duplicated parameter names: 'eta.Cl'"))
    })

    test_that("Un-estimated paramteres raise errors", {
      uif.ode <- function() {
        ini({
          lCL <- 1.37
          lV <- 4.19
          lCLD <- 1.37
          lVT <- 3.87
          ## Prop error isn't estimated
          prop.err <- 1
          add.err <- 0.1
          eta.Cl + eta.V ~ c(
            0.1,
            0.01, 0.01
          )
        })
        model({
          CL <- exp(lCL + eta.Cl)
          V <- exp(lV + eta.V)
          CLD <- exp(lCLD)
          VT <- exp(lVT)
          K10 <- CL / V
          K12 <- CLD / V
          K21 <- CLD / VT
          d / dt(centr) <- K21 * periph - K12 * centr - K10 * centr
          d / dt(periph) <- -K21 * periph + K12 * centr
          cp <- centr / V
          cp ~ add(add.err)
        })
      }

      expect_error(nlmixr(uif.ode), rex::rex("The following parameter(s) were in the ini block but not in the model block: prop.err"))

      uif <- function() {
        ini({
          tka <- exp(0.5)
          tcl <- exp(-3.2)
          tv <- exp(1)
          eta.ka ~ 0.1
          ## Should be eta.cl
          eta.v ~ 0.2
          add.err ~ 0.1
        })
        model({
          ka <- tka + eta.ka
          cl <- tcl + eta.cl
          v <- tv
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp ~ add(add.err)
        })
      }

      expect_error(nlmixr(uif), rex::rex("The following parameter(s) were in the ini block but not in the model block: eta.v"))
    })


    test_that("Residuals are population parameters", {
      uif <- function() {
        ini({
          tka <- exp(0.5)
          tcl <- exp(-3.2)
          tv <- exp(1)
          eta.ka ~ 0.1
          eta.cl ~ 0.2
          ## Should be assign since it is a THETa, should I support it....?
          add.err ~ 0.1
        })
        model({
          ka <- tka + eta.ka
          cl <- tcl + eta.cl
          v <- tv
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp ~ add(add.err)
        })
      }

      expect_error(nlmixr(uif), rex::rex("Residual error component(s) need to be defined with assignment ('=' or '<-') in ini block (not '~'): add.err"))
    })

    test_that("Parameters need to be named", {
      uif <- function() {
        ini({
          tka <- exp(0.5)
          tcl <- exp(-3.2)
          tv <- exp(1)
          eta.ka ~ 0.1
          eta.cl ~ 0.2
          ## Should be assign since it is a THETa, should I support it....?
          0.1
        })
        model({
          ka <- tka + eta.ka
          cl <- tcl + eta.cl
          v <- tv
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp ~ add(add.err)
        })
      }
      ## , rex::rex("The following THETAs are unnamed: THETA[4]")

      expect_error(nlmixr(uif))

      uif <- function() {
        ini({
          tka <- exp(0.5)
          tcl <- exp(-3.2)
          tv <- exp(1)
          eta.ka ~ 0.1
          ~0.2
          ## Should be assign since it is a THETa, should I support it....?
          add.err <- 0.1
        })
        model({
          ka <- tka + eta.ka
          cl <- tcl + eta.cl
          v <- tv
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp ~ add(add.err)
        })
      }

      ## rex::rex("The following ETAs are unnamed: ETA[2]")
      expect_error(nlmixr(uif))
    })

    test_that("Parameters cannot be missing or Infinite", {
      uif <- function() {
        ini({
          tka <- 1 / 0
          tcl <- exp(-3.2)
          tv <- exp(1)
          eta.ka ~ 0.1
          eta.cl ~ 0.2
          add.err <- 4
        })
        model({
          ka <- tka + eta.ka
          cl <- tcl + eta.cl
          v <- tv
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp ~ add(add.err)
        })
      }
      expect_error(nlmixr(uif), rex::rex("The following parameters initial estimates are infinite: tka"))

      uif <- function() {
        ini({
          tka <- NA
          tcl <- exp(-3.2)
          tv <- exp(1)
          eta.ka ~ 0.1
          eta.cl ~ 0.2
          add.err <- 4
        })
        model({
          ka <- tka + eta.ka
          cl <- tcl + eta.cl
          v <- tv
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp ~ add(add.err)
        })
      }

      expect_error(nlmixr(uif), rex::rex("The following parameters initial estimates are NA: tka"))

      uif <- function() {
        ini({
          tka <- 3
          tcl <- exp(-3.2)
          tv <- exp(1)
          eta.ka ~ 0.1
          eta.cl ~ 0.2
          add.err
        })
        model({
          ka <- tka + eta.ka
          cl <- tcl + eta.cl
          v <- tv
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp ~ add(add.err)
        })
      }

      expect_error(nlmixr(uif), rex::rex("Residual distribution parameter(s) estimates were not found in ini block"))
    })

    test_that("Parameters cannot be missing or Infinite", {
      uif <- function() {
        ini({
          tka <- NA
          tcl <- exp(-3.2)
          tv <- exp(1)
          eta.ka ~ 0.1
          eta.cl ~ 0.2
          add.err <- 1
        })
        model({
          ka <- tka + eta.ka
          cl <- tcl + eta.cl
          v <- tv
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          cp <- add(add.err)
        })
      }
      expect_error(nlmixr(uif), "There must be at least one prediction")
    })
  },
 test="cran"
)
