nlmixrTest({
    context("Make sure SAEM/nlme throws an error with time varying covariates")
    test_that("Error w/time-varying covariates", {
      d <- theo_sd
      ## Make this time-varying
      d$WT[d$TIME > 12] <- d$WT[d$TIME > 12] + 0.01
      one.cmt <- function() {
        ini({
          tka <- 0.45 # Log Ka
          tcl <- 1 # Log Cl
          tv <- 3.45 # Log V
          all.cl <- 0
          eta.ka ~ 0.6
          eta.cl ~ 0.3
          eta.v ~ 0.1
          add.err <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl + all.cl * wt)
          v <- exp(tv + eta.v)
          linCmt() ~ add(add.err)
        })
      }

      expect_error({
        f <- nlmixr(one.cmt, d, "nlme")
      })

      expect_error({
        f <- nlmixr(one.cmt, d, "saem")
      })

      f <- suppressWarnings(nlmixr(one.cmt, d, "focei", control = list(sigdig = 2, print = 0)))

      expect_true(inherits(f, "nlmixrFOCEi"))
    })
  },
  test = "cran"
)
