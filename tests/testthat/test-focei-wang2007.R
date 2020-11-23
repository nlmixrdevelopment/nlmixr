nlmixrTest(
  {
    context("Test FO, FOCE, and FOCEi objective functions")

    ## For some reason the ODE and solved FOCE proportional models
    ## give quite different results.  However, FOCE doesn't work as
    ## well with prop models.
    .foceiFit <- function(...) {
      suppressWarnings(foceiFit(...))
    }

    mypar1 <- function() {
      ke <- theta[1] * exp(eta[1])
    }

    mypar2 <- function() {
      k <- theta[1] * exp(eta[1])
      v <- 1
    }

    mod <- RxODE({
      ipre <- 10 * exp(-ke * t)
    })

    out.focei.prop <- readRDS(test_path("out.focei.prop.rds"))

    pred <- function() ipre

    inits <- list(THTA = c(0.5))
    inits$OMGA <- list(ETA[1] ~ .04)

    dat <- Wang2007
    dat$DV <- dat$Y

    m1 <- RxODE({
      d / dt(ipre) <- -ke * ipre
    })

    ## Enhance data frame to include dosing records.
    dat2 <- dat[dat$Time == 0, ]
    dat2$EVID <- 101
    dat2$AMT <- 10
    dat2 <- rbind(dat2, data.frame(dat, EVID = 0, AMT = 0))
    dat2 <- dat2[(order(dat2$ID, -dat2$EVID, dat2$Time)), ]

    testErr <- function(type, fun, val = rep(NA_real_, 6), addProp = 2) {
      fit1 <- .foceiFit(dat2, inits, mypar1, m1, pred, fun,
        control = foceiControl(
          maxOuterIterations = 0, covMethod = "",
          addProp = paste0("combined", addProp)
        )
      )
      fit2 <- .foceiFit(dat, inits, mypar1, mod, pred, fun,
        control = foceiControl(
          maxOuterIterations = 0, covMethod = "",
          addProp = paste0("combined", addProp)
        )
      )
      fit3 <- .foceiFit(dat2, inits, mypar1, m1, pred, fun,
        control = foceiControl(
          maxOuterIterations = 0, covMethod = "",
          addProp = paste0("combined", addProp),
          interaction = FALSE
        )
      )
      fit4 <- .foceiFit(dat, inits, mypar1, mod, pred, fun,
        control = foceiControl(
          maxOuterIterations = 0, covMethod = "",
          addProp = paste0("combined", addProp), interaction = FALSE
        ),
        interaction = FALSE
      )
      fit5 <- .foceiFit(dat2, inits, mypar1, m1, pred, fun,
        control = foceiControl(
          maxOuterIterations = 0, covMethod = "",
          addProp = paste0("combined", addProp),
          interaction = FALSE, fo = TRUE
        )
      )
      fit6 <- .foceiFit(dat, inits, mypar1, mod, pred, fun,
        control = foceiControl(
          maxOuterIterations = 0, covMethod = "",
          addProp = paste0("combined", addProp), interaction = FALSE
        ),
        interaction = FALSE, fo = TRUE
      )
      .n <- paste(type, c("focei ode", "focei", "foce ode", "foce", "fo ode", "fo"), paste0("combined", addProp))
      ret <- c(fit1$objective, fit2$objective, fit3$objective, fit4$objective, fit5$objective, fit6$objective)
      ret <- setNames(ret, .n)
      val <- setNames(val, .n)
      ## Now test
      if (!all(is.na(val))) {
        test_that(
          type,
          expect_equal(round(ret, 3), round(val, 3))
        )
      }
      return(ret)
    }

    fit.prop <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
      return(prop(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "")
    )

    fit.prop2 <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
      return(prop(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "")
    )


    test_that("Matches NONMEM objective proportional function; (Based on Wang2007)", {
      expect_equal(round(fit.prop$objective, 3), 39.458) # Matches Table 2 Prop FOCEI for NONMEM
      expect_equal(round(fit.prop$`ETA[1]`, 4), round(out.focei.prop$ETA1, 4)) # match NONMEM output
      ## Individual properties
      expect_equal(round(fit.prop$IPRED, 4), round(out.focei.prop$IPRE, 4))
      expect_equal(round(fit.prop$IRES, 4), round(out.focei.prop$IRES, 4))
      expect_equal(round(fit.prop$IWRES, 4), round(out.focei.prop$IWRES, 4))
      ## WRES variants
      expect_equal(round(fit.prop$PRED, 4), round(out.focei.prop$NPRED, 4)) # matches output of PRED from NONMEM
      expect_equal(round(fit.prop$PRED, 4), round(out.focei.prop$PRED, 4)) # matches output of PRED from NONMEM
      expect_equal(round(fit.prop$RES, 4), round(out.focei.prop$RES, 4)) # match NONMEM output
      expect_equal(round(fit.prop$RES, 4), round(out.focei.prop$NRES, 4)) # match NONMEM output
      ## FOI equivalents
      expect_equal(round(fit.prop$PRED, 4), round(out.focei.prop$PREDI, 4)) # matches output of PRED from NONMEM
      ## CWRES variants
      expect_equal(round(fit.prop$CRES, 4), round(out.focei.prop$CRES, 4)) # match NONMEM output
      expect_equal(round(fit.prop$CPRED, 4), round(out.focei.prop$CPRED, 4)) # match NONMEM output
      expect_equal(round(fit.prop$CWRES, 4), round(out.focei.prop$CWRES, 4)) # match NONMEM output
      ## Note that E[x] for CPRED and CPREDI are equal
      expect_equal(round(fit.prop$CRES, 4), round(out.focei.prop$CRESI, 4)) # match NONMEM output
      expect_equal(round(fit.prop$CPRED, 4), round(out.focei.prop$CPREDI, 4)) # match NONMEM output
    })
    test_that("Matches NONMEM objective proportional function; (Based on Wang2007; unoptimized)", {
      # Check unoptimized expression
      expect_equal(round(fit.prop2$objective, 3), 39.458) # Matches Table 2 Prop FOCEI for NONMEM
      expect_equal(round(fit.prop2$`ETA[1]`, 4), round(out.focei.prop$ETA1, 4)) # match NONMEM output
      ## Individual properties
      expect_equal(round(fit.prop2$IPRED, 4), round(out.focei.prop$IPRE, 4))
      expect_equal(round(fit.prop2$IRES, 4), round(out.focei.prop$IRES, 4))
      expect_equal(round(fit.prop2$IWRES, 4), round(out.focei.prop$IWRES, 4))
      ## WRES variants
      expect_equal(round(fit.prop2$PRED, 4), round(out.focei.prop$NPRED, 4)) # matches output of PRED from NONMEM
      expect_equal(round(fit.prop2$PRED, 4), round(out.focei.prop$PRED, 4)) # matches output of PRED from NONMEM
      expect_equal(round(fit.prop2$RES, 4), round(out.focei.prop$RES, 4)) # match NONMEM output
      expect_equal(round(fit.prop2$RES, 4), round(out.focei.prop$NRES, 4)) # match NONMEM output
      ## FOI equivalents
      expect_equal(round(fit.prop2$PRED, 4), round(out.focei.prop$PREDI, 4)) # matches output of PRED from NONMEM
      ## CWRES variants
      expect_equal(round(fit.prop2$CRES, 4), round(out.focei.prop$CRES, 4)) # match NONMEM output
      expect_equal(round(fit.prop2$CPRED, 4), round(out.focei.prop$CPRED, 4)) # match NONMEM output
      expect_equal(round(fit.prop2$CWRES, 4), round(out.focei.prop$CWRES, 4)) # match NONMEM output
      ## Note that E[x] for CPRED and CPREDI are equal
      expect_equal(round(fit.prop2$CRES, 4), round(out.focei.prop$CRESI, 4)) # match NONMEM output
      expect_equal(round(fit.prop2$CPRED, 4), round(out.focei.prop$CPREDI, 4)) # match NONMEM output
    })

    testErr("prop", function() {
      prop(.1)
    }, c(39.458, 39.458, 39.275, 39.207, 39.213, 39.207))

    testErr("add", function() {
      add(.1)
    }, c(-2.059, -2.059, -2.059, -2.059, 0.026, -2.059))

    testErr("add+prop", function() {
      add(.1) + prop(.1)
    }, c(39.735, 39.735, 39.562, 39.499, 39.505, 39.499), addProp = 2)

    testErr("add+prop", function() {
      add(.1) + prop(.1)
    }, c(43.554, 43.554, 43.416, 43.394, 43.398, 43.394), addProp = 1)

    testErr("add+pow", function() {
      add(.1) + pow(.1, 0.5)
    }, c(21.108, 21.108, 21.065, 20.725, 20.788, 20.725), addProp = 2)

    testErr("add+pow", function() {
      add(.1) + pow(.1, 0.5)
    }, c(26.065, 26.065, 26.011, 25.873, 25.897, 25.873), addProp = 1)

    testErr("lnorm", function() {
      lnorm(0.1)
    }, c(40.039, 40.039, 40.039, 40.039, 40.055, 40.039), addProp = 1)

    testErr("lnorm(NA)+prop", function() {
      lnorm(NA) + prop(0.1)
    }, c(118.419, 118.419, 118.279, 118.311, 118.311, 118.311), addProp = 1)

    testErr("lnorm(NA)+pow", function() {
      lnorm(NA) + pow(0.1, 0.5)
    }, c(94.535, 94.535, 94.461, 94.478, 94.478, 94.478), addProp = 1)

    testErr("lnorm+prop", function() {
      lnorm(0.1) + prop(0.1)
    }, c(123.318, 123.318, 123.219, 123.24, 123.24, 123.24), addProp = 1)

    testErr("lnorm+prop", function() {
      lnorm(0.1) + prop(0.1)
    }, c(118.777, 118.777, 118.646, 118.676, 118.676, 118.676), addProp = 2)

    testErr("lnorm+pow", function() {
      lnorm(0.1) + pow(0.1, 0.5)
    }, c(102.899, 102.899, 102.855, 102.865, 102.865, 102.865), addProp = 1)

    testErr("lnorm+pow", function() {
      lnorm(0.1) + pow(0.1, 0.5)
    }, c(95.634, 95.634, 95.57, 95.585, 95.585, 95.585), addProp = 2)

    ## Box Cox
    testErr("add+boxCox", function() {
      return(add(.1) + boxCox(.5))
    }, c(2.06, 2.06, 2.06, 2.06, 3.529, 2.06))

    testErr("prop+boxCox", function() {
      return(prop(.1) + boxCox(.5))
    }, c(61.473, 61.473, 61.298, 61.324, 61.325, 61.324))

    testErr("pow+boxCox", function() {
      return(pow(.1, 0.5) + boxCox(.5))
    }, c(39.567, 39.567, 39.477, 39.402, 39.411, 39.402))

    testErr("add+prop+boxCox", function() {
      return(add(0.1) + prop(.1) + boxCox(.5))
    }, c(66.075, 66.075, 65.949, 65.973, 65.974, 65.973), addProp = 1)

    testErr("add+prop+boxCox", function() {
      return(add(0.1) + prop(.1) + boxCox(.5))
    }, c(61.802, 61.802, 61.636, 61.661, 61.662, 61.661), addProp = 2)

    testErr("add+pow+boxCox", function() {
      return(add(0.1) + pow(.1, 0.5) + boxCox(.5))
    }, c(46.768, 46.768, 46.709, 46.693, 46.695, 46.693), addProp = 1)

    testErr("add+pow+boxCox", function() {
      return(add(0.1) + pow(.1, 0.5) + boxCox(.5))
    }, c(40.451, 40.451, 40.372, 40.313, 40.32, 40.313), addProp = 2)

    # Now yeoJohnson
    testErr("add+yeoJohnson", function() {
      return(add(.1) + yeoJohnson(.5))
    }, c(2.339, 2.339, 2.339, 2.339, 3.575, 2.339))

    testErr("prop+yeoJohnson", function() {
      return(prop(.1) + yeoJohnson(.5))
    }, c(62.821, 62.821, 62.647, 62.676, 62.677, 62.676))

    testErr("pow+yeoJohnson", function() {
      return(pow(.1, 0.5) + yeoJohnson(.5))
    }, c(40.724, 40.724, 40.632, 40.575, 40.581, 40.575))

    testErr("add+prop+yeoJohnson", function() {
      return(add(0.1) + prop(.1) + yeoJohnson(.5))
    }, c(67.453, 67.453, 67.329, 67.353, 67.354, 67.353), addProp = 1)

    testErr("add+prop+yeoJohnson", function() {
      return(add(0.1) + prop(.1) + yeoJohnson(.5))
    }, c(63.152, 63.152, 62.988, 63.016, 63.017, 63.016), addProp = 2)

    testErr("add+pow+yeoJohnson", function() {
      return(add(0.1) + pow(.1, 0.5) + yeoJohnson(.5))
    }, c(48.036, 48.036, 47.978, 47.967, 47.969, 47.967), addProp = 1)

    testErr("add+pow+yeoJohnson", function() {
      return(add(0.1) + pow(.1, 0.5) + yeoJohnson(.5))
    }, c(41.628, 41.628, 41.548, 41.503, 41.509, 41.503), addProp = 2)

    ## logitNorm
    testErr("logitNorm", function() {
      return(logitNorm(.1, 0, 12))
    }, c(0.612, 0.612, 0.612, 0.612, 0.786, 0.612))

    testErr("logitNorm(NA)+prop", function() {
      return(logitNorm(NA, 0, 12) + prop(0.1))
    }, c(67.882, 67.882, 67.731, 67.765, 67.765, 67.765))

    testErr("logitNorm(NA)+pow", function() {
      return(logitNorm(NA, 0, 12) + pow(0.1, 0.5))
    }, c(44.632, 44.632, 44.542, 44.556, 44.556, 44.556))

    testErr("logitNorm+prop", function() {
      return(logitNorm(.1, 0, 12) + prop(0.1))
    }, c(72.699, 72.699, 72.591, 72.615, 72.615, 72.615), addProp = 1)

    testErr("logitNorm+prop", function() {
      return(logitNorm(.1, 0, 12) + prop(0.1))
    }, c(68.233, 68.233, 68.09, 68.122, 68.123, 68.122), addProp = 2)

    testErr("logitNorm+pow", function() {
      return(logitNorm(.1, 0, 12) + pow(0.1, 0.5))
    }, c(52.641, 52.641, 52.589, 52.6, 52.6, 52.6), addProp = 1)

    testErr("logitNorm+pow", function() {
      return(logitNorm(.1, 0, 12) + pow(0.1, 0.5))
    }, c(45.668, 45.668, 45.59, 45.603, 45.603, 45.603), addProp = 2)

    ## logitNorm + yeoJohnson
    testErr("logitNorm+yeoJohnson", function() {
      return(logitNorm(.1, 0, 12) + yeoJohnson(0.5))
    }, c(5.127, 5.127, 5.127, 5.127, 5.484, 5.127))

    testErr("logitNorm(NA)+prop+yeoJohnson", function() {
      return(logitNorm(NA, 0, 12) + prop(0.1) + yeoJohnson(0.5))
    }, c(73.881, 73.881, 73.73, 73.764, 73.764, 73.764))

    testErr("logitNorm(NA)+pow+yeoJohnson", function() {
      return(logitNorm(NA, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
    }, c(50.631, 50.631, 50.54, 50.551, 50.552, 50.551))

    testErr("logitNorm+prop+yeoJohnson", function() {
      return(logitNorm(.1, 0, 12) + prop(0.1) + yeoJohnson(0.5))
    }, c(78.693, 78.693, 78.585, 78.609, 78.609, 78.609), addProp = 1)

    testErr("logitNorm+prop+yeoJohnson", function() {
      return(logitNorm(.1, 0, 12) + prop(0.1) + yeoJohnson(0.5))
    }, c(74.231, 74.231, 74.088, 74.12, 74.12, 74.12), addProp = 2)

    testErr("logitNorm+pow+yeoJohnson", function() {
      return(logitNorm(.1, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
    }, c(58.628, 58.628, 58.575, 58.586, 58.586, 58.586), addProp = 1)

    testErr("logitNorm+pow+yeoJohnson", function() {
      return(logitNorm(.1, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
    }, c(51.662, 51.662, 51.584, 51.595, 51.595, 51.595), addProp = 2)

    ## probitNorm
    testErr("probitNorm", function() {
      return(probitNorm(.1, 0, 12))
    }, c(12.827, 12.827, 12.827, 12.827, 12.847, 12.827))

    testErr("probitNorm(NA)+prop", function() {
      return(probitNorm(NA, 0, 12) + prop(0.1))
    }, c(88.875, 88.875, 88.733, 88.766, 88.766, 88.766))

    testErr("probitNorm(NA)+pow", function() {
      return(probitNorm(NA, 0, 12) + pow(0.1, 0.5))
    }, c(65.098, 65.098, 65.02, 65.037, 65.037, 65.037))

    testErr("probitNorm+prop", function() {
      return(probitNorm(0.1, 0, 12) + prop(0.1))
    }, c(93.761, 93.761, 93.661, 93.682, 93.682, 93.682), addProp = 1)

    testErr("probitNorm+prop", function() {
      return(probitNorm(0.1, 0, 12) + prop(0.1))
    }, c(89.232, 89.232, 89.098, 89.129, 89.129, 89.129), addProp = 2)

    testErr("probitNorm+pow", function() {
      return(probitNorm(0.1, 0, 12) + pow(0.1, 0.5))
    }, c(73.405, 73.405, 73.359, 73.37, 73.37, 73.37), addProp = 1)

    testErr("probitNorm+pow", function() {
      return(probitNorm(0.1, 0, 12) + pow(0.1, 0.5))
    }, c(66.187, 66.187, 66.12, 66.135, 66.135, 66.135), addProp = 2)

    ## probitNorm + yeoJohnson

    testErr("probitNorm+yeoJohnson", function() {
      return(probitNorm(.1, 0, 12) + yeoJohnson(0.5))
    }, c(16.768, 16.768, 16.768, 16.768, 16.799, 16.768))

    testErr("probitNorm(NA)+prop+yeoJohnson", function() {
      return(probitNorm(NA, 0, 12) + prop(0.1) + yeoJohnson(0.5))
    }, c(93.071, 93.071, 92.929, 92.962, 92.962, 92.962))

    testErr("probitNorm(NA)+pow+yeoJohnson", function() {
      return(probitNorm(NA, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
    }, c(69.295, 69.295, 69.217, 69.234, 69.234, 69.234))

    testErr("probitNorm(0.1)+prop+yeoJohnson", function() {
      return(probitNorm(0.1, 0, 12) + prop(0.1) + yeoJohnson(0.5))
    }, c(97.957, 97.957, 97.856, 97.878, 97.878, 97.878), addProp = 1)

    testErr("probitNorm(0.1)+prop+yeoJohnson", function() {
      return(probitNorm(0.1, 0, 12) + prop(0.1) + yeoJohnson(0.5))
    }, c(93.429, 93.429, 93.295, 93.326, 93.326, 93.326), addProp = 2)

    testErr("probitNorm(0.1)+pow+yeoJohnson", function() {
      return(probitNorm(0.1, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
    }, c(77.599, 77.599, 77.553, 77.564, 77.564, 77.564), addProp = 1)

    testErr("probitNorm(0.1)+pow+yeoJohnson", function() {
      return(probitNorm(0.1, 0, 12) + pow(0.1, 0.5) + yeoJohnson(0.5))
    }, c(70.383, 70.383, 70.316, 70.331, 70.331, 70.331), addProp = 2)

    ## lognormal -- equivalent to add on log-space and back-transformed.

    ## Next run on the log-transformed space
    datl <- dat
    datl$DV <- log(datl$DV)
    datl2 <- dat2
    datl2$DV <- log(datl2$DV)
    predl <- function() log(ipre)

    fit.lnorm <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
      return(lnorm(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "")
    )

    test_that("Matches NONMEM objective lognormal function; (Based on Wang2007)", {
      expect_equal(round(fit.lnorm$objective, 3), 40.039)
    })

    fit.lnorm0 <- .foceiFit(datl, inits, mypar1, mod, predl, function() {
      return(add(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "")
    )

    test_that("Matches NONMEM objective lognormal function; (Based on Wang2007)", {
      expect_equal(round(fit.lnorm$objective, 3), 40.039)
      expect_equal(fit.lnorm0$objective + 2 * sum(datl$DV), fit.lnorm$objective)
      expect_equal(round(fit.lnorm0$objective, 3), -42.106)
    })

    fit.lnorm2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
      return(lnorm(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "")
    )

    fit.lnorm20 <- .foceiFit(datl2, inits, mypar1, m1, predl, function() {
      return(add(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "")
    )

    test_that("Matches NONMEM objective lognormal function; ODE (Based on Wang2007)", {
      expect_equal(round(fit.lnorm2$objective, 3), 40.039)
      expect_equal(fit.lnorm20$objective + 2 * sum(datl$DV), fit.lnorm2$objective)
      expect_equal(round(fit.lnorm20$objective, 3), -42.106)
    })

    fit.lnorm <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
      return(lnorm(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
    )

    fit.lnorm0 <- .foceiFit(datl, inits, mypar1, mod, predl, function() {
      return(add(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
    )

    test_that("Matches NONMEM objective lognormal error FOCE (Based on Wang2007)", {
      expect_equal(round(fit.lnorm$objective, 3), 40.039)
      expect_equal(fit.lnorm0$objective + 2 * sum(datl$DV), fit.lnorm$objective)
      expect_equal(round(fit.lnorm0$objective, 3), -42.106)
    })

    fit.lnorm2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
      return(lnorm(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
    )

    fit.lnorm20 <- .foceiFit(datl2, inits, mypar1, m1, predl, function() {
      return(add(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
    )

    test_that("Matches NONMEM objective lognormal error FOCE; ODE (Based on Wang2007)", {
      expect_equal(round(fit.lnorm2$objective, 3), 40.039)
      expect_equal(fit.lnorm20$objective + 2 * sum(datl$DV), fit.lnorm2$objective)
      expect_equal(round(fit.lnorm20$objective, 3), -42.106)
    })

    fit.lnorm <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
      return(lnorm(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
    )

    fit.lnorm0 <- .foceiFit(datl, inits, mypar1, mod, predl, function() {
      return(add(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
    )

    test_that("Matches NONMEM objective lognormal error FO (Based on Wang2007)", {
      expect_equal(round(fit.lnorm$objective, 3), 40.055)
      expect_equal(fit.lnorm0$objective + 2 * sum(datl$DV), fit.lnorm$objective)
      expect_equal(round(fit.lnorm0$objective, 3), -42.09)
    })

    fit.lnorm2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
      return(lnorm(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
    )

    fit.lnorm20 <- .foceiFit(datl2, inits, mypar1, m1, predl, function() {
      return(add(.1))
    },
    control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
    )

    test_that("Matches NONMEM objective lognormal error FO; ODE (Based on Wang2007)", {
      expect_equal(round(fit.lnorm2$objective, 3), 40.055)
      expect_equal(fit.lnorm20$objective + 2 * sum(datl$DV), fit.lnorm2$objective)
      expect_equal(round(fit.lnorm20$objective, 3), -42.09)
    })
  },
  test = "focei"
)
