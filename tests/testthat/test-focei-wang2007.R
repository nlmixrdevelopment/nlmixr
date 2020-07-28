nlmixrTest({
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

  fit.prop <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  fit.prop2 <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", optExpression = FALSE)
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

  m1 <- RxODE({
    d / dt(ipre) <- -ke * ipre
  })

  ## Enhance data frame to include dosing records.
  dat2 <- dat[dat$Time == 0, ]
  dat2$EVID <- 101
  dat2$AMT <- 10
  dat2 <- rbind(dat2, data.frame(dat, EVID = 0, AMT = 0))
  dat2 <- dat2[(order(dat2$ID, -dat2$EVID, dat2$Time)), ]

  fit.prop2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )


  test_that("Matches NONMEM objective proportional function; ODE (Based on Wang2007)", {
    expect_equal(round(fit.prop2$objective, 3), 39.458)
  })


  fit.prop <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective proportional error FOCE (Based on Wang2007)", {
    expect_equal(round(fit.prop$objective, 3), 39.207)
  })

  etaMat <- as.matrix(ranef(fit.prop)[, -1, drop = FALSE])

  fit.prop2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective proportional error FOCE; ODE (Based on Wang2007)", {
    expect_equal(round(fit.prop2$objective, 3), 39.279)
  })


  fit.prop <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective proportional error FO (Based on Wang2007)", {
    expect_equal(round(fit.prop$objective, 3), 39.213)
  })

  fit.prop2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective proportional error FO; ODE (Based on Wang2007)", {
    expect_equal(round(fit.prop2$objective, 3), 39.213) # Matches Table 2 Prop FOCEI for NONMEM
  })


  #### Add

  fit.add <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(add(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  test_that("Matches NONMEM objective additive error; (Based on Wang2007)", {
    expect_equal(round(fit.add$objective, 3), -2.059) # Matches Table 2 Add FOCEI for NONMEM
  })

  fit.add2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(add(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  test_that("Matches NONMEM objective additive error; ODE (Based on Wang2007)", {
    expect_equal(round(fit.add2$objective, 3), -2.059)
  })

  fit.add <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(add(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective additive error FOCE (Based on Wang2007)", {
    expect_equal(round(fit.add$objective, 3), -2.059)
  })


  fit.add2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(add(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective additive error FOCE; ODE (Based on Wang2007)", {
    expect_equal(round(fit.add2$objective, 3), -2.059)
  })


  fit.add <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(add(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective additive error FO (Based on Wang2007)", {
    expect_equal(round(fit.add$objective, 3), 0.026)
  })

  fit.add2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(add(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective additive error FO; ODE (Based on Wang2007)", {
    expect_equal(round(fit.add2$objective, 3), 0.026)
  })

  ### Add+Prop

  fit.addprop <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(add(.1) + prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  test_that("Matches NONMEM objective additive+proportional function; (Based on Wang2007)", {
    expect_equal(round(fit.addprop$objective, 3), 39.735)
  })

  fit.addprop2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(add(.1) + prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  test_that("Matches NONMEM objective additive+proportional function; ODE (Based on Wang2007)", {
    expect_equal(round(fit.addprop2$objective, 3), 39.735)
  })

  fit.addprop <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(add(.1) + prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective additive+proportional error FOCE (Based on Wang2007)", {
    expect_equal(round(fit.addprop$objective, 3), 39.499)
  })

  fit.addprop2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(add(.1) + prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective additive+proportional error FOCE; ODE (Based on Wang2007)", {
    expect_equal(round(fit.addprop2$objective, 3), 39.563)
  })

  fit.addprop <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(add(.1) + prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective additive+proportional error FO (Based on Wang2007)", {
    expect_equal(round(fit.addprop$objective, 3), 39.505)
  })

  fit.addprop2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(add(.1) + prop(.1))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective additive+proportional error FO; ODE (Based on Wang2007)", {
    expect_equal(round(fit.addprop2$objective, 3), 39.505)
  })

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

  ### Now try  TBS boxCox

  fit.boxCox <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1) + boxCox(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  test_that("Matches NONMEM objective Box-Cox function; (Based on Wang2007)", {
    expect_equal(round(fit.boxCox$objective, 3), 61.473)
  })

  fit.boxCox2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1) + boxCox(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  test_that("Matches NONMEM objective Box-Cox function; ODE (Based on Wang2007)", {
    expect_equal(round(fit.boxCox2$objective, 3), 61.473)
  })

  fit.boxCox <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1) + boxCox(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective Box-Cox error FOCE (Based on Wang2007)", {
    expect_equal(round(fit.boxCox$objective, 3), 61.324)
  })

  fit.boxCox2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1) + boxCox(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective Box-Cox error FOCE; ODE (Based on Wang2007)", {
    expect_equal(round(fit.boxCox2$objective, 3), 61.28)
  })

  fit.boxCox <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1) + boxCox(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective Box-Cox error FO (Based on Wang2007)", {
    expect_equal(round(fit.boxCox$objective, 3), 61.325)
  })

  fit.boxCox2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1) + boxCox(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective Box-Cox error FO; ODE (Based on Wang2007)", {
    expect_equal(round(fit.boxCox2$objective, 3), 61.325)
  })

  ### Now try  TBS

  fit.yeoJohnson <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1) + yeoJohnson(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  test_that("Matches NONMEM objective Yeo-Johnson function; (Based on Wang2007)", {
    expect_equal(round(fit.yeoJohnson$objective, 3), 62.821)
  })

  fit.yeoJohnson2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1) + yeoJohnson(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "")
  )

  test_that("Matches NONMEM objective Yeo-Johnson function; ODE (Based on Wang2007)", {
    expect_equal(round(fit.yeoJohnson2$objective, 3), 62.821)
  })

  fit.yeoJohnson <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1) + yeoJohnson(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective Yeo-Johnson error FOCE (Based on Wang2007)", {
    expect_equal(round(fit.yeoJohnson$objective, 3), 62.676)
  })

  fit.yeoJohnson2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1) + yeoJohnson(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE)
  )

  test_that("Matches NONMEM objective Yeo-Johnson error FOCE; ODE (Based on Wang2007)", {
    expect_equal(round(fit.yeoJohnson2$objective, 3), 62.627)
  })

  fit.yeoJohnson <- .foceiFit(dat, inits, mypar1, mod, pred, function() {
    return(prop(.1) + yeoJohnson(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective Yeo-Johnson error FO (Based on Wang2007)", {
    expect_equal(round(fit.yeoJohnson$objective, 3), 62.677)
  })

  fit.yeoJohnson2 <- .foceiFit(dat2, inits, mypar1, m1, pred, function() {
    return(prop(.1) + yeoJohnson(.5))
  },
  control = foceiControl(maxOuterIterations = 0, covMethod = "", interaction = FALSE, fo = TRUE)
  )

  test_that("Matches NONMEM objective Yeo-Johnson error FO; ODE (Based on Wang2007)", {
    expect_equal(round(fit.yeoJohnson2$objective, 3), 62.677)
  })
}, test="focei")
